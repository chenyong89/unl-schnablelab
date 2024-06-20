"""
Call SNPs on high throughput sequencing data using GATK, Freebayes, SAMtools
"""

import os
import re
import sys
import pandas as pd
from pathlib import Path

from .base import find_sm
from schnablelab.apps.tools import create_df_from_path
from schnablelab.apps.base import ActionDispatcher, OptionParser, create_slurm


def main():
    actions = (
        ('gen_gvcf',
         'generate gvcf for each sample using GATK HaplotypeCaller'),
        ('agg_gvcf',
         'aggregate gvcf files to a GenomicsDB for each genomic interval'),
        ('gen_vcf', 'create the raw vcf from GenomicsDB datastores'),
        ('freebayes', 'call SNPs using freebayes'),
        ('gatk', 'call SNPs using gatk'),)
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def gen_gvcf(args):
    """
    %prog gen_gvcf ref.fa bams.csv region.txt out_dir

    run GATK HaplotypeCaller in GVCF mode
    one g.vcf file for one smaple may contain multiple replicates
    args:
        ref.fa: reference sequence file
        bams.csv: csv file containing all bam files and their sample names
        region.txt:
            genomic intervals defined by each row to speed up GVCF calling
            e.g.: Chr01, Chr01:1-100
        out_dir: where the gVCF files save to
    """
    p = OptionParser(gen_gvcf.__doc__)
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='do not convert commands to slurm jobs')
    p.add_slurm_opts(job_prefix=gen_gvcf.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    ref, bams_csv, region_txt, out_dir, = args
    out_dir_path = Path(out_dir)
    if not out_dir_path.exists():
        print(f'output directory {out_dir_path} does not exist, creating...')
        out_dir_path.mkdir()
    regions = []
    with open(region_txt) as f:
        for i in f:
            regions.append(i.rstrip())
    mem = int(opts.memory)//1024
    df_bam = pd.read_csv(bams_csv)
    # check if bai files exist
    for bam in df_bam['fnpath']:
        if not Path(bam+'.bai').exists():
            print(f'no index file for {bam}...')
            sys.exit('Index your bam files first!')

    cmds = []
    for sm, grp in df_bam.groupby('sm'):
        print(f'{grp.shape[0]} bam files for sample {sm}')
        input_bam = '-I ' + ' -I '.join(grp['fnpath'].tolist())
        for region in regions:
            output_fn = f'{sm}_{region}.g.vcf'
            cmd = (f"gatk --java-options '-Xmx{mem}g' "
                   f"HaplotypeCaller -R {ref} "
                   f"{input_bam} -O {out_dir_path/output_fn} "
                   f"--sample-name {sm} "
                   f"--emit-ref-confidence GVCF -L {region}")
            cmds.append(cmd)
    cmd_sh = f'{opts.job_prefix}.cmds{len(cmds)}.sh'
    pd.DataFrame(cmds).to_csv(cmd_sh, index=False, header=None)
    print(f'check {cmd_sh} for all the commands!')

    cmd_header = 'ml gatk4/4.1'
    if not opts.disable_slurm:
        slurm_dict = vars(opts)
        slurm_dict['cmd_header'] = cmd_header
        create_slurm(cmds, slurm_dict)


def agg_gvcf(args):
    """
    %prog agg_gvcf input_dir out_dir

    aggregate GVCF files in the input_dir to a GenomicsDB datastore for each
    genomic interval
    args:
        intput_dir:
            the directory containing all gvcf files
        out_dir:
            the output directory where subdir will be created for each
            genomic interval
    """
    p = OptionParser(agg_gvcf.__doc__)
    p.add_option('--gvcf_fn_pattern', default='*.g.vcf',
                 help='file extension of gvcf files')
    p.add_option('--sm_re_pattern',
                 default=r"^P[0-9]{3}[_-]W[A-Z][0-9]{2}[^a-z0-9]",
                 help='RE pattern to pull sample name from filename')
    p.add_option('--gatk_tmp_dir', default='./gatk_tmp',
                 help='temporary directory for genomicsDBImport')
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='do not convert commands to slurm jobs')
    p.add_slurm_opts(job_prefix=agg_gvcf.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args
    in_dir_path = Path(in_dir)
    out_dir_path = Path(out_dir)
    if not in_dir_path.exists():
        sys.exit(f'input directory {in_dir_path} does not exist!')
    if not out_dir_path.exists():
        print(f'output directory {out_dir_path} does not exist, creating...')
        out_dir_path.mkdir()
    tmp_dir = Path(opts.gatk_tmp_dir)
    if not tmp_dir.exists():
        print('tmp directory does not exist, creating...')
        tmp_dir.mkdir()

    # The -Xmx value the tool is run with should be less than the total amount
    # of physical memory available by at least a few GB
    mem = int(opts.memory)//1024-2

    # set the environment variable TILEDB_DISABLE_FILE_LOCKING=1
    try:
        os.environ['TILEDB_DISABLE_FILE_LOCKING']
    except KeyError:
        sys.exit('Set the environment variable TILEDB_DISABLE_FILE_LOCKING=1 '
                 'before running gatk!')

    df = create_df_from_path(in_dir_path, pattern=opts.gvcf_fn_pattern)
    df['interval'] = df['fn'].apply(lambda x: x.split('.')[0].split('_')[1])
    prog = re.compile(opts.sm_re_pattern)
    df['sm'] = df['fn'].apply(lambda x: find_sm(x, prog))

    cmds = []
    for interval, grp in df.groupby('interval'):
        interval_dir = out_dir_path/(interval.replace(':', '_'))
        # The --genomicsdb-workspace-path must point to a non-existent or
        # empty directory
        if interval_dir.exists():
            if len(interval_dir.glob('*')) != 0:
                sys.exit(f'{interval_dir} is not an empty directory!')
        gvcf_map = str(interval) + '.map'
        print(f'{grp.shape[0]} gvcf files found for interval {interval}, '
              f'generating the corresponding map file {gvcf_map}...')
        grp[['sm', 'fnpath']].to_csv(gvcf_map, header=None, index=False,
                                     sep='\t')

        cmd = (f"gatk --java-options '-Xmx{mem}g -Xms{mem}g' GenomicsDBImport "
               f"--sample-name-map {gvcf_map} "
               f"--genomicsdb-workspace-path {interval_dir} "
               f"--batch-size 50 --intervals {interval} "
               f"--reader-threads {opts.ncpus_per_node} --tmp-dir {tmp_dir}")
        cmds.append(cmd)

    cmd_sh = f'{opts.job_prefix}.cmds{len(cmds)}.sh'
    pd.DataFrame(cmds).to_csv(cmd_sh, index=False, header=None)
    print(f'check {cmd_sh} for all the commands!')

    cmd_header = 'ml gatk4/4.1'
    if not opts.disable_slurm:
        slurm_dict = vars(opts)
        slurm_dict['cmd_header'] = cmd_header
        create_slurm(cmds, slurm_dict)


def gen_vcf(args):
    """
    %prog gen_vcf ref.fa genomicDB_dir out_dir

    create the raw VCFs from GenomicsDB datastores
    args:
        ref.fa: the reference sequence fasta file
        genomicDB_dir: the root directory of genomicDB workspace
        out_dir: where the vcf files will be saved
    """
    p = OptionParser(gen_vcf.__doc__)
    p.add_option('--gatk_tmp_dir', default='./gatk_tmp',
                 help='temporary directory to use')
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='do not convert commands to slurm jobs')
    p.add_slurm_opts(job_prefix=gen_vcf.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    ref, db_dir, out_dir, = args
    out_dir_path = Path(out_dir)
    if not out_dir_path.exists():
        print(f'output directory {out_dir_path} does not exist, creating...')
        out_dir_path.mkdir()
    mem = int(opts.memory)//1024-1

    cmds = []
    for db in Path(db_dir).glob('*'):
        if db.is_dir():
            region = db.name
            vcf_fn = f"{region}.vcf.gz"
            cmd = (f"gatk --java-options '-Xmx{mem}g' GenotypeGVCFs "
                   f"-R {ref} -V gendb://{db} -O {out_dir_path/vcf_fn} "
                   f"--tmp-dir={opts.gatk_tmp_dir}")
            cmds.append(cmd)
    cmd_sh = f'{opts.job_prefix}.cmds{len(cmds)}.sh'
    pd.DataFrame(cmds).to_csv(cmd_sh, index=False, header=None)
    print(f'check {cmd_sh} for all the commands!')

    cmd_header = 'ml gatk4/4.1'
    if not opts.disable_slurm:
        slurm_dict = vars(opts)
        slurm_dict['cmd_header'] = cmd_header
        create_slurm(cmds, slurm_dict)


def gatk(args):
    """
    %prog gatk ref.fa bam_list.txt region.txt out_dir

    run GATK HaplotypeCaller
    """
    p = OptionParser(gatk.__doc__)
    p.add_slurm_opts(job_prefix=gatk.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    ref, bams, regions, out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')
    with open(bams) as f:
        inputs = ''.join(['-I %s \\\n' % (i.rstrip()) for i in f])
    with open(regions) as f:
        for reg in f:
            reg = reg.strip()
            if ':0-' in reg:
                reg = reg.replace(':0-', ':1-')
            reg_fn = reg.replace(':', '_')
            reg_fn_vcf = '%s.gatk.vcf' % reg_fn
            reg_fn_vcf_path = out_path/reg_fn_vcf
            cmd = (f"gatk --java-options '-Xmx13G' HaplotypeCaller \\\n"
                   f"-R {ref} -L {reg} \\\n{inputs}-O {reg_fn_vcf_path}")
            cmd_header = 'ml gatk4/4.1\n'
            slurm_dict = vars(opts)
            slurm_dict['job_prefix'] = slurm_dict['job_prefix']+f'_{reg}'
            slurm_dict['cmd_header'] = cmd_header
            create_slurm([cmd], slurm_dict)


def freebayes(args):
    """
    %prog freebayes region.txt ref.fa bam_list.txt out_dir

    create freebayes slurm jobs for each split region defined in region.txt
    file
    """
    p = OptionParser(freebayes.__doc__)
    p.add_option('--max_depth', default=10000,
                 help=('loci where the mapping depth higher than this value '
                       'will be ignored'))
    p.add_slurm_opts(job_prefix=freebayes.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    region, ref, bams, out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')

    with open(region) as f:
        for reg in f:
            reg = reg.strip()
            reg_fn = reg.replace(':', '_')
            reg_fn_vcf = f'{reg_fn}.fb.vcf'
            reg_fn_vcf_path = out_path/reg_fn_vcf
            cmd = (f'freebayes -r {reg} -f {ref} -C 1 -F 0.05 -L {bams} '
                   f'-u -n 2 -g {opts.max_depth} > {reg_fn_vcf_path}')
            cmd_header = 'ml freebayes/1.3\n'
            slurm_dict = vars(opts)
            slurm_dict['job_prefix'] = slurm_dict['job_prefix']+f'_{reg}'
            slurm_dict['cmd_header'] = cmd_header
            create_slurm([cmd], slurm_dict)


if __name__ == "__main__":
    main()
