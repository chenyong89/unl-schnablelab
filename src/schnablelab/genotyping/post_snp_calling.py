"""
Post SNP calling functions
"""
import sys
import subprocess
import numpy as np
import pandas as pd
import os.path as op

from tqdm import tqdm
from pathlib import Path
from subprocess import run
from .base import ParseVCF
from schnablelab.apps.base import (ActionDispatcher,
                                   OptionParser,
                                   glob,
                                   create_slurm)

# the location of executable linkimpute, beagle
exe_dir = op.abspath(op.dirname(__file__))
lkipt = exe_dir + '/../apps/LinkImpute.jar'
begle = exe_dir + '/../apps/beagle.24Aug19.3e8.jar'
tassel = exe_dir + '/../apps/tassel-5-standalone/run_pipeline.pl'


def main():
    actions = (
        ('qc_missing', 'Remove SNPs with high missing rate'),
        ('qc_maf', 'Remove SNPs with rare MAF'),
        ('qc_het', 'Remove SNPs with high heterizygous rate'),
        ('index_vcf', 'index vcf using bgzip and tabix'),
        ('split_vcf', 'split a vcf to smaller files with equal size'),
        ('merge_vcf', 'merge split vcf files follow specified pattern to one'),
        ('cat_fq', 'combine fq files follow specified pattern into one file'),
        ('impute_beagle', 'impute missing data in vcf using beagle'),
        ('fix_indel', 'Fix the InDels problems in hmp file from Tassel'),
        ('FilterVCF', 'remove bad snps using bcftools'),
        ('qc_ALT', 'filter number of ALTs using bcftools'),
        ('fix_GT_sep', 'fix separator in vcf from freebayes'),
        ('cal_LD', 'calculate LD r2 using Plink'),
        ('summarize_LD', 'summarize LD results from Plink')
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def qc_missing(args):
    """
    %prog filter_missing input_vcf

    Remove SNPs with high missing rate
    """
    p = OptionParser(qc_missing.__doc__)
    p.add_option('--cutoff', default=0.7, type='float',
                 help='SNPs with missing rate > cutoff will be removed')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputvcf, = args
    outputvcf = Path(inputvcf).name.replace(
        '.vcf', '_mis%s.vcf' % opts.missing_cutoff)

    vcf = ParseVCF(inputvcf)
    n = 0
    with open(outputvcf, 'w') as f:
        f.writelines(vcf.HashChunk)
        pbar = tqdm(vcf.missing_rate,
                    total=vcf.num_SNPs,
                    desc='Filter Missing',
                    position=0)
        for i, miss in pbar:
            if miss <= opts.missing_cutoff:
                f.write(i)
            else:
                n += 1
            pbar.set_description('processing chromosome %s' % i.split()[0])
    print('Done! %s SNPs removed! check output %s...' % (n, outputvcf))


def qc_maf(args):
    """
    %prog qc_maf input_vcf

    Remove rare MAF SNPs
    """
    p = OptionParser(qc_maf.__doc__)
    p.add_option('--maf_cutoff', default=0.01, type='float',
                 help='SNPs with MAF lower than this cutoff will be removed')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputvcf, = args
    outputvcf = Path(inputvcf).name.replace(
        '.vcf', '_maf%s.vcf' % opts.maf_cutoff)

    vcf = ParseVCF(inputvcf)
    n = 0
    with open(outputvcf, 'w') as f:
        f.writelines(vcf.HashChunk)
        pbar = tqdm(vcf.maf, total=vcf.num_SNPs, desc='Filter MAF', position=0)
        for i, maf in pbar:
            if maf >= opts.maf_cutoff:
                f.write(i)
            else:
                n += 1
            pbar.set_description('processing chromosome %s' % i.split()[0])
    print('Done! %s SNPs removed! check output %s...' % (n, outputvcf))


def qc_het(args):
    """
    %prog qc_het input_vcf

    Remove SNPs with high heterizygous rate
    """
    p = OptionParser(qc_het.__doc__)
    p.add_option('--het_cutoff', default=0.1, type='float',
                 help='SNPs with heterozygous rate > cutoff will be removed')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputvcf, = args
    outputvcf = Path(inputvcf).name.replace(
        '.vcf', '_het%s.vcf' % opts.het_cutoff)

    vcf = ParseVCF(inputvcf)
    n = 0
    with open(outputvcf, 'w') as f:
        f.writelines(vcf.HashChunk)
        pbar = tqdm(vcf.hetero_rate, total=vcf.num_SNPs,
                    desc='Filter Heterozygous', position=0)
        for i, het in pbar:
            if het <= opts.het_cutoff:
                f.write(i)
            else:
                n += 1
            pbar.set_description('processing chromosome %s' % i.split()[0])
    print('Done! %s SNPs removed! check output %s...' % (n, outputvcf))


def qc_ALT(args):
    """
    %prog in_dir out_dir

    filter number of ALTs using bcftools
    """
    p = OptionParser(qc_ALT.__doc__)
    p.add_option('--pattern', default='*.vcf',
                 help='file pattern for vcf files in dir_in')
    p.add_slurm_opts(job_prefix=qc_ALT.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')
    dir_path = Path(in_dir)
    vcfs = dir_path.glob(opts.pattern)
    for vcffile in vcfs:
        prefix = '.'.join(vcffile.name.split('.')[0:-1])
        new_f = prefix + '.alt1.vcf'
        cmd = "bcftools view -i 'N_ALT=1' %s > %s" % (vcffile, new_f)
        cmd_header = 'ml bacftools\n'
        slurm_dict = vars(opts)
        slurm_dict['job_prefix'] = slurm_dict['job_prefix']+prefix
        slurm_dict['cmd_header'] = cmd_header
        create_slurm([cmd], slurm_dict)


def fix_GT_sep(args):
    """
    %prog fix_GT_sep in_dir out_dir

    replace the allele separator `.` in freebayes vcf file to `/` which is
    required for beagle
    """
    p = OptionParser(fix_GT_sep.__doc__)
    p.add_option('--pattern', default='*.vcf',
                 help='file pattern for vcf files in dir_in')
    p.add_slurm_opts(job_prefix=fix_GT_sep.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')
    dir_path = Path(in_dir)
    vcfs = dir_path.glob(opts.pattern)
    for vcf in vcfs:
        sm = '.'.join(vcf.name.split('.')[0:-1])
        out_fn = sm+'.fixGT.vcf'
        out_fn_path = out_path/out_fn
        cmd = "perl -pe 's/\s\.:/\t.\/.:/g' %s > %s" % (vcf, out_fn_path)
        slurm_dict = vars(opts)
        slurm_dict['job_prefix'] = slurm_dict['job_prefix']+sm
        create_slurm([cmd], slurm_dict)


def index_vcf(args):
    """
    %prog index_vcf in_dir out_dir

    index vcf using bgzip and tabix
    """
    p = OptionParser(index_vcf.__doc__)
    p.add_option('--pattern', default='*.vcf',
                 help='file pattern for vcf files in dir_in')
    p.add_slurm_opts(job_prefix=index_vcf.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')
    dir_path = Path(in_dir)
    vcfs = dir_path.glob(opts.pattern)
    for vcf in vcfs:
        sm = '.'.join(vcf.name.split('.')[0:-1])
        out_fn = vcf.name+'.gz'
        out_fn_path = out_path/out_fn
        cmd1 = 'bgzip -c %s > %s\n' % (vcf, out_fn_path)
        cmd2 = 'tabix -p vcf %s\n' % (out_fn_path)
        cmd_header = 'ml tabix\n'
        slurm_dict = vars(opts)
        slurm_dict['job_prefix'] = slurm_dict['job_prefix']+sm
        slurm_dict['cmd_header'] = cmd_header
        create_slurm([cmd1, cmd2], slurm_dict)


def split_vcf(args):
    """
    %prog split_vcf N vcf

    split vcf to N smaller files with equal size
    """
    p = OptionParser(split_vcf.__doc__)
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    N, vcffile, = args
    N = int(N)
    prefix = vcffile.split('.')[0]
    cmd_header = "sed -ne '/^#/p' %s > %s.header" % (vcffile, prefix)
    subprocess.call(cmd_header, shell=True)
    child = subprocess.Popen('wc -l %s' % vcffile, shell=True,
                             stdout=subprocess.PIPE)
    total_line = int(child.communicate()[0].split()[0])
    print('total %s lines' % total_line)
    step = total_line / N
    print(1)
    cmd_first = "sed -n '1,%sp' %s > %s.1.vcf" % (step, vcffile, prefix)
    subprocess.call(cmd_first, shell=True)
    for i in range(2, N):
        print(i)
        st = (i - 1) * step + 1
        ed = i * step
        cmd = "sed -n '%s,%sp' %s > %s.%s.tmp.vcf" % (
            st, ed, vcffile, prefix, i)
        subprocess.call(cmd, shell=True)
    print(i + 1)
    cmd_last = "sed -n '%s,%sp' %s > %s.%s.tmp.vcf" % (
        (ed + 1), total_line, vcffile, prefix, (i + 1))
    subprocess.call(cmd_last, shell=True)
    for i in range(2, N + 1):
        cmd_cat = 'cat %s.header %s.%s.tmp.vcf > %s.%s.vcf' % (
            prefix, prefix, i, prefix, i)
        subprocess.call(cmd_cat, shell=True)


def cat_fq(args):
    """
    %prog cat_fq pattern(with quotation) fn_out

    combine fq files follow specified pattern into one file
    """
    p = OptionParser(cat_fq.__doc__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    fq_pattern, fn_out, = args
    fns = glob(fq_pattern)
    cmd = 'cat %s > %s' % (' '.join(fns), fn_out)
    print(cmd)
    run(cmd, shell=True)


def merge_vcf(args):
    """
    %prog merge_vcf pattern out_fn

    merge split vcf files follow specified pattern to one
    Pattern example:
        'hmp321_agpv4_chr9.%s.beagle.vcf'
    """

    p = OptionParser(merge_vcf.__doc__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    pattern, out_fn, = args

    fns = [str(i) for i in list(Path('.').glob(pattern))]
    fns_sorted = sorted(fns, key=lambda x: int(x.split('.')[0][3:]))
    print(fns_sorted)
    print('%s files found!' % len(fns_sorted))

    f = open(out_fn, 'w')
    print(fns_sorted[0])
    with open(fns_sorted[0]) as f1:
        for i in f1:
            f.write(i)
    for i in fns_sorted[1:]:
        print(i)
        with open(i) as f2:
            for j in f2:
                if not j.startswith('#'):
                    f.write(j)


def impute_beagle(args):
    """
    %prog impute_beagle dir_in dir_out

    impute missing data in vcf using beagle
    """
    p = OptionParser(impute_beagle.__doc__)
    p.add_option('--pattern', default='*.vcf',
                 help='file pattern for vcf files in dir_in')
    p.add_option('--parameter_file',
                 help='file including window, overlap parameters')
    p.add_slurm_opts(job_prefix=impute_beagle.__name__)
    p.set_cpus(cpus=10)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')

    if opts.parameter_file:
        df = pd.read_csv(opts.parameter_file)
        df = df.set_index('chr')

    dir_path = Path(in_dir)
    vcfs = dir_path.glob(opts.pattern)
    for vcf in vcfs:
        print(vcf.name)
        chrom = int(vcf.name.split('.')[0].split('hr')[-1])
        print('chr : %s' % chrom)
        sm = '.'.join(vcf.name.split('.')[0:-1])
        out_fn = sm+'.BG'
        out_fn_path = out_path/out_fn
        if opts.parameter_file:
            window = df.loc[chrom, 'marker_10cm']
            overlap = df.loc[chrom, 'marker_2cm']
            print('window: %s; overlap: %s' % (window, overlap))
            cmd = (f'beagle -Xmx60G gt={vcf} out={out_fn_path} '
                   f'window={window} overlap={overlap} nthreads=10')
        else:
            cmd = f'beagle -Xmx60G gt={vcf} out={out_fn_path} nthreads=10'
        cmd_header = 'ml beagle/4.1\n'
        slurm_dict = vars(opts)
        slurm_dict['job_prefix'] = slurm_dict['job_prefix']+sm
        slurm_dict['cmd_header'] = cmd_header
        create_slurm([cmd], slurm_dict)


def fix_indel(args):
    """
    %prog fix_index hmp

    Fix the InDels problems in hmp file generated from Tassel
    """
    p = OptionParser(fix_indel.__doc__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    hmpfile, = args
    prefix = '.'.join(hmpfile.split('.')[0:-1])

    f = open(hmpfile)
    f1 = open('%s.Findel.hmp' % prefix, 'w')

    bict = {'A': 'T', 'T': 'C', 'C': 'A', 'G': 'A'}
    for i in f:
        j = i.split()
        if '+' in j[1]:
            nb = bict[j[1][0]] if j[1][0] != '+' else bict[j[1][-1]]
            tail = '\t'.join(j[11:]) + '\n'
            n_tail = tail.replace('+', nb)
            head = '\t'.join(j[0:11])
            n_head = head.replace('/+', '/%s' % nb) \
                if j[1][0] != '+' \
                else head.replace('+/', '%s/' % nb)
            n_line = n_head + '\t' + n_tail
            f1.write(n_line)
        elif '-' in j[1]:
            nb = bict[j[1][0]] if j[1][0] != '-' else bict[j[1][-1]]
            tail = '\t'.join(j[11:]) + '\n'
            n_tail = tail.replace('-', nb)
            head = '\t'.join(j[0:11])
            n_head = head.replace('/-', '/%s' % nb) \
                if j[1][0] != '-' \
                else head.replace('-/', '%s/' % nb)
            n_line = n_head + '\t' + n_tail
            f1.write(n_line)
        else:
            f1.write(i)
    f.close()
    f1.close()


def cal_LD(args):
    """
    %prog cal_LD vcf_fn/plink_prefix genome_size(Mb) num_SNPs

    calculate LD using Plink
    args:
        vcf_fn/plink_prefix:
            specify of vcf/vcf.gz file or plink bed/bim/fam files
        genome_size(Mb):
            the size of the reference genome in Mb. use 684 for sorghum
        num_SNPs:
            the number of SNPs in the genotype file.
    """
    p = OptionParser(cal_LD.__doc__)
    p.add_option('--maf_cutoff', default='0.01',
                 help='only use SNP with the MAF > cutoff for LD calculation')
    p.add_option('--max_distance', type='int', default=1000000,
                 help='max distance (bp) of a pair of SNPs to calcualte LD')
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='do not convert commands to slurm jobs')
    p.add_slurm_opts(job_prefix=cal_LD.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_fn, g_size, n_snps, = args
    in_fn, g_size, n_snps = Path(in_fn), int(g_size)*1000000, int(n_snps)

    if in_fn.name.endswith('.vcf') or in_fn.name.endswith('.vcf.gz'):
        input = f'--vcf {in_fn}'
    else:
        input = f'--bfile {in_fn}'
    n = 10
    ld_window, ld_window_bp = [], []
    while True:
        ld_window.append(n)
        dist = g_size//n_snps*n
        ld_window_bp.append(dist)
        n *= 10
        if dist >= 1000000:
            break

    out_fn = Path(in_fn).name.split('.')[0]
    cmds = []
    cmd = (f'plink {input} --r2 --ld-window 10 '
           f'--ld-window-kb {ld_window_bp[0]//1000} --ld-window-r2 0 '
           f'--maf {opts.maf_cutoff} --out {out_fn}')
    cmds.append(cmd)
    for win_snp, win_bp in zip(ld_window[1:], ld_window_bp[1:]):
        prob = 10/win_snp
        cmd = (f'plink {input} --thin {prob} --r2 --ld-window 10 '
               f'--ld-window-kb {win_bp//1000} --ld-window-r2 0 '
               f'--maf {opts.maf_cutoff} --out {out_fn}.thin{prob}')
        cmds.append(cmd)
        print(cmd)
    cmd_sh = '%s.cmds%s.sh' % (opts.job_prefix, len(cmds))
    pd.DataFrame(cmds).to_csv(cmd_sh, index=False, header=None)
    print(f'check {cmd_sh} for all the commands!')

    if not opts.disable_slurm:
        cmd_header = 'ml plink'
        slurm_dict = vars(opts)
        slurm_dict['cmd_header'] = cmd_header
        create_slurm(cmds, slurm_dict)


def summarize_LD(args):
    """
    %prog summarizeLD ld1 ld2 ... output_prefix

    summarize LD results from Plink
    args:
        ld*:
            the output files from Plink
        output_prefix:
            the output prefix of summarized results
    """
    p = OptionParser(summarize_LD.__doc__)
    p.add_option('--order', type='int', default=4,
                 help='degree of the fitting polynomial used in numpy.polyfit')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    *lds, out_prefix = args
    df = pd.concat([pd.read_csv(ld, delim_whitespace=True, usecols=[1, 4, 6])
                    for ld in lds])
    df = df.drop_duplicates(['BP_A', 'BP_B'])
    df['Dist_bp'] = df['BP_B']-df['BP_A']
    log_bin = [10**i for i in np.arange(0.1, float(6)+0.1, 0.1)]
    inds = np.digitize(df['Dist_bp'].values, bins=log_bin)
    df['x'] = inds
    df_plot = df.groupby('x')['R2'].mean().reset_index()
    df_plot.to_csv(f'{out_prefix}.LDsummary.csv', index=False, sep='\t')

    import seaborn as sns
    import matplotlib.pyplot as plt
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    plt.rcParams['xtick.major.pad'] = 1
    plt.rcParams['ytick.major.pad'] = 1
    _, ax = plt.subplots(figsize=(5.5, 4))
    ax = sns.regplot(x='x', y='R2', data=df_plot, order=opts.order,
                     scatter=True,
                     scatter_kws={'edgecolor': 'k', 'facecolor': 'white',
                                  's': 20})

    line1 = plt.Line2D(range(1), range(1), linestyle='none', color="k",
                       marker='o',
                       markerfacecolor="white",
                       markersize=5)
    line2 = plt.Line2D(range(1), range(1), color="#1f77b4")
    plt.legend((line1, line2), ('Original', 'Fitted'), frameon=False,
               numpoints=2)

    ax.set_xlim(0, 60)
    ax.set_ylim(bottom=0)

    ax.set_xlabel('Distance', fontsize=12)
    ax.set_ylabel(r'$r^2$', fontsize=12)

    ax.set_xticklabels(['1bp', '10bp', '100bp', '1Kb', '10Kb', '100Kb', '1Mb'])

    ax.spines['bottom'].set_position(('axes', -0.01))
    ax.spines['left'].set_position(('axes', -0.01))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.savefig(f'{out_prefix}.png', dpi=300)


if __name__ == "__main__":
    main()
