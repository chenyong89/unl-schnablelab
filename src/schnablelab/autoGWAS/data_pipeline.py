import sys
import pandas as pd
import os.path as op
from tqdm import tqdm
from pathlib import Path
from subprocess import call
from schnablelab.apps.base import ActionDispatcher, OptionParser, create_slurm
from .base import ParseHapmap, geno_to_value


def main():
    actions = (
        ('qc_missing', 'filter out SNPs with high missing rate'),
        ('qc_MAF', 'filter out SNP with extremely low MAF'),
        ('qc_hetero', 'filter out SNPs with high heterozygous rates'),
        ('hapmap_to_map_ped', 'convert hmp file to plink map and ped file'),
        ('ped_to_bed', 'convert plink ped format to binary bed format'),
        ('hapmap_to_vcf', 'transform hapmap format to vcf format'),
        ('hapmap_to_bimbam',
         'transform hapmap format to BIMBAM format for GEMMA'),
        ('hapmap_to_numeric',
         'transform hapmap format to numeric format for GAPIT and FarmCPU'),
        ('hapmap_to_mvp', 'transform hapmap format to MVP genotypic format'),
        ('generate_kinship',
         'using gemma to generate centered kinship matrix'),
        ('generate_pca', 'using tassel to generate the first N PCs'),
        ('independent_SNPs', 'estimate the number of independent SNPs'),
        ('single_code_to_double_code',
         'convert single hmp to double type hmp'),
        ('data_info', 'get basic info for a hmp file'),
        ('cal_MAFs', 'calculate the MAF for all/specified SNPs in hmp'),
        ('sort_hapmap',
         'sort hmp file based on chromosome order and position'),
        ('combine_hapmap',
         'combine split chromosome Hmps to a single large one'),
        ('modify_sample_name', 'modify sample names in table header'),
        ('extract_SNPs', 'extract a subset of specified SNPs from a hmp file'),
        ('extract_samples', 'extract a subset of samples from a hmp file'),
        ('sampling_SNPs', 'randomly sample SNPs from a large hapmap file'),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


plink = op.abspath(op.dirname(__file__)) + '/../apps/plink'
gec = op.abspath(op.dirname(__file__)) + '/../apps/gec.jar'
gemma = op.abspath(op.dirname(__file__)) + '/../apps/gemma'
tassel = op.abspath(op.dirname(__file__)) + \
    '/../apps/tassel-5-standalone/run_pipeline.pl'


def modify_sample_name(args):
    """
    %prog modify_sample_name input_hmp names.csv

    substitute the sample names in hmp header
    name.csv:
        comma separated without header line
        1st column is old name
        2nd column is the new name
    """
    p = OptionParser(modify_sample_name.__doc__)
    _, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, names_csv, = args
    outputhmp = Path(inputhmp).name.replace('.hmp', '_reheader.hmp')

    hmp = ParseHapmap(inputhmp)

    cmd = 'sed '
    for _, row in pd.read_csv(names_csv, header=None).iterrows():
        old_nm, new_nm = row[0], row[1]
        if old_nm not in hmp.SMs:
            print(f'{id} was not found in hmp...')
        else:
            cmd += f"-e '1s/{old_nm}/{new_nm}/' "
    cmd += f'{inputhmp} > {outputhmp}'
    print(f'command:\n{cmd}')
    choice = input("Run the above command? (yes/no) ")
    if choice == 'yes':
        call(cmd, shell=True)
        print(f'Done! check output file: {outputhmp}')


def qc_missing(args):
    """
    %prog qc_missing input_hmp
    Remove SNPs with high missing rate (1st step in QC pipeline)
    """
    p = OptionParser(qc_missing.__doc__)
    p.add_option('--missing_cutoff', default=0.7, type='float', 
                 help='remove SNPs higher than specified missing rate cutoff')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, = args
    outputhmp = Path(inputhmp).name.replace('.hmp', f'_mis{opts.missing_cutoff}.hmp')

    hmp = ParseHapmap(inputhmp)
    n = 0
    with open(outputhmp, 'w') as f:
        f.write(hmp.headerline)
        pbar = tqdm(hmp.Missings, total=hmp.numSNPs)
        for i, miss in pbar:
            if miss <= opts.missing_cutoff:
                f.write(i)
            else:
                n += 1
            pbar.set_description(f'processing chromosome {i.split()[2]}')
    print(f'Done! {n} SNPs removed! check output file: {outputhmp}')


def qc_MAF(args):
    """
    %prog qc_MAF input_hmp
    Remove rare MAF SNPs (2nd step in QC pipeline)
    """
    p = OptionParser(qc_MAF.__doc__)
    p.add_option('--MAF_cutoff', default=0.01, type='float',
                 help='SNPs lower than specified cutoff rate will be removed.')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, = args
    outputhmp = Path(inputhmp).name.replace('.hmp', f'_maf{opts.MAF_cutoff}.hmp')

    hmp = ParseHapmap(inputhmp)
    n = 0
    with open(outputhmp, 'w') as f:
        f.write(hmp.headerline)
        pbar = tqdm(hmp.MAFs, total=hmp.numSNPs)
        for i, maf in pbar:
            if maf >= opts.MAF_cutoff:
                f.write(i)
            else:
                n += 1
            pbar.set_description(f'processing chromosome {i.split()[2]}')
    print(f'Done! {n} SNPs removed! check output {outputhmp}...')


def qc_hetero(args):
    """
    %prog qc_hetero input_hmp

    Remove high heterozygous rate SNPs (3rd step in QC pipeline)
    """
    p = OptionParser(qc_hetero.__doc__)
    p.add_option('--het_cutoff', default=0.1, type='float',
                 help='SNPs higher than specified heterozygous rate will be removed')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, = args
    outputhmp = Path(inputhmp).name.replace('.hmp', f'_het{opts.het_cutoff}.hmp')

    hmp = ParseHapmap(inputhmp)
    n = 0
    with open(outputhmp, 'w') as f:
        f.write(hmp.headerline)
        pbar = tqdm(hmp.Heteros, total=hmp.numSNPs)
        for i, het in pbar:
            if het <= opts.het_cutoff:
                f.write(i)
            else:
                n += 1
            pbar.set_description(f'processing chromosome {i.split()[2]}')
    print(f'Done! {n} SNPs removed! check output {outputhmp}...')


def extract_SNPs(args):
    """
    %prog extract_SNPs input_hmp SNPs.csv

    extract a subset of SNPs defined in `SNPs.csv` from the input_hmp
    put one snp id per row without header in the SNPs.csv file
    """
    p = OptionParser(extract_SNPs.__doc__)
    _, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, SNPcsv, = args
    outputhmp = Path(inputhmp).name.replace('.hmp', '_subSNPs.hmp')

    hmp = ParseHapmap(inputhmp)
    df_hmp = hmp.AsDataframe()

    IDs = pd.read_csv(SNPcsv, header=None)[0].values
    num_IDs = IDs.shape[0]
    print(f'number of specified SNPs: {num_IDs}')
    df_hmp = df_hmp[df_hmp['rs#'].isin(IDs)]
    print(f'{df_hmp.shape[0]} out of {num_IDs} found in Hmp')
    df_hmp.to_csv(outputhmp, sep='\t', index=False, na_rep='NA')
    print('Done! check output {outputhmp}...')


def extract_samples(args):
    """
    %prog extract_samples input_hmp SMs.csv

    extract a subset of samples defined in SMs.csv from the input_hmp
    put one sample name per row without header in SMs.csv
    """
    p = OptionParser(extract_samples.__doc__)
    _, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, SMcsv, = args
    outputhmp = Path(inputhmp).name.replace('.hmp', '_subSMs.hmp')

    hmp = ParseHapmap(inputhmp)
    df_hmp = hmp.to_df()

    IDs = pd.read_csv(SMcsv, header=None)[0].values
    num_IDs = IDs.shape[0]
    print(f'number of specified Samples: {num_IDs}')

    subsm = hmp.attributes
    for id in IDs:
        if id not in hmp.SMs:
            print(f'{id} was not found in hmp...')
        else:
            subsm.append(id)
    print(f'{len(subsm)-11} out of {num_IDs} found in Hmp')

    df_hmp = df_hmp[subsm]
    df_hmp.to_csv(outputhmp, sep='\t', index=False, na_rep='NA')
    print('Done! check output {outputhmp}...')


def sampling_SNPs(args):
    """
    %prog sampling_SNPs input_hmp

    efficiently pick up some SNPs from a huge hmp file
    """
    p = OptionParser(sampling_SNPs.__doc__)
    p.add_option('--downscale', default=10,
                 help='specify the downscale level')
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='do not convert commands to slurm jobs')
    p.add_slurm_opts(job_prefix=sampling_SNPs.__name__)
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    inputhmp, = args
    outputhmp = Path(inputhmp).name.replace('.hmp', f'_ds{opts.downsize}.hmp')
    cmd = "sed -n '1~%sp' %s > %s" % (opts.downsize, inputhmp, outputhmp)
    print('cmd:\n%s\n' % cmd)
    if not opts.disable_slurm:
        slurm_dict = vars(opts)
        create_slurm([cmd], slurm_dict)


def hapmap_to_map_ped(args):
    """
    %prog hapmap_to_map_ped input_hmp

    Convert hmp file to Plink map and ped files
    """
    p = OptionParser(hapmap_to_map_ped.__doc__)
    _, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, = args
    output_prefix = Path(inputhmp).name.split('.hmp')[0]

    hmp = ParseHapmap(inputhmp)
    df_map, df_ped = hmp.to_map_ped(missing=False)
    print('saving map file...')
    df_map.to_csv(f'{output_prefix}.map', sep='\t', index=False, header=None)
    print('saving ped file...')
    df_ped.to_csv('{output_prefix}.ped', sep='\t', index=False, header=None)


def ped_to_bed(args):
    """
    %prog ped_to_bed ped_prefix

    Convert plink ped/map to binary bed/bim/fam format using Plink
    """
    p = OptionParser(ped_to_bed.__doc__)
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='disable converting commands to slurm jobs')
    p.add_slurm_opts(job_prefix=ped_to_bed.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    ped_prefix, = args
    cmd_local = f'{plink} --noweb --file {ped_prefix} --make-bed --out {ped_prefix}'
    print(f'cmd on local desktop:\n{cmd_local}\n')
    if not opts.disable_slurm:
        cmd_header = 'ml plink'
        cmd = f'plink --noweb --file {ped_prefix} --make-bed --out {ped_prefix}'
        slurm_dict = vars(opts)
        slurm_dict['cmd_header'] = cmd_header
        create_slurm([cmd], slurm_dict)


def hapmap_to_vcf(args):
    """
    %prog hapmap_to_vcf input_hmp

    convert hmp to vcf format using tassel
    """
    p = OptionParser(hapmap_to_vcf.__doc__)
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='disable converting commands to slurm jobs')
    p.add_slurm_opts(job_prefix=hapmap_to_vcf.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    hmpfile, = args
    cmd_local = f'{tassel} -Xms512m -Xmx10G -fork1 -h {hmpfile} -export -exportType VCF\n'
    print('cmd to run on local:\n', cmd_local)
    if not opts.disable_slurm:
        cmd_header = 'ml tassel/5.2\n'
        cmd = f'run_pipeline.pl -Xms512m -Xmx10G -fork1 -h {hmpfile} -export -exportType VCF\n'
        slurm_dict = vars(opts)
        slurm_dict['cmd_header'] = cmd_header
        create_slurm([cmd], slurm_dict)


def hapmap_to_bimbam(args):
    """
    %prog hapmap_to_bimbam hmp_file bimbam_file_prefix

    Convert hmp genotypic data to GEMMA bimbam files (*.mean and *.annotation).
    """
    p = OptionParser(hapmap_to_bimbam.__doc__)
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    hmp, bim_prefix = args
    f1 = open(hmp)
    f1.readline()
    f2 = open(bim_prefix + '.mean', 'w')
    f3 = open(bim_prefix + '.annotation', 'w')
    for i in f1:
        j = i.split()
        rs = j[0]
        try:
            ref, alt = j[1].split('/')
        except:
            print('omit rs...')
            continue
        newNUMs = geno_to_value(ref, alt, j[11:])
        newline = '%s,%s,%s,%s\n' % (rs, ref, alt, ','.join(newNUMs))
        f2.write(newline)
        pos = j[3]
        chro = j[2]
        f3.write('%s,%s,%s\n' % (rs, pos, chro))
    f1.close()
    f2.close()
    f3.close()


def hapmap_to_numeric(args):
    """
    %prog hapmap_to_numeric hmp numeric_file_prefix

    Convert hapmap genotypic data to numeric format for GAPIT and FarmCPU 
    output two files: prefix.GD and prefix.GM
    """
    p = OptionParser(hapmap_to_numeric.__doc__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())

    hmp, num_pre = args
    f1 = open(hmp)
    header = f1.readline()
    SMs = '\t'.join(header.split()[11:]) + '\n'
    f2 = open(num_pre + '.GD', 'w')
    f2.write(SMs)
    f3 = open(num_pre + '.GM', 'w')
    f3.write('SNP\tChromosome\tPosition\n')
    for i in f1:
        j = i.split()
        rs = j[0]
        try:
            ref, alt = j[1].split('/')
        except:
            print('omit rs...')
            continue
        newNUMs = geno_to_value(ref, alt, j[11:])
        newline = '\t'.join(newNUMs) + '\n'
        f2.write(newline)
        pos = j[3]
        chro = j[2]
        f3.write('%s\t%s\t%s\n' % (rs, chro, pos))
    f1.close()
    f2.close()
    f3.close()


def hapmap_to_mvp(args):
    """
    %prog hapmap_to_mvp hmp output_prefix

    Convert hmp genotypic data to numeric format for MVP
    output two fiels: prefix.numeric and prefix.map).
    """
    p = OptionParser(hapmap_to_mvp.__doc__)
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    hmp, mvp_pre = args
    f1 = open(hmp)
    f1.readline()
    f2 = open(mvp_pre + '.numeric', 'w')
    f3 = open(mvp_pre + '.map', 'w')
    f3.write('SNP\tChrom\tBP\n')
    for i in f1:
        j = i.split()
        rs = j[0]
        ref, alt = j[1].split('/')[0], j[1].split('/')[1]
        newNUMs = geno_to_value(ref, alt, j[11:])
        newline = '\t'.join(newNUMs) + '\n'
        f2.write(newline)
        chro, pos = j[2], j[3]
        f3.write('%s\t%s\t%s\n' % (rs, chro, pos))
    f1.close()
    f2.close()
    f3.close()


def generate_kinship(args):
    """
    %prog generate_kinship genotype.mean

    Calculate kinship matrix file using gemma
    """
    p = OptionParser(generate_kinship.__doc__)
    p.add_option('--type', default='1', choices=('1', '2'),
                 help='specify the way to calculate the relateness'\
                      '1: centered; 2: standardized')
    p.add_option('--output_dir', default='.',
                 help='specify the output dir')
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='disable converting commands to slurm jobs')
    p.add_slurm_opts(job_prefix=generate_kinship.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    geno_mean, = args
    # generate a fake bimbam phenotype based on genotype file
    with open(geno_mean) as f:
        num_SMs = len(f.readline().split(',')[3:])
    mean_prefix = geno_mean.replace('.mean', '')
    tmp_pheno = '%s.tmp.pheno' % mean_prefix
    with open(tmp_pheno, 'w') as f1:
        for i in range(num_SMs):
            f1.write('sm%s\t%s\n' % (i, 20))
    cmd = '%s -g %s -p %s -gk %s -outdir %s -o gemma.centered.%s' \
        % (gemma, geno_mean, tmp_pheno, opts.type, opts.out_dir,
           Path(mean_prefix).name)
    print('command to run on local:\n', cmd)
    if not opts.disable_slurm:
        slurm_dict = vars(opts)
        create_slurm([cmd], slurm_dict)


def generate_pca(args):
    """
    %prog generate_pca input_hmp N

    Generate first N PCs using tassel
    """
    p = OptionParser(generate_pca.__doc__)
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='disable converting commands to slurm jobs')
    p.add_slurm_opts(job_prefix=generate_pca.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    hmpfile, N, = args
    out_prefix = Path(hmpfile).name.replace('.hmp', '')
    cmd_local = f'{tassel} -Xms28g -Xmx29g -fork1 -h {hmpfile} '\
                f'-PrincipalComponentsPlugin -ncomponents {N} '\
                f'-covariance true -endPlugin -export {out_prefix}_{N}PCA '\
                f'-runfork1\n'
    print('cmd to run on local:\n%s\n%s' % (cmd_local))

    if not opts.disable_slurm:
        cmd_header = 'ml java/1.8\nml tassel/5.2'
        cmd_slurm = f'run_pipeline.pl -Xms28g -Xmx29g -fork1 -h {hmpfile} '\
                    f'-PrincipalComponentsPlugin -ncomponents {N} '\
                    '-covariance true -endPlugin -export '\
                    f'{out_prefix}_{N}PCA -runfork1\n'
        slurm_dict = vars(opts)
        slurm_dict['memory'] = 30000
        slurm_dict['cmd_header'] = cmd_header
        create_slurm([cmd_slurm], slurm_dict)


def independent_SNPs(args):
    """
    %prog independent_SNPs bed_prefix output_fn

    Estimate number of idenpendent SNPs using GEC tool
    """
    p = OptionParser(independent_SNPs.__doc__)
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='add this option to disable converting commands to slurm jobs')
    p.add_slurm_opts(job_prefix=independent_SNPs.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    bed_prefix, output_fn = args
    cmd = 'java -Xmx18g -jar %s --noweb --effect-number --plink-binary %s --genome --out %s' % (gec, bed_prefix, output_fn)
    print('cmd:\n%s\n' % cmd)

    if not opts.disable_slurm:
        slurm_dict = vars(opts)
        slurm_dict['memory'] = 20000
        create_slurm([cmd], slurm_dict)


def single_code_to_double_code(args):
    """
    %prog single_code_to_double_code input_single_hmp 
    convert single type hmp file to double type hmp file
    """
    p = OptionParser(single_code_to_double_code.__doc__)
    _, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, = args
    outputhmp = Path(inputhmp).name.replace('.hmp', '_db.hmp')

    hmp = ParseHapmap(inputhmp)
    df_hmp = hmp.to_df()
    df_hmp.to_csv(outputhmp, sep='\t', index=False, na_rep='NA')
    print(f'Done! check output file: {outputhmp}')


def data_info(args):
    """
    %prog data_info input_hmp

    get basic info for a hmp file
    """
    p = OptionParser(data_info.__doc__)
    _, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, = args
    hmp = ParseHapmap(inputhmp)

    print(f'Genotype type: {hmp.code_type}')
    print('Number of samples: {val:,}'.format(val=hmp.num_sample))
    print('Number of SNPs: {val:,}'.format(val=hmp.num_snp))
    samples = '\n  '.join(hmp.samples)
    print(f'Sample names: \n  {samples}')


def cal_MAFs(args):
    """
    %prog cal_MAFs input_hmp

    calculate MAF for each SNP in hapmap file
    """
    p = OptionParser(cal_MAFs.__doc__)
    _, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, = args
    outputcsv = Path(inputhmp).name.replace('.hmp', '.maf.csv')
    hmp = ParseHapmap(inputhmp)
    with open(outputcsv, 'w') as f:
        pbar = tqdm(hmp.MAFs, total=hmp.numSNPs, desc='get MAF', position=0)
        for i, maf in pbar:
            f.write(f'{maf}\n')
            pbar.set_description('calculating chromosome {i.split()[2]}')
    print('Done! check output file: {outputcsv}')


def sort_hapmap(args):
    """
    %prog sort_hapmap input_hmp
    Sort hmp based on chromosome order and position using python pandas. 
    Can also use tassel command:
     'run_pipeline.pl -Xms16g -Xmx18g -SortGenotypeFilePlugin -inputFile in_fn -outputFile out_fn -fileType Hapmap'
    """
    p = OptionParser(sort_hapmap.__doc__)
    _, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, = args
    outputhmp = Path(inputhmp).name.replace('.hmp', '_sorted.hmp')

    hmp = ParseHapmap(inputhmp)
    df_sorted_hmp = hmp.to_df(sorting=True)
    df_sorted_hmp.to_csv(outputhmp, sep='\t', index=False, na_rep='NA')
    print(f'Done! check output file: {outputhmp}')


def combine_hapmap(args):
    """
    %prog combine_hapmap N pattern output

    combine multiple hapmap files into single one

    args:
    N: number of hapmap files
    pattern example: hmp321_agpv4_chr%s.hmp
    """
    p = OptionParser(combine_hapmap.__doc__)
    p.add_option('--header', default='yes', choices=('yes', 'no'),
                 help='choose whether add header or not')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    N, hmp_pattern, new_f, = args
    N = int(N)

    f = open(new_f, 'w')

    fn1 = open(hmp_pattern % 1)
    print('adding file 1...')
    if opts.header == 'yes':
        for i in fn1:
            f.write(i)
    else:
        fn1.readline()
        for i in fn1:
            f.write(i)
    fn1.close()
    for i in range(2, N + 1):
        print(f'adding file {i}...')
        fn = open(hmp_pattern % i)
        fn.readline()
        for j in fn:
            f.write(j)
        fn.close()
    f.close()


def reorg_TasselPCA_for_gwas(fn_tassel_pca, output_prefix):
    '''
    reorg PCA from tassel for GWAS tools
    '''
    df = pd.read_csv(fn_tassel_pca, delim_whitespace=True, header=2)
    df1 = df[df.columns[1:]]
    gapit_pca = output_prefix + '.gapit'
    farm_pca = output_prefix + '.farmcpu'
    gemma_pca = output_prefix + '.gemmaMVP'
    df.to_csv(gapit_pca, sep='\t', index=False)
    df1.to_csv(farm_pca, sep='\t', index=False)
    df1.to_csv(gemma_pca, sep='\t', index=False, header=False)
    print('%s, %s, %s have been generated.' % (gapit_pca, farm_pca, gemma_pca))


def reorg_GemmaKinship_for_gapit(gemmaKinship, hmp_file):
    '''
    Reorganize kinship result from GEMMA for GAPIT.

    The hmp file only provides the order of the smaple names.
    '''
    f = open(hmp_file)
    SMs = f.readline().split()[11:]
    f.close()
    f1 = open(gemmaKinship)
    f2 = open('GAPIT.' + gemmaKinship, 'w')
    for i, j in zip(SMs, f1):
        newline = i + '\t' + j
        f2.write(newline)
    f1.close()
    f2.close()
    print("'GAPIT.%s' has been generated." % gemmaKinship)


def reorg_pheno_for_gemma(in_dir, out_dir, file_pattern,
                          header=True, sep='\t'):
    """
    reorg phenotype table for GEMMA where missing value will be changed to NA

    args:
        file_pattern (str): pattern of the normal phenotype files. e.g: '*.csv'
        header (bool): whether a header exist in your normal phenotype file
        sep (str): separator in the input phenotype files
    """
    out_path = Path(out_dir)
    if not out_path.exists():
        raise Exception('output dir does not exist...')
    dir_path = Path(in_dir)
    input_files = dir_path.glob(file_pattern)
    for fn in input_files:
        df = pd.read_csv(fn, sep=sep) if header \
            else pd.read_csv(fn, sep=sep, header=None)
        output_fn = fn.name+'.gemma'
        df.iloc[:, 1].to_csv(out_path/output_fn, index=False, header=False,
                             na_rep='NA')
        print('Finished! %s has been generated.' % output_fn)


def explore_residuals(phenotype_file):
    """
    estimate the residual phenotypes from two origianl phenotypes

    args:
        phenotype_file: input phenotype with first three columns:
            name
            pheno1
            pheno2
    """
    from scipy.stats import linregress
    import matplotlib.pyplot as plt
    df = pd.read_csv(phenotype_file)
    pheno1, pheno2 = df.iloc[:, 1], df.iloc[:, 2]

    # plotting
    fig, ax = plt.subplots()
    ax.scatter(pheno1, pheno2, color='lightblue', s=50, alpha=0.8,
               edgecolors='0.3', linewidths=1)
    slope, intercept, r_value, p_value, std_err = linregress(pheno1, pheno2)
    y_pred = intercept + slope * pheno1
    ax.plot(pheno1, y_pred, 'red', linewidth=1, label='Fitted line')
    text_x = max(pheno1) * 0.8
    text_y = max(y_pred)
    ax.text(text_x, text_y, r'${\mathrm{r^2}}$' + ': %.2f' % r_value**2,
            fontsize=15, color='red')
    xlabel, ylabel = df.columns[1], df.columns[2]
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.savefig('%s_%s_r2.png' % (xlabel, ylabel))

    # find residuals
    df['y_1'] = y_pred
    df['residuals'] = df[df.columns[2]] - df['y_1']
    residual_pheno = df[[df.columns[0], df.columns[-1]]]
    residual_pheno.to_csv('residual.csv', sep='\t', na_rep='NaN', index=False)


if __name__ == "__main__":
    main()
