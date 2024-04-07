import sys
import numpy as np
import pandas as pd
import os.path as op
import matplotlib.pyplot as plt
from tqdm import tqdm
from pathlib import Path
from subprocess import call

from .base import ParseGWASresults, cal_MAF_row
from schnablelab.apps.base import (ActionDispatcher,
                                   OptionParser,
                                   create_slurm,
                                   SLURM_header)

plink = op.abspath(op.dirname(__file__))+'/../apps/plink'
faOneRecord = op.abspath(op.dirname(__file__))+'/../apps/faOneRecord'


def main():
    actions = (
        ('plot_manhattan', 'make manhattan plot'),
        ('count_peaks', 'count unique peaks identified in a GWAS run'),
        ('cal_MAF', 'calculate the MAFs of specified significant SNPs'),
        ('fetch_EVs_from_FarmCPU',
         'extract effect sizes of selected SNPs from FarmCPU results'),
        ('shared_sig_SNPs',
         'find shared significant SNPs between GEMMA and FarmCPU'),
        ('sig_SNPs',
         'extract significant SNPs from GWAS results'),
        ('linked_SNPs', 'extract genetically linked SNPs using plink'),
        ('extract_geno_from_vcf', 'fetch genotypes for SNPs from vcf file'),
        ('extract_gene_from_gff3',
         'extract gene names from gff3 for specified SNP list ')
        ('extract_gene_func',
         'extract functions of candidate genes'),
        ('extract_protein_seq',
         'extract protein sequences of condidated genes'),
        ('plot_effect_size', 'plot histgram of effect sizes'),
        ('plot_MAF', 'make density plot of MAF'),
        ('pdf_to_png', 'convert pdf image to png format using imagemagick'),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def plot_manhattan(args):
    """
    %prog plot_manhattan GWAS_result Figure_title

    make Mamhattan plot
    """
    p = OptionParser(plot_manhattan.__doc__)
    p.add_option('--software', default='gemma',
                 choices=('gemma', 'other', 'farmcpu', 'gapit'),
                 help='softare where the GWAS result came from')
    p.add_option('--pvalue', default=0.05, choices=(0.05, 0.01),
                 help='the pvalue cutoff')
    p.add_option('--MeRatio', type='float', default=1.0,
                 help="ratio of independent SNPs (maize:0.32, sorghum:0.53)")
    p.add_option('--sort', default=False, action='store_true',
                 help="sort GWAS results based on chromosome and position")
    p.add_option('--idx_cols', default=None,
                 help="index of snp,chr,pos,pvalue if softare is other")
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    gwasfile, title = args
    if opts.software == 'other':
        if opts.usecols is None:
            sys.exit('--usecols must be specified if software is other')
        else:
            opts.usecols = list(map(int, opts.usecols.split(',')))
            print('indics of columns to be read: %s' % opts.usecols)

    gwas0 = ParseGWASresults(gwasfile, opts.software, sorting=opts.sort,
                             idx_cols=opts.idx_cols)
    df_plot = gwas0.df
    chr_pos_df = df_plot.groupby('chr').max()
    chr_lens = chr_pos_df['pos'].cumsum().tolist()[0:-1]
    chr_lens.insert(0, 0)
    chr_len_df = pd.DataFrame(chr_lens, index=chr_pos_df.index,
                              columns=['chr_offset']).reset_index()
    df_plot = df_plot.merge(right=chr_len_df, on='chr')
    df_plot['x'] = df_plot['pos'] + df_plot['chr_offset']
    df_plot = df_plot[['chr', 'pos', 'x', '-log10Pvalue']]

    plt.rcParams['xtick.labelsize'] = 8
    plt.rcParams['ytick.labelsize'] = 8
    plt.rcParams['xtick.major.pad'] = 1.5
    plt.rcParams['ytick.major.pad'] = 1.5
    colors = {1: '#a6cee3', 2: '#1f78b4', 3: '#b2df8a', 4: '#33a02c',
              5: '#fb9a99', 6: '#e31a1c', 7: '#fdbf6f', 8: '#ff7f00',
              9: '#cab2d6', 10: '#6a3d9a'}
    _, ax = plt.subplots(figsize=(6.8, 1.9))
    labels, label_POS = [], []
    color_idx = 1
    for chr, grp in df_plot.groupby(by='chr'):
        labels.append(chr)
        grp.plot(kind='scatter', x='x', y='-log10Pvalue', s=3.5,
                 color=colors[color_idx], ax=ax)
        label_pos = grp['x'].iloc[-1] - (grp['x'].iloc[-1] - grp['x'].iloc[0])/2
        label_POS.append(label_pos)
        color_idx += 1

    cutoff = -np.log10(opts.pvalue/(opts.MeRatio * gwas0.numberofSNPs))
    ax.axhline(cutoff, linestyle=':', linewidth=0.6, color='k', alpha=0.7)

    x_limit_right = df_plot['x'].iloc[-1]
    ax.set_xticks(label_POS)
    ax.set_xticklabels(labels)
    ax.set_xlim(-8000000, x_limit_right+8000000)

    y_limit_top = np.ceil(df_plot['-log10Pvalue'].max())
    y_ticks = np.arange(0, y_limit_top+1, 3)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_ticks.astype('int'))
    ax.set_ylim(-y_limit_top*0.04, y_limit_top+0.5)

    ax.set_ylabel(r'$\rm -log_{10}(\it{p})$', fontsize=12, labelpad=1)
    ax.set_xlabel('Chromosome', fontsize=12, labelpad=1)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title(title, fontsize=12)
    plt.tight_layout()
    plt.savefig(f'{title}_Manhattan.png', dpi=300)


def count_peaks(args):
    """
    %prog count_peaks GWASfile pvalue_cutoff output_prefix

    Identify peaks in the output of a GWAS run
    Args:
        GWASfile: GWAS file to be parsed
        pvalue: the pvalue cutoff in -log10
        output_prefix: The prefix of output files
    """
    p = OptionParser(count_peaks.__doc__)
    p.add_option('--software', default='gemma',
                 choices=('gemma', 'farmcpu', 'gapit', 'other'),
                 help='softare where the GWAS result came from. '
                      'If other, specify --usecols option')
    p.add_option('--WindowSize', type='int', default=150_000,
                 help='Maximum distance between two significant SNPs within '
                      'the same peak in base pairs')
    p.add_option('--sort', default=False, action='store_true',
                 help="sort GWAS file based on chromosome and pos")
    p.add_option('--idx_cols', default=None,
                 help="index of snp,chr,pos,pvalue if softare is other")
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    gwasfile, cutoff, outprefix = args

    if opts.software == 'other':
        if opts.usecols is None:
            sys.exit('--idx_cols option must be specified if software is other')
        else:
            idx_cols = list(map(int, opts.idx_cols.split(',')))
            print('indics of columns to be read: ', idx_cols)

    gwas0 = ParseGWASresults(gwasfile, opts.software, sorting=opts.sort,
                             idx_cols=idx_cols)
    df = gwas0.SignificantSNPs(sig_cutoff=float(cutoff))
    print('number of significant SNPs: ', df.shape[0])

    # find peaks in each chromosome
    peaks = []
    n = 0
    for idx, grp in df.groupby('chr'):
        print('chr: ', idx)
        last_pos = (-1 * opts.WindowSize) - 1
        for idx, row in grp.iterrows():
            pos = row['pos']
            if pos - last_pos < opts.WindowSize:
                peaks.append(n)
                last_pos = pos
            else:
                n += 1
                last_pos = pos
                peaks.append(n)
    df['peaks'] = peaks
    # summarize information in each peak
    f0 = open(f'{outprefix}.summary.csv', 'w')
    f0.write('Peak,#_of_SNPs,PeakStartPos,PeakStopPos,PeakLength,MostSigSNP,-log10Pvalue\n')
    for idx, grp in df.groupby('peaks'):
        st, ed = grp['pos'].min(), grp['pos'].max()
        distance = ed - st
        most_sig_snp = grp.loc[grp['-log10Pvalue'].idxmax(), 'snp']
        largest_pvalue = grp['-log10Pvalue'].max()
        snp_number = grp.shape[0]
        f0.write(f'{idx},{snp_number}, {st},{ed},{distance},{most_sig_snp},{largest_pvalue:.3f}\n')
    f0.close()
    df['-log10Pvalue'] = df['-log10Pvalue'].apply(lambda x: f'{x:%.3f}')
    df.to_csv(f'{outprefix}.csv', index=False)

    print(f'results have been saved to {outprefix}.csv and {outprefix}.summary.csv!')


def cal_MAF(args):
    """
    %prog cal_MAF SNPlist_file hmp_file

    Calculate MAF of SNPs lised in a tabular file
    """
    p = OptionParser(cal_MAF.__doc__)
    p.add_option('--header', default='no', choices=('yes', 'no'),
                 help='if there is a header in your SNP list file')
    p.add_option('--snp_col_idx', default='0',
                 help='specify the SNP column')
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    SNPlist_fn, hmp_fn = args
    df = pd.read_csv(SNPlist_fn, delim_whitespace=True, header=None) \
        if opts.header == 'no' \
        else pd.read_csv(SNPlist_fn, delim_whitespace=True)
    SNPs = df.iloc[:, int(opts.snp_col_idx)]
    SNPs.to_csv('SNPs_list.csv', index=False)
    cmd = f'grep -f SNPs_list.csv {hmp_fn} > Genotypes_list.csv'
    call(cmd, shell=True)
    f = open('Genotypes_list.csv')
    f1 = open(f'MAF.{SNPlist_fn}', 'w')
    f1.write('SNPs\tMAF\tCount\tMinorAlleleSMs\tMajorAlleleSMs\tHeteroSMs\n')
    header = np.array(open(hmp_fn).readline().split())
    for i in f:
        snp, maf, count, minor_idx, major_idx, hetero_idx = cal_MAF_row(i)
        minor_SMs = ','.join(list(header[minor_idx]))
        major_SMs = ','.join(list(header[major_idx]))
        hetero_SMs = ','.join(list(header[hetero_idx]))
        newi = '%s\t%s\t%s\t%s\t%s\t%s\n'\
            % (snp, maf, count, minor_SMs, major_SMs, hetero_SMs)
        f1.write(newi)
    f.close()
    f1.close()
    print(f'results have been saved to MAF.{SNPlist_fn}')


def fetch_EVs_from_FarmCPU(args):
    """
    %prog SNPlist_csv FarmCPUresult

    extract effect size of specified SNPs from FarmCPU result
    """
    p = OptionParser(fetch_EVs_from_FarmCPU.__doc__)
    p.add_option('--header', default='no', choices=('yes', 'no'),
                 help='specify if there is a header in your SNP list file')
    p.add_option('--col_idx', default='0',
                 help='specify the SNP column')
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    SNPlist_csv, farmResult = args
    df = pd.read_csv(SNPlist_csv, delim_whitespace=True, header=None) \
        if opts.header == 'no' \
        else pd.read_csv(SNPlist_csv, delim_whitespace=True)
    SNPs = df.iloc[:, int(opts.col_idx)]
    SNPs.to_csv('SNPs_list.csv', index=False)
    cmd = 'grep -f SNPs_list.csv {farmResult} > FarmCPU_list.csv'
    call(cmd, shell=True)
    f = open('FarmCPU_list.csv')
    f1 = open(f'EVs.{SNPlist_csv}', 'w')
    f1.write('SNPs\tEVs\n')
    for i in f:
        j = i.strip().split(',')
        snp, ev = j[0], j[-1]
        newi = f'{snp}\t{ev}\n'
        f1.write(newi)
    f.close()
    f1.close()
    print(f'results have been saved to EVs.{SNPlist_csv}')


def shared_sig_SNPs(args):
    """
    %prog SigSNPsFromGEMMA SigSNPsFromFarmcpu output

    find shared SNPs between GEMMA and FarmCPU
    """
    p = OptionParser(shared_sig_SNPs.__doc__)
    if len(args) == 0:
        sys.exit(not p.print_help())

    SigSNPsFromGEMMA, SigSNPsFromFarmcpu, output, = args
    df1 = pd.read_csv(SigSNPsFromFarmcpu, delim_whitespace=True)
    df2 = pd.read_csv(SigSNPsFromGEMMA, delim_whitespace=True)
    df = df2[df2['rs'].isin(df1['SNP'])]
    df.to_csv(output, index=False, sep='\t')
    print(f'shared significantly SNPs have saved to {output}')


def sig_SNPs(args):
    """
    %prog gwas_results output_fn

    extract the first top N significant SNPs from GWAS result
    """
    p = OptionParser(sig_SNPs.__doc__)
    p.add_option('--MeRatio', type='float', default=1.0,
                 help="specify the ratio of independent SNPs, maize is 0.32,"
                      " sorghum is 0.53")
    p.add_option('--chrom', default='all',
                 help="specify chromosome, e.g. 1, 2, 'all' for genome level")
    p.add_option('--software', default='gemma',
                 choices=('gemma', 'gapit', 'farmcpu', 'other'),
                 help='specify which software generates the GWAS result')
    p.add_option('--usecols', default=None,
                 help="specify index (0-based) of snp,chr,pos,pvalue "
                      "(comma separated without space) if softare is other")
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    gwasfile, output_fn, = args

    if opts.software == 'other':
        if opts.usecols is None:
            sys.exit('--usecols must be specified if software is other')
        else:
            opts.usecols = list(map(int, opts.usecols.split(',')))
    gwas0 = ParseGWASresults(gwasfile, opts.software, usecols=opts.usecols)
    df_significant = gwas0.SignificantSNPs(p_cutoff=0.05, MeRatio=opts.MeRatio)
    if opts.chrom != 'all':
        df_significant = df_significant[df_significant['chr'] == opts.chrom]
    df_significant.to_csv(output_fn, index=False, sep='\t')
    print(f'significant SNPs have saved to {output_fn}')


def linked_SNPs(args):
    """
    %prog input_SNPlist_file bed_prefix r2_cutoff output_prefix

    extract linked SNPs using plink.
    """
    p = OptionParser(linked_SNPs.__doc__)
    p.add_option('--col_idx', type='int', default=0,
                 help='specify which column contains SNP ID (0-based)')
    p.add_option('--header', default='yes', choices=('yes', 'no'),
                 help='add this option if there is no header in the input SNPlist file')
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='add this option to disable converting commands to slurm jobs')
    p.add_slurm_opts(job_prefix=linked_SNPs.__name__)
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    SNPlist_fn, bedprefix, cutoff, output_prefix, = args
    if opts.header == 'yes':
        df = pd.read_csv(SNPlist_fn, delim_whitespace=True,
                         usecols=[opts.col_idx])
    else:
        df = pd.read_csv(SNPlist_fn, delim_whitespace=True,
                         usecols=[opts.col_idx], header=None)
    pre = Path(SNPlist_fn).name.split('.')[0]
    df.to_csv(f'{pre}.SNPs_list.csv', index=False, header=None)
    cmd_local = f'{plink} --bfile {bedprefix} --r2 --ld-snp-list '\
                f'{pre}.SNPs_list.csv --ld-window-kb 5000 --ld-window 99999 '\
                f'--ld-window-r2 {cutoff} --noweb --out {output_prefix}\n'
    print('cmd on local:\n', cmd_local)

    cmd_header = 'ml plink'
    cmd_hcc = f'plink --bfile {bedprefix} --r2 --ld-snp-list '\
              f'{pre}.SNPs_list.csv --ld-window-kb 5000 --ld-window 99999 '\
              f'--ld-window-r2 {cutoff} --noweb --out {output_prefix}\n'
    print(f'cmd on HCC:\n{cmd_header}\n{cmd_hcc}')

    if not opts.disable_slurm:
        put2slurm_dict = vars(opts)
        put2slurm_dict['cmd_header'] = cmd_header
        create_slurm([cmd_hcc], put2slurm_dict)


def extract_geno_from_vcf(args):
    """
    %prog SNP_list_file VCF output

    extract genotypes for a bunch of SNPs from VCF
    """
    p = OptionParser(extract_geno_from_vcf.__doc__)
    p.add_option('--header', default='yes', choices=('yes', 'no'),
                 help='specify if there is a header in your SNP list file')
    p.add_option('--column', default='0',
                 help='specify which column is your SNP column 0-based')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    snplist, vcf, output, = args

    df = pd.read_csv(snplist, delim_whitespace=True) \
        if opts.header == 'yes' \
        else pd.read_csv(snplist, delim_whitespace=True, header=None)
    SNPs = df.iloc[:, int(opts.column)]
    SNP_keys = '\t'+SNPs+'\t'
    SNP_keys.to_csv('SNP_keys.csv', index=False)
    print('grep keys generated: SNP_keys.csv')

    cmd1 = f'zgrep -f SNP_keys.csv {vcf} > SNPs_keys.tmp.vcf' \
        if vcf.endswith('gz')\
        else f'grep -f SNP_keys.csv {vcf} > SNPs_keys.tmp.vcf'
    call(cmd1, shell=True)
    print('grep vcf done: SNPs_keys.tmp.vcf')

    cmd2 = f"zgrep -m 1 -P '#CHROM\tPOS' {vcf} > header.vcf" \
        if vcf.endswith('gz')\
        else f"zgrep -m 1 -P '#CHROM\tPOS' {vcf} > header.vcf"
    call(cmd2, shell=True)
    vcf_header = open('header.vcf')
    df_header = vcf_header.readline().split()
    print('header done: header.vcf')

    df_geno = pd.read_csv('SNPs_keys.tmp.vcf', delim_whitespace=True,
                          header=None)
    df_geno.columns = df_header
    df_geno0 = df_geno[['#CHROM','POS','ID','REF','ALT']]
    df_geno1 = df_geno[df_geno.columns[9:]]
    df_geno2 = df_geno1.applymap(lambda x: x.split(':')[0])
    df_geno_final = pd.concat([df_geno0, df_geno2], axis=1)
    df_geno_final.to_csv(output, index=False)
    print(f'genotype processing done: {output}')


def extract_gene_from_gff3(args):
    """
    %prog SNPlist gene.gff3

    extract gene names from gff3 file
    remember to customize the rule to get gene name, and the chr name in gff3
    """
    p = OptionParser(extract_gene_from_gff3.__doc__)
    p.add_option('--header', default='yes', choices=('yes', 'no'),
                 help='specify if there is a header in your SNP list file')
    p.add_option('--col_idx', default='1,2,0',
                 help='specify the index of Chr, Pos, SNP columns')
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    SNPlist, gff, = args

    df0 = pd.read_csv(SNPlist, header=None) \
        if opts.header == 'no' \
        else pd.read_csv(SNPlist)
    cols = [df0.columns[int(i)] for i in opts.col_idx.split(',')]
    df0 = df0[cols]
    df0.columns = ['chr', 'pos', 'snp']
    df0 = df0.sort_values(['chr', 'pos'])
    df0 = df0.reset_index(drop=True)
    print(df0)

    df1 = pd.read_csv(gff, sep='\t', header=None)
    # customize the rule right below
    df1['gene'] = df1.iloc[:, 8].apply(lambda x: x.split(';')[1].split('=')[1])
    df1 = df1[[0, 3, 4, 'gene']]
    df1.columns = ['chr', 'start', 'end', 'gene']
    print(df1.head())

    for g in df0.groupby('chr'):
        chrom = g[0]
        print('chr: ', chrom)
        SNPs = list(g[1]['snp'].unique())
        print('sig SNPs: ', SNPs)
        for pos in g[1]['pos']:
            print('SNP position: ', pos)
            # customize the chr name in gff3
            chrom_gff = 'Chr%02d' % chrom
            df2 = df1[df1['chr'] == chrom_gff]
            df2['gene_length'] = np.abs(df2['end'] - df2['start'])
            df2['st_dist'] = np.abs(pos - df2['start'])
            df2['ed_dist'] = np.abs(pos - df2['end'])
            df2['min_dist'] = df2[['st_dist', 'ed_dist']].min(axis=1)
            df2[df2['min_dist'] < 5000].to_csv(f'{chrom}_{pos}_genes.csv',
                                               index=False, sep='\t')


def extract_gene_func(args):
    """
    %prog GeneList FunctionFile output

    extract gene functions
    """
    p = OptionParser(extract_gene_func.__doc__)
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    genelist, FuncFile, output, = args
    cmd = f'grep -f {genelist} {FuncFile} > {output}'
    call(cmd, shell=True)
    print('Done! Check file: ', output)


def extract_protein_seq(args):
    """
    %prog GeneList seq_file output_prefix

    extract protein sequences of candidate genes
    """
    p = OptionParser(extract_protein_seq.__doc__)
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    genelist, SeqFile, out_prefix, = args
    cmd = f"grep '>' {SeqFile}|cut -d ' ' -f 1|cut -d '>' -f 2 > AllGene.names"
    call(cmd, shell=True)

    df_Genes = pd.read_csv(genelist, header=None)
    df_Trans = pd.read_csv('AllGene.names', header=None)
    df_Trans['gene'] = df_Trans[0].str.split('_').str.get(0)
    df1 = df_Trans[df_Trans['gene'].isin(df_Genes[0])]
    df1['gene'] = df1['gene'].astype('category')
    df1['gene'].cat.set_categories(df_Genes[0].tolist(), inplace=True)
    df2 = df1.sort_values(['gene', 0]).reset_index(drop=True)
    df2[0].to_csv(f'{out_prefix}.ProSeq.names', index=False, header=False)

    for i in list(df2[0]):
        print('fetching ', i)
        cmd = f"{faOneRecord} {SeqFile} {i} >> {out_prefix}.seqs"
        call(cmd, shell=True)
    print('Done!')


def plot_effect_size(args):
    """
    %prog EVlist(FarmCPU result) output_prefix
    plot the histogram of effect sizes
    """
    p = OptionParser(plot_effect_size.__doc__)
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    EVlist, output_prefix = args
    df = pd.read_csv(EVlist)
    EVs = df.iloc[:, -1]
    xlim = min(max(EVs), abs(min(EVs)))
    ax = EVs.plot(kind='hist', bins=60, grid=True, alpha=0.75, edgecolor='k')
    ax.set_xlim(-xlim, xlim)
    ax.set_xlabel('Effect size')
    ax.set_ylabel('Counts')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}.pdf')
    plt.savefig(f'{output_prefix}.png')


def plot_MAF(args):
    """
    %prog MAF_file1 MAF_file2 ... Label_of_MAF_file1 Label_of_MAF_file2 ...
    make density plot of MAFs (maximum of four MAF files) with the file name 'MAF_density.png'
    """
    p = OptionParser(plot_MAF.__doc__)
    p.add_option('--header', default='no', choices=('yes', 'no'),
                 help='specify if there is a header in your MAF file')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())

    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    fig, ax = plt.subplots(figsize=(4, 3.8))
    colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3']

    n = int(len(args)/2)
    fns, labels = args[:n], args[n:]
    print('MAF files: %s' % ' '.join(fns))
    print('Labels: %s' % ' '.join(labels))

    with tqdm(total=n) as pbar:
        for idx, (fn, label) in enumerate(zip(fns, labels)):
            df = pd.read_csv(fn, header=None) \
                if opts.header == 'no' else pd.read_csv(fn)
            df.columns = [label]
            ax = sns.kdeplot(
                df[label],
                shade=True,
                label=label,
                color=colors[idx],
                ax=ax,
                markerfacecolor='black',
                markersize=2)
            pbar.set_description(f'ploting {label}...')
            pbar.update(1)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('axes', -0.015))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('axes', -0.015))

    ax.set_xlim(0, 0.5)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Minor Allele Frequency', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    plt.tight_layout()
    plt.savefig('MAF_density.png', dpi=300)


def pdf_to_png(args):
    """
    %prog pdf2png dir_in dir_out

    convert pdf to png using imagemagick
    """
    p = OptionParser(pdf_to_png.__doc__)
    p.set_slurm_opts(jn=True)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())

    in_dir, out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')
    dir_path = Path(in_dir)
    pdfs = dir_path.glob('*.pdf')
    for pdf in pdfs:
        print(pdf)
        prf = pdf.name.replace('.pdf', '')
        png = pdf.name.replace('.pdf', '.png')
        header = SLURM_header % (100, 15000, prf, prf, prf)
        header += 'ml imagemagick\n'
        cmd = f'convert -density 300 {pdf} -resize 25% {out_path}/{png}\n'
        header += cmd
        with open(f'pdf2png.{prf}.slurm', 'w') as f:
            f.write(header)


if __name__ == '__main__':
    main()
