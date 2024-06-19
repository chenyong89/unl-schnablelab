"""
Perform data QC before running Genotype-Corrector
"""
import sys
import logging
import numpy as np
import pandas as pd

from pathlib import Path
from scipy.stats import chisquare

from .base import (eprint,
                   getChunk,
                   get_blocks,
                   sort_merge_sort,
                   bin_markers)
from optparse import OptionGroup, SUPPRESS_HELP
from schnablelab.apps.base import OptionParser, ActionDispatcher


def qc_missing(args):
    """
    %prog qc_missing input.matrix output.matrix

    QC on missing data
    """
    p = OptionParser(qc_missing.__doc__)
    p.add_option("-i", "--input", help=SUPPRESS_HELP)
    p.add_option("-o", "--output", help=SUPPRESS_HELP)
    p.add_option('--cutoff_snp', default=0.5, type='float',
                 help=('SNPs with missing rate higher than this cutoff will '
                       'be removed'))
    p.add_option('--rm_bad_samples', default=False, action="store_true",
                 help='whether remove samples with high missing ratio')
    p.add_option('--cutoff_sample', type='float',
                 help=('if --rm_bad_samples is True, sample whose missing '
                       'rate higher than this cutoff will be removed. Sample '
                       'missing rate will be calculated after flitering out '
                       'low quality SNPs'))
    q = OptionGroup(p, "format options")
    p.add_option_group(q)
    q.add_option('--homo1', default="A",
                 help='character for homozygous genotype')
    q.add_option("--homo2", default='B',
                 help="character for alternative homozygous genotype")
    q.add_option('--hete', default='X',
                 help='character for heterozygous genotype')
    q.add_option('--missing', default='-',
                 help='character for missing value')
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    if opts.rm_bad_samples and not opts.cutoff_sample:
        eprint('specify --cutoff_sample option if --rm_bad_samples is enabled')
        sys.exit(1)

    inmap, outmap = args
    inputmatrix = opts.input or inmap
    outputmatrix = opts.output or outmap

    chr_order, chr_nums = getChunk(inputmatrix)
    map_reader = pd.read_csv(inputmatrix, delim_whitespace=True,
                             index_col=[0, 1],  iterator=True)
    Good_SNPs = []
    for chrom in chr_order:
        print('{}...'.format(chrom))
        chunk = chr_nums[chrom]
        df_chr_tmp = map_reader.get_chunk(chunk)
        df_chr_tmp_num = df_chr_tmp.replace(
            [opts.homo1, opts.homo2, opts.hete, opts.missing], [0, 2, 1, 9])
        sample_num = df_chr_tmp_num.shape[1]
        good_rates = df_chr_tmp_num.apply(
            lambda x: (x == 9).sum()/sample_num, axis=1) <= opts.cutoff_snp
        good_snp = df_chr_tmp.loc[good_rates, :]
        Good_SNPs.append(good_snp)
    df1 = pd.concat(Good_SNPs)
    n0_snp = sum(chr_nums.values())
    n1_snp, n0_sm = df1.shape
    pct = n1_snp/float(n0_snp)*100
    print(f'{n0_snp} SNP markers before quality control.')
    print(f'{n1_snp}({pct:.1f}%) markers left after QC')

    if opts.rm_bad_samples:
        print('start quality control on samples')
        good_samples = df1.apply(
            lambda x: (x == opts.missing).sum()/n1_snp, axis=0) \
            <= opts.cutoff_sample
        df2 = df1.loc[:, good_samples]
        n1_sm = df2.shape[1]
        pct_sm = n1_sm/float(n0_sm)*100
        print(f'{n0_sm} samples before quality control.')
        print(f'{n1_snp}({pct_sm:.1f}%) markers left after QC')
        df2.to_csv(outputmatrix, sep='\t', index=True)
    else:
        df1.to_csv(outputmatrix, sep='\t', index=True)


def qc_sd(args):
    """
    %prog qc_sd input.matrix output.matrix

    run quality control on segregation distortions in each SNP.
    """
    p = OptionParser(qc_sd.__doc__)
    p.add_option("-i", "--input", help=SUPPRESS_HELP)
    p.add_option("-o", "--output", help=SUPPRESS_HELP)
    p.add_option('--population', default='RIL', choices=('RIL', 'F2', 'BCFn'),
                 help="population type")
    p.add_option('--sig_cutoff', default=1e-2, type='float',
                 help=('set the chi square test cutoff between 0 '
                       '(less strigent) to 1 (more strigent)'))
    q = OptionGroup(p, "format options")
    p.add_option_group(q)
    q.add_option('--homo1', default="A",
                 help='character for homozygous genotype')
    q.add_option("--homo2", default='B',
                 help="character for alternative homozygous genotype")
    q.add_option('--hete', default='X',
                 help='character for heterozygous genotype')
    q.add_option('--missing', default='-',
                 help='character for missing value')
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    inmap, outmap = args
    inputmatrix = opts.input or inmap
    outputmatrix = opts.output or outmap

    if opts.sig_cutoff >= 1 or opts.sig_cutoff <= 0:
        eprint('the sig_cutoff must be smaller than 1 and larger than 0')
        sys.exit(1)

    chr_order, chr_nums = getChunk(inputmatrix)
    map_reader = pd.read_csv(inputmatrix, delim_whitespace=True,
                             index_col=[0, 1], iterator=True)
    Good_SNPs = []
    for chrom in chr_order:
        print('{}...'.format(chrom))
        chunk = chr_nums[chrom]
        df_chr_tmp = map_reader.get_chunk(chunk)
        df_chr_tmp_num = df_chr_tmp.replace(
            [opts.homo1, opts.homo2, opts.hete, opts.missing], [0, 2, 1, 9])
        ob0, ob2 = (df_chr_tmp_num == 0).sum(axis=1), (df_chr_tmp_num == 2)\
            .sum(axis=1)
        obsum = ob0 + ob2
        exp0, exp2 = (obsum*0.75, obsum*0.25) if opts.population == 'BCFn' \
            else (obsum*0.5, obsum*0.5)
        df_chi = pd.DataFrame(dict(zip(['ob0', 'ob2', 'exp0', 'exp2'],
                                       [ob0, ob2, exp0, exp2])))
        min_cond = ((df_chi['ob0'] > 5) & (df_chi['ob2'] > 5)).values
        pval_cond = chisquare(
            [df_chi['ob0'], df_chi['ob2']],
            [df_chi['exp0'], df_chi['exp2']]).pvalue >= opts.sig_cutoff
        good_snp = df_chr_tmp.loc[(min_cond & pval_cond), :]
        Good_SNPs.append(good_snp)
    df1 = pd.concat(Good_SNPs)
    n0_snp = sum(chr_nums.values())
    n1_snp = df1.shape[0]
    pct = n1_snp/float(n0_snp)*100
    print(f'{n0_snp} SNP markers before quality control.')
    print(f'{n1_snp}({pct:.1f}%) markers left after the quality control.')
    df1.to_csv(outputmatrix, sep='\t', index=True)


def qc_hetero(args):
    """
    %prog qc_hetero input.matrix output.matrix

    quality control on the continuous same homozygous in heterozygous region.
    """
    p = OptionParser(qc_hetero.__doc__)
    p.add_option("-i", "--input", help=SUPPRESS_HELP)
    p.add_option("-o", "--output", help=SUPPRESS_HELP)
    p.add_option("--read_len", default=100, type='int',
                 help="read length for SNP calling")
    p.add_option("--logfile", default='GC.qc_hetero.info',
                 help="specify the file saving binning info")
    q = OptionGroup(p, "format options")
    p.add_option_group(q)
    q.add_option('--homo1', default="A",
                 help='character for homozygous genotype')
    q.add_option("--homo2", default='B',
                 help="character for alternative homozygous genotype")
    q.add_option('--hete', default='X',
                 help='character for heterozygous genotype')
    q.add_option('--missing', default='-',
                 help='character for missing value')
    r = OptionGroup(p, 'advanced options')
    p.add_option_group(r)
    r.add_option('--nonhetero_lens', default=8, type='int',
                 help=('number of homozygous between two heterozygous in '
                       'a heterozygous region'))
    r.add_option('--min_homo', default=2, type='int',
                 help=('number of continuous homozygous within the read '
                       'length in a heterozygous region'))
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    inmap, outmap = args
    inputmatrix = opts.input or inmap
    outputmatrix = opts.output or outmap

    logging.basicConfig(filename=opts.logfile,
                        level=logging.DEBUG,
                        format="%(asctime)s:%(levelname)s:%(message)s")

    chr_order, chr_nums = getChunk(inputmatrix)
    map_reader = pd.read_csv(inputmatrix, delim_whitespace=True,
                             index_col=[0, 1], iterator=True)
    Good_SNPs = []
    for chrom in chr_order:
        print('{}...'.format(chrom))
        logging.debug(chrom)
        chunk = chr_nums[chrom]
        df_chr_tmp = map_reader.get_chunk(chunk)
        _ = df_chr_tmp.index
        df_chr_tmp_num = df_chr_tmp.replace(
            [opts.homo1, opts.homo2, opts.hete, opts.missing],
            [0, 2, 1, 9]).loc[chrom]
        chr_bin_ids = []
        for sm in df_chr_tmp_num:
            geno_grouper = (df_chr_tmp_num[sm].diff(1) != 0)\
                .astype('int').cumsum()
            idx, geno, lens = [], [], []
            for __, grp_geno in df_chr_tmp_num[sm].groupby(geno_grouper):
                idx.append(grp_geno.index)
                geno.append(grp_geno.unique()[0])
                lens.append(grp_geno.shape[0])
            df_grp_geno = pd.DataFrame(dict(zip(['idx', 'geno', 'lens'],
                                                [idx, geno, lens])))
            # 1: hetero genotype 0: others(homo1, homo2, missing)
            df_grp_geno['type'] = df_grp_geno['geno'].apply(
                lambda x: 1 if x == 1 else 0)
            type_grouper = (df_grp_geno['type'].diff(1) != 0)\
                .astype('int').cumsum()
            for __, grp_type in df_grp_geno.groupby(type_grouper):
                if grp_type['type'].unique()[0] == 0:
                    nonhetero_lens = grp_type['lens'].sum()
                    if nonhetero_lens <= opts.nonhetero_lens:
                        for __, row in grp_type.iterrows():
                            if row.geno == 0 or row.geno == 2:
                                bin_ids = get_blocks(row['idx'].values,
                                                     dist=opts.read_len,
                                                     block_size=opts.min_homo)
                                if bin_ids:
                                    for bin_index in bin_ids:
                                        if bin_index not in chr_bin_ids:
                                            chr_bin_ids.append(bin_index)
        if chr_bin_ids:
            dropping_ids = []
            merged_bin_ids = sort_merge_sort(chr_bin_ids)
            for idx_block in merged_bin_ids:
                logging.debug('positions: {}'.format(idx_block))
                genos_block = df_chr_tmp_num.loc[idx_block, :]
                missings = genos_block.apply(lambda x: (x == 9).sum(), axis=1)
                heteros = genos_block.apply(lambda x: (x == 1).sum(), axis=1)
                dropping_index = list(pd.concat([missings, heteros], axis=1)
                                      .sort_values([0, 1]).index[1:])
                dropping_ids.extend(dropping_index)
            df_chr_tmp = df_chr_tmp.drop(dropping_ids, level=1)
        Good_SNPs.append(df_chr_tmp)
    df1 = pd.concat(Good_SNPs)
    n0_snp = sum(chr_nums.values())
    n1_snp = df1.shape[0]
    pct = n1_snp/float(n0_snp)*100
    print(f'{n0_snp} SNP markers before quality control.')
    print(f'{n1_snp}({pct:.1f}%) markers left after the quality control.')
    df1.to_csv(outputmatrix, sep='\t', index=True)
    print('Done! Check {} for running details.'.format(opts.logfile))


def sort_pos(args):
    """
    %prog sort_pos input.map output.sorted.map

    sort markers based on position 
    """
    p = OptionParser(sort_pos.__doc__)
    opts, args = p.parse_args(args)
    if len(args) != 2:
        sys.exit(not p.print_help())

    inmap, outmap, = args
    if Path(outmap).exists():
        eprint(f"ERROR: The output file `{outmap}` already exists")
        sys.exit(1)
    df = pd.read_csv(inmap, delim_whitespace=True)
    idx_col = list(df.columns[0:2])
    df.sort_values(idx_col).to_csv(outmap, sep='\t', index=False)
    print(f'Done! Check {outmap} file.')


def vcf2map(args):
    """
    %prog vcf2map input.vcf output.matrix

    convert vcf format to genotype matrix format
    """
    p = OptionParser(vcf2map.__doc__)
    p.add_option("-i", "--input", help=SUPPRESS_HELP)
    p.add_option("-o", "--output", help=SUPPRESS_HELP)
    p.add_option('--homo1', default="A",
                 help='character for homozygous genotype')
    p.add_option("--homo2", default='B',
                 help="character for alternative homozygous genotype")
    p.add_option('--hete', default='X',
                 help='character for heterozygous genotype')
    p.add_option('--missing', default='-',
                 help='character for missing value')
    p.add_option("--logfile", default='GC.vcf2map.info',
                 help="specify the log file")
    opts, args = p.parse_args(args)
    if len(args) != 2:
        sys.exit(not p.print_help())

    invcf, outmap = args
    inputvcf = opts.input or invcf
    outputmatrix = opts.output or outmap

    if Path(outputmatrix).exists():
        eprint(f"ERROR: The output file `{outputmatrix}` already exists")
        sys.exit(1)

    logging.basicConfig(filename=opts.logfile,
                        level=logging.DEBUG,
                        format="%(asctime)s:%(levelname)s:%(message)s")

    right_gt = {'0|0': opts.homo1, '0/0': opts.homo1,
                '0|1': opts.hete, '1|0': opts.hete,
                '0/1': opts.hete, '1/0': opts.hete,
                '1|1': opts.homo2, '1/1': opts.homo2,
                '.|.': opts.missing, './.': opts.missing, '.': opts.missing}
    useless_cols = ['ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    index_cols = ['#CHROM', 'POS']
    vcffile = open(inputvcf)
    n = 0
    for i in vcffile:
        if i.startswith('##'):
            n += 1
        else:
            break
    vcffile.close()
    chr_order, chr_nums = getChunk(inputvcf, ignore=n+1)
    vcf_reader = pd.read_csv(inputvcf, header=n, delim_whitespace=True,
                             usecols=lambda x: x not in useless_cols,
                             iterator=True)
    tmp_chr_list = []
    for chrom in chr_order:
        logging.debug('{}...'.format(chrom))
        print('{}...'.format(chrom))
        chunk = chr_nums[chrom]
        df_chr_tmp = vcf_reader.get_chunk(chunk)
        df_chr_tmp = df_chr_tmp.set_index(index_cols)
        df_chr_tmp = df_chr_tmp.applymap(lambda x: x.split(':')[0])
        df_chr_tmp = df_chr_tmp.applymap(
            lambda x: right_gt[x] if x in right_gt else np.nan)
        df_chr_tmp.dropna(inplace=True)
        tmp_chr_list.append(df_chr_tmp)
    df1 = pd.concat(tmp_chr_list)
    df1.to_csv(outputmatrix, sep='\t', index=True)
    vcffile.close()
    print('Done!')


def qc_dup(args):
    """
    %prog qc_dup corrected.matrix output.matrix

    QC on duplicated markers by by merging consecutive markers with same
    genotype
    """
    p = OptionParser(qc_dup.__doc__)
    p.add_option("-i", "--input", help=SUPPRESS_HELP)
    p.add_option("-o", "--output", help=SUPPRESS_HELP)
    p.add_option('--diff_num', default=0, type='int',
                 help=('number of different genotypes between two consecutive '
                       'markers less than or equal to this value will be '
                       'merged. missing values will not be counted.'))
    p.add_option('--missing', default='-',
                 help='character for missing value in genotype matrix file')
    p.add_option("--logfile", default='GC.bin.log',
                 help="specify the file saving running info")
    opts, args = p.parse_args(args)
    if len(args) != 2:
        sys.exit(not p.print_help())

    inmap, outmap = args
    inputmatrix = opts.input or inmap
    outputmatrix = opts.output or outmap

    if Path(outputmatrix).exists():
        eprint(f"ERROR: The output file `{outputmatrix}` already exists")
        sys.exit(1)

    chr_order, chr_nums = getChunk(inputmatrix)
    map_reader = pd.read_csv(inputmatrix, delim_whitespace=True,
                             index_col=[0, 1],  iterator=True)
    Good_SNPs = []
    binning_info = []
    for chrom in chr_order:
        print('{}...'.format(chrom))
        chunk = chr_nums[chrom]
        df_chr_tmp = map_reader.get_chunk(chunk)
        if df_chr_tmp.shape[0] == 1:
            Good_SNPs.append(df_chr_tmp)
        else:
            represent_idx, block_idx, results = bin_markers(
                df_chr_tmp.loc[chrom], diff=opts.diff_num,
                missing_value=opts.missing)
            good_snp = df_chr_tmp.loc[(chrom, results), :]
            Good_SNPs.append(good_snp)
            if represent_idx:
                df_binning_info = pd.DataFrame(dict(zip(
                    ['chr', 'representative_marker', 'markers'],
                    [chrom, represent_idx, block_idx])))
                binning_info.append(df_binning_info)
    df1 = pd.concat(Good_SNPs)
    df1.to_csv(outputmatrix, sep='\t', index=True)
    n0_snp = sum(chr_nums.values())
    n1_snp = df1.shape[0]
    pct = n1_snp/float(n0_snp)*100
    print(f'{n0_snp} SNP markers before compression.')
    print(f'{n1_snp}({pct:.1f}%) markers left after compression.')

    if binning_info:
        df2 = pd.concat(binning_info)
        df2.to_csv(opts.logfile, sep='\t', index=False)
        print('Check {} for binning details.'.format(opts.logfile))


def main():
    actions = (
        ('vcf2map', 'convert vcf to genotype matrix file'),
        ('sort_pos', 'sort makers based on position'),
        ('qc_missing', 'QC on missing gneotypes'),
        ('qc_sd', 'QC on segregation distortions'),
        ('qc_hetero', 'QC on heterozygous region'),
        ('qc_dup', 'QC on duplicated markers with same genotype'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == '__main__':
    main()
