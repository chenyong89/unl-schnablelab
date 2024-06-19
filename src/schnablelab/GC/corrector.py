"""
Correct Genotype Calls in biparental populations
"""
import re
import sys
import logging
import pandas as pd
import numpy as np

from pathlib import Path
from collections import Counter
from scipy.stats import binom
from optparse import OptionGroup, SUPPRESS_HELP

from .base import eprint, getChunk, ParseConfig
from schnablelab.apps.base import OptionParser, ActionDispatcher


class FixSliding():
    def __init__(self, seqnum, initial_corrected_seq, h_island_idx):
        self.homo_pattern1 = re.compile('9*09*09*|9*29*29*')
        self.hete_pattern = re.compile('1|9*09*29*|9*29*09*')
        self.homo_pattern2 = re.compile(
            '9*29*29*2*9*2*9*2*9*2*9*2*9*|9*09*09*0*9*0*9*0*9*0*9*0*9*')
        self.seqnum = seqnum
        self.initial_corrected_seq = initial_corrected_seq
        self.h_island_idx = h_island_idx

    def case1(self):
        """
        original seq: hhhh,aaaa
        after initial correction: hhh,aaaaa
        """
        st = self.h_island_idx[-1]+1
        ed = st+6
        if ed <= self.seqnum.index[-1]:
            indent_genos = ''.join(map(str, self.seqnum.loc[st: ed].values))
            result = self.homo_pattern1.search(indent_genos)
            try:
                i = result.start()
                if i > 0:
                    self.initial_corrected_seq.loc[st:st+i-1] = 1
            except AttributeError:
                pass

    def case2(self):
        """
        original seq: hhhh,aaaa
        after initial correction: hhhhh,aaa
        """
        ed = self.h_island_idx[-1]
        st = ed - 6
        if st >= self.seqnum.index[0]:
            indent_genos = ''.join(map(str, self.seqnum.loc[st: ed].values))
            result = self.homo_pattern1.search(indent_genos)
            try:
                i = result.start()
                if i > 0:
                    if self.hete_pattern.search(indent_genos, i) is None:
                        geno = 0 if '0' in result.group() else 2
                        self.initial_corrected_seq.loc[st+i:ed] = geno
            except AttributeError:
                pass

    def case3(self):
        """
        original seq: aaaa,hhhh
        after initial correction: aaa,hhhhh
        """
        st = self.h_island_idx[0]
        ed = st+6
        if ed <= self.h_island_idx[-1]:
            indent_genos = ''.join(map(str, self.seqnum.loc[st: ed].values))
            result = self.homo_pattern2.match(indent_genos)
            try:
                i = result.end()
                if i > 0:
                    geno = 0 if '0' in result.group() else 2
                    self.initial_corrected_seq.loc[st:st+i-1] = geno
            except AttributeError:
                pass

    def case4(self):
        """
        original seq: aaaa,hhhh
        after initial correction: aaaaa,hhh
        """
        ed = self.h_island_idx[0]-1
        st = ed - 6
        if st >= self.seqnum.index[0]:
            indent_genos = ''.join(map(str, self.seqnum.loc[st: ed].values))
            result = self.homo_pattern2.search(indent_genos)
            try:
                i = result.end()
                if i < 7:
                    self.initial_corrected_seq.loc[st+i:ed] = 1
            except AttributeError:
                pass


def get_corrected_num(seqnum, corrected_seq):
    """
    count number of genotpes corrected
    """
    return (seqnum != corrected_seq).sum()


class Correction(object):
    """
    This class contains steps to correct the original seq per sample
    """
    def __init__(self, config_params, orig_seq_without_idx_num,
                 random_seed=None):
        self.cargs = config_params
        self.seq_num = orig_seq_without_idx_num
        self.seq_num_no1 = self.seq_num.copy()
        if random_seed is not None:
            np.random.seed(random_seed)
        self.seq_num_no1[self.seq_num_no1 == 1] = np.random.choice(
            [0, 2], (self.seq_num_no1 == 1).sum())

    def get_rolling_geno(self):
        rolling_geno = self.seq_num_no1.rolling(
            self.cargs.win_size, center=True)\
            .apply(self.get_mid_geno, raw=True, args=(self.cargs,))
        return rolling_geno

    def get_rolling_score(self):
        rolling_score = self.seq_num_no1.rolling(
            self.cargs.win_size, center=True)\
            .apply(self.get_score, raw=True, args=(self.cargs,))
        return rolling_score

    def get_corrected(self):
        '''
        fix the h island
        '''
        rolling_geno = self.get_rolling_geno()
        rolling_score = self.get_rolling_score()
        grouper = (rolling_geno.diff(1) != 0).astype('int').cumsum()
        for __, grp in rolling_geno.groupby(grouper):
            geno = grp.unique()[0]
            if geno == 1:
                h_score_island = rolling_score[grp.index]
                if self.judge_h_island(h_score_island):
                    # adjust the sliding problem of the real h island.
                    # h_island_len = grp.shape[0]
                    fix_slide = FixSliding(self.seq_num,
                                           rolling_geno,
                                           grp.index)
                    fix_slide.case1()
                    fix_slide.case2()
                    fix_slide.case3()
                    fix_slide.case4()
                else:
                    # call back the fake h island to the origianl genotypes
                    real_genos = self.callback(h_score_island)
                    rolling_geno[grp.index] = real_genos

        # substitute NaNs with original genotypes
        corrected = rolling_geno.where(
            ~rolling_geno.isna(), other=self.seq_num).astype('int')
        return corrected

    def judge_h_island(self, h_scores):
        """
        judge if the h island in the corrected seqs is real or fake.
        """
        length = h_scores.shape[0]
        if length >= 3:
            trends = h_scores.diff().iloc[1:]
            ups, downs = (trends > 0).sum(), (trends < 0).sum()
            if ups == 0 or downs == 0:
                # fake h island because the score curve changes monotonically
                return False
            else:
                # real h island
                return True
        else:
            # real h island
            return True

    def _count_genos(self, np_array):
        """
        count genotypes in a given seq
        """
        counts = Counter(np_array)
        return counts[0], counts[2], counts[9]

    def get_mid_geno(self, np_array, cargs_obj):
        """
        return the genotype with highest probability in the central.
        """
        a_count, b_count, miss_count = self._count_genos(np_array)
        ab_count = a_count + b_count
        if ab_count > cargs_obj.win_size//2:
            a_ex_prob = binom.pmf(b_count, ab_count, cargs_obj.error_a)
            h_ex_prob = binom.pmf(b_count, ab_count, 0.5+cargs_obj.error_h/2)
            b_ex_prob = binom.pmf(b_count, ab_count, 1-cargs_obj.error_b)
            d = {key: value for (key, value) in
                 zip([0, 1, 2], [a_ex_prob, h_ex_prob, b_ex_prob])}
            return max(d, key=d.get)
        else:
            return np.nan

    def get_score(self, np_array, cargs_obj):
        """
        calculate the score for each sliding window in the seq_num_no1
        """
        a_count, b_count, __ = self._count_genos(np_array)
        if a_count+b_count > cargs_obj.win_size//2:
            return a_count/float(b_count) if b_count != 0 \
                else a_count/(b_count+0.1)
        else:
            return np.nan

    def callback(self, h_scores):
        """
        call back the h island to the origianl genotypes based on the socre
        """
        realgenos = []
        for val in h_scores:
            if val > 1:
                realgenos.append(0)
            elif val < 1:
                realgenos.append(2)
            else:
                realgenos.append(np.nan)
        return realgenos


def correct(args):
    """
    %prog correct config.txt input.matrix

    Correct genotype calls and impute missing values in biparental populations
    """
    p = OptionParser(correct.__doc__)
    p.add_option("-c", "--configfile", help=SUPPRESS_HELP)
    p.add_option("-m", "--matrixfile", help=SUPPRESS_HELP)
    p.add_option('--itertimes', default=7, type='int',
                 help='maximum rounds of corrections')
    p.add_option('--random_seed', type='int',
                 help='specify random seed')
    q = OptionGroup(p, "output options")
    p.add_option_group(q)
    q.add_option('--opp',
                 help='the prefix of the output filenames')
    q.add_option("--logfile", default='GC.correct.log',
                 help="specify the file saving running info")
    q.add_option('--debug', default=False, action="store_true",
                 help=('generate a tmp file containing both original and '
                       'corrected genotypes for further manual correction'))

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    configfile, mapfile = args
    inputmatrix = opts.matrixfile or mapfile
    inputconfig = opts.configfile or configfile

    opf = inputmatrix.rsplit(".", 1)[0]+'.corrected.map' \
        if opts.opp is None else '{}.map'.format(opts.opp)
    if Path(opf).exists():
        eprint(f"ERROR: The output file `{opf}` already exists")
        sys.exit(1)

    logging.basicConfig(filename=opts.logfile,
                        level=logging.DEBUG,
                        format="%(asctime)s:%(levelname)s:%(message)s")

    cargs = ParseConfig(inputconfig)
    if cargs.win_size % 2 == 0:
        eprint("ERROR: The slding window value cannot be even")
        sys.exit(1)
    logging.debug("Parameters in config file: {0}".format(cargs.__dict__))

    chr_order, chr_nums = getChunk(inputmatrix)
    map_reader = pd.read_csv(inputmatrix, delim_whitespace=True,
                             index_col=[0, 1],  iterator=True)
    tmp_chr_list = []
    for chrom in chr_order:
        logging.debug('{}...'.format(chrom))
        print('{}...'.format(chrom))
        chunk = chr_nums[chrom]
        df_chr_tmp = map_reader.get_chunk(chunk)
        marker_num, sample_num = df_chr_tmp.shape
        logging.debug('{} contains {} markers and {} samples.'.format(
            chrom, marker_num, sample_num))
        tmp_sm_list = []
        for sm in df_chr_tmp:
            logging.debug('Start correcting {}...'.format(sm))
            orig_seq = df_chr_tmp[sm]
            orig_idx = orig_seq.index
            seq_no_idx = orig_seq.reset_index(drop=True)
            seq_no_idx_num = seq_no_idx.replace(
                [cargs.gt_a, cargs.gt_b, cargs.gt_h, cargs.gt_miss],
                [0, 2, 1, 9])
            if seq_no_idx_num.shape[0] <= cargs.win_size:
                logging.debug('num of markers smaller than predefined window '
                              'size, skip...')
                final_seq_no_idx = seq_no_idx
            else:
                logging.debug('correction round 1...')
                _correction = Correction(cargs, seq_no_idx_num,
                                         opts.random_seed)
                corrected_seq = _correction.get_corrected()
                corrected_n = get_corrected_num(seq_no_idx_num,
                                                corrected_seq)
                round_n = 2
                while round_n <= opts.itertimes:
                    logging.debug(f'correction round {round_n}...')
                    _correction = Correction(cargs, seq_no_idx_num,
                                             opts.random_seed)
                    corrected_seq = _correction.get_corrected()
                    corrected_n1 = get_corrected_num(
                        seq_no_idx_num, corrected_seq)
                    round_n += 1
                    if (corrected_n1 - corrected_n)/(corrected_n+0.01) <= 0.01:
                        break
                    else:
                        corrected_n = corrected_n1
                final_seq_no_idx = corrected_seq.replace(
                    [0, 2, 1, 9],
                    [cargs.gt_a, cargs.gt_b, cargs.gt_h, cargs.gt_miss])
            final_seq_no_idx.index = orig_idx
            final_seq = final_seq_no_idx
            tmp_sm_list.append(final_seq)
        df_sm_tmp = pd.concat(tmp_sm_list, axis=1)
        tmp_chr_list.append(df_sm_tmp)
    df_corrected = pd.concat(tmp_chr_list)
    df_corrected.to_csv(opf, sep='\t', index=True)

    if opts.debug:
        logging.debug('generating the tmp file for debug use...')
        df_uncorrected = pd.read_csv(inputmatrix, delim_whitespace=True,
                                     index_col=[0, 1])
        df_debug = df_corrected.where(
            df_corrected == df_uncorrected,
            other=df_corrected+'('+df_uncorrected+')')
        df_debug.to_csv(opf+'.debug', sep='\t', index=True)
    print('Done!')


def cleanup(args):
    """
    %prog cleanup tmp.matrix out.matrix

    remove redundant info for Debug in the temporary matrix file
    """
    p = OptionParser(cleanup.__doc__)
    p.add_option("-i", "--input", help=SUPPRESS_HELP)
    p.add_option("-o", "--output", help=SUPPRESS_HELP)
    opts, args = p.parse_args(args)
    if len(args) != 2:
        sys.exit(not p.print_help())

    inmap, outmap = args
    inputmatrix = opts.input or inmap
    outputmatrix = opts.output or outmap

    df = pd.read_csv(inputmatrix, delim_whitespace=True, index_col=[0, 1])
    df.applymap(lambda x: x.split('(')[0]).to_csv(outputmatrix, sep='\t',
                                                  index=True)
    print('Done!')


MST_header = """population_type {}
population_name {}
distance_function {}
cut_off_p_value {}
no_map_dist {}
no_map_size {}
missing_threshold {}
estimation_before_clustering {}
detect_bad_data {}
objective_function {}
number_of_loci {}
number_of_individual {}

"""

mst_homos, mst_miss, mst_hete = ('a', 'A', 'b', 'B'), ('-', 'U'), 'X'


def format(args):
    """
    %prog format corrected.matrix

    convert corrected genotype matrix file to other formats
    Example:
    command:
    `python -m schnablelab.genotype_corrector.gc format test.map --mstmap --mstmap_pop_type RIL2`
    will generate `test.mstmap` for MSTmap.
    """
    p = OptionParser(format.__doc__)
    p.add_option("-i", "--input", help=SUPPRESS_HELP)
    p.add_option("--mstmap",  default=False, action="store_true",
                 help='convert to MSTmap format')
    p.add_option("--rqtl",  default=False, action="store_true",
                 help='convert to R/qtl format')
    p.add_option("--joinmap",  default=False, action="store_true",
                 help='convert to JoinMap format')

    q = OptionGroup(p, "format options for input matrix file")
    p.add_option_group(q)
    q.add_option('--homo1', default="A", choices=('a', 'A'),
                 help='character for homozygous genotype')
    q.add_option("--homo2", default='B', choices=('b', 'B'),
                 help="character for alternative homozygous genotype")
    q.add_option('--hete', default='X', choices=('h', 'H', 'X'),
                 help='character for heterozygous genotype')
    q.add_option('--missing', default='-', choices=('-', 'U'),
                 help='character for missing value')

    r = OptionGroup(p, "parameters for MSTmap")
    p.add_option_group(r)
    r.add_option('--mstmap_pop_type',
                 help=('values can be DH or RILd where d is a natural number. '
                       'e.g., RIL6 means a RIL population at generation 6.'
                       'Use RIL2 for F2. Use DH for BC1, DH and Hap.'))
    r.add_option("--population_name", default='LinkageGroup',
                 help='name for the mapping population')
    r.add_option('--distance_function', default='kosambi',
                 choices=('kosambi', 'haldane'),
                 help="choose Kosambi's and Haldane's distance functions")
    r.add_option('--cut_off_p_value', default=0.000001,
                 help=('specifies the threshold to be used for clustering the '
                       'markers into LGs'))
    r.add_option('--no_map_dist', default=15,
                 help='check mstmap manual for details')
    r.add_option('--no_map_size', default=5,
                 help='check mstmap manual for details')
    r.add_option('--missing_threshold', default=0.4,
                 help=('any marker with more than this value will be removed '
                       'completely without being mapped'))
    r.add_option('--estimation_before_clustering', default='no',
                 choices=('yes', 'no'),
                 help=('if yes, MSTmap will try to estimate missing data '
                       'before clustering the markers into linkage groups'))
    r.add_option('--detect_bad_data', default='yes', choices=('yes', 'no'),
                 help='if yes turn on the error detection feature in MSTmap')
    r.add_option('--objective_function', default='COUNT',
                 choices=('COUNT', 'ML'),
                 help='specifies the objective function')

    s = OptionGroup(p, "parameters for JoinMap or R/qtl")
    p.add_option_group(s)
    s.add_option('--pop_type', default='RIL', choices=('RIL', 'F2'),
                 help='specify mapping population type')

    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(not p.print_help())

    inmap, = args
    inputmatrix = opts.input or inmap

    if (not opts.rqtl) and (not opts.joinmap) and (not opts.mstmap):
        eprint("ERROR: add at least one output format option.")
        sys.exit(1)

    if opts.mstmap:
        if not opts.mstmap_pop_type:
            eprint("ERROR: please choose population type for mstmap format.")
            sys.exit(1)
        if not (opts.mstmap_pop_type.startswith('RIL') or
                opts.mstmap_pop_type == 'DH'):
            eprint('ERROR: only RILd and DH supported in MSTmap')
            sys.exit(1)

        opf = inputmatrix.rsplit(".", 1)[0]+'.mstmap'  # output file
        if Path(opf).exists():
            eprint(f"ERROR: Output file `{opf}` already exists!")
            sys.exit(1)

        df = pd.read_csv(inputmatrix, delim_whitespace=True)
        cols = list(df.columns[2:])
        cols.insert(0, 'locus_name')
        df['locus_name'] = df.iloc[:, 0].astype('str') + '_' + \
            df.iloc[:, 1].astype('str')
        df = df[cols]
        print(df.head())
        snp_num, sm_num = df.shape[0], df.shape[1]-1
        f1 = open(opf, 'w')
        f1.write(MST_header.format(
            opts.mstmap_pop_type, opts.population_name, opts.distance_function,
            opts.cut_off_p_value, opts.no_map_dist, opts.no_map_size,
            opts.missing_threshold, opts.estimation_before_clustering,
            opts.detect_bad_data, opts.objective_function, snp_num, sm_num))
        f1.close()

        df.to_csv(opf, sep='\t', index=False, mode='a')
        print('Done, check file {}!'.format(opf))

    if opts.joinmap:
        opf = inputmatrix.rsplit(".", 1)[0]+'.joinmap.xlsx'  # output file
        if Path(opf).exists():
            eprint(f"ERROR: Output file `{opf}` already exists!")
            sys.exit(1)

        df = pd.read_csv(inputmatrix, delim_whitespace=True)
        need_reps, reps = [], []
        if opts.homo1 != 'a':
            need_reps.append(opts.homo1)
            reps.append('a')
        if opts.homo2 != 'b':
            need_reps.append(opts.homo2)
            reps.append('b')
        if opts.hete != 'h':
            need_reps.append(opts.hete)
            reps.append('h')
        if opts.missing != '-':
            need_reps.append(opts.missing)
            reps.append('-')
        if need_reps:
            df = df.replace(need_reps, reps)

        cols = list(df.columns[2:])
        cols.insert(0, 'Classification')
        cols.insert(0, 'locus_name')
        df['locus_name'] = df.iloc[:, 0].astype('str') + '_' + \
            df.iloc[:, 1].astype('str')
        df['Classification'] = '(a,h,b)'
        df = df[cols]
        df.to_excel(opf)
        print('Done! Now you can load the genotype data into a JoinMap '
              f'project from the spreadsheet {opf} to a dataset node.')

    if opts.rqtl:
        opf = inputmatrix.rsplit(".", 1)[0]+'.rqtl.csv'  # output file
        if Path(opf).exists():
            eprint(f"ERROR: Output file `{opf}` already exists!")
            sys.exit(1)

        df = pd.read_csv(inputmatrix, delim_whitespace=True)
        need_reps, reps = [], []
        if opts.homo1 != 'A':
            need_reps.append(opts.homo1)
            reps.append('A')
        if opts.homo2 != 'B':
            need_reps.append(opts.homo2)
            reps.append('B')
        if opts.hete != 'H':
            need_reps.append(opts.hete)
            reps.append('H')
        if opts.missing != '-':
            need_reps.append(opts.missing)
            reps.append('-')
        if need_reps:
            df = df.replace(need_reps, reps)

        cols = list(df.columns[2:])
        cols.insert(0, 'id')
        cols.insert(1, 'group')
        df['id'] = df.iloc[:, 0].astype('str') + '_' + \
            df.iloc[:, 1].astype('str')
        df['group'] = 1
        df = df[cols]

        df = df.set_index('id')
        df1 = df.transpose()
        df1 = df1.reset_index()
        columns = list(df1.columns)
        columns[0] = 'id'
        df1.columns = columns

        df1.iloc[0, 0] = np.nan
        df1.to_csv(opf, index=False, na_rep='')
        print('Done, check file {}!'.format(opf))


def main():
    actions = (
        ('correct', 'correct wrong genotype calls'),
        ('cleanup', 'clean redundant info in the tmp matirx file'),
        ('format', 'convert genotype matix file to other formats'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == '__main__':
    main()
