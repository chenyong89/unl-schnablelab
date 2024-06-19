"""
Base utilties for genotype corrector
"""
import sys
import numpy as np
import pandas as pd
from collections import defaultdict
from configparser import ConfigParser


class ParseConfig(object):
    """
    Parse the configure file using configparser
    """
    def __init__(self, configfile):
        config = ConfigParser()
        assert config.read(configfile), 'config file does not exist!'
        self.po_type = config.get('Section1', 'Population_type')
        self.gt_a = config.get('Section2', 'Letter_for_homo1')
        self.gt_h = config.get('Section2', 'Letter_for_hete')
        self.gt_b = config.get('Section2', 'Letter_for_homo2')
        self.gt_miss = config.get('Section2', 'Letter_for_missing_data')
        self.error_a = config.getfloat('Section3', 'error_rate_for_homo1')
        self.error_b = config.getfloat('Section3', 'error_rate_for_homo2')
        self.error_h = abs(self.error_a - self.error_b)
        self.win_size = config.getint('Section4', 'Sliding_window_size')


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def getChunk(fn, ignore=1):
    '''
    ignore: rows starts with pound sign
    '''
    df0_chr = defaultdict(int)
    chr_order = []
    with open(fn) as f:
        for dash_line in range(ignore):
            f.readline()
        for i in f:
            j = i.split()[0].split('-')[0]
            df0_chr[j] += 1
            if j in chr_order:
                pass
            else:
                chr_order.append(j)
    if len(chr_order) != len(set(chr_order)):
        sys.exit('Please check your marker name and sort them by chr name.')
    return chr_order, df0_chr


def get_blocks(np_1d_array, dist=150, block_size=2):
    """
    group values to a block with specified distance
    Examples:
    >>> a = np.array([1,2,4,10,12,13,15])
    >>> test(a, dist=1)
    [[1, 2], [12, 13]]
    >>> test(a, dist=2)
    [[1, 2, 4], [10, 12, 13, 15]]
    """
    first_val = np_1d_array[0]
    temp = [first_val]  # save temp blocks
    pre_val = first_val
    results = []
    for val in np_1d_array[1:]:
        if (val - pre_val) <= dist:
            temp.append(val)
        else:
            if len(temp) >= block_size:
                results.append(temp)
            temp = [val]
        pre_val = val
    if len(temp) >= block_size:
        results.append(temp)
    return results


def sort_merge_sort(arrays):
    """
    get redundant lists by merging lists with overlaping region.
    Example:
    >>> a = [[1,3], [3, 5], [6,10], [7, 9], [11,15], [11,12],[16,30]]
    >>> sort_merge_sort(a)
    >>> [array([1, 3, 5]), array([ 6,  7,  9, 10]), array([11, 12, 15]), [16, 30]]
    """
    val_start = [i[0] for i in arrays]
    val_end = [i[-1] for i in arrays]
    df = pd.DataFrame(dict(zip(
        ['array', 'val_start', 'val_end'],
        [arrays, val_start, val_end]))).sort_values(['val_start', 'val_end'])\
        .reset_index(drop=True)
    first_arr = df.loc[0, 'array']
    temp = first_arr
    pre_arr = first_arr
    results = []
    for arr in df.loc[1:, 'array']:
        if arr[0] <= pre_arr[-1]:
            temp.extend(arr)
        else:
            if len(temp) == len(pre_arr):
                results.append(pre_arr)
            else:
                temp_sorted_unique = pd.Series(temp).sort_values().unique()
                results.append(temp_sorted_unique)
            temp = arr
        pre_arr = arr
    if len(temp) == len(pre_arr):
        results.append(pre_arr)
    else:
        temp_sorted_unique = pd.Series(temp).sort_values().unique()
        results.append(temp_sorted_unique)
    return results


def bin_markers(df, diff=0, missing_value='-'):
    """
    merge consecutive markers with same genotypes
    return slelected row index
    """
    df = df.replace(missing_value, np.nan)
    first_row = df.iloc[0, :]
    temp = [df.index[0]]  # save temp index
    pre_row = first_row
    df_rest = df.iloc[1:, :]
    result_ids = []
    for idx, row in df_rest.iterrows():
        df_tmp = pd.concat([pre_row, row], axis=1).dropna()
        diff_num = (df_tmp.iloc[:, 0] != df_tmp.iloc[:, 1]).sum()
        if diff_num <= diff:  # get ride of it
            temp.append(idx)
        else:
            result_ids.append(temp)
            temp = [idx]
        pre_row = row
    result_ids.append(temp)

    results = []
    represent_idx, block_idx = [], []
    for index in result_ids:
        if len(index) > 1:
            df_tmp = df.loc[index, :]
            good_idx = df_tmp.isnull().sum(axis=1).idxmin()
            results.append(good_idx)
            represent_idx.append(good_idx)
            block_idx.append(index)
        else:
            results.append(index[0])
    return represent_idx, block_idx, results
