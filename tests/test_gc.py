import os.path as op
import pandas as pd
import matplotlib.pyplot as plt

from schnablelab.GC.corrector import Correction

data_dir = op.abspath(op.dirname(__file__))+'/../src/schnablelab/GC/data'
df0 = pd.read_csv(f'{data_dir}/F2_simulated.map', sep='\t')
df1 = pd.read_csv(f'{data_dir}/F2_X0.45miss0.35.map', sep='\t')

win_size = 15
error_a = 0.03
error_b = 0.01
error_h = abs(error_a - error_b)
chr_mapping = {'A': 0, 'B': 2, 'X': 1, '-': 9}

# ground truth seq
b = df0.loc[df1['chr'] == 'chr1', 'SM-1'].to_list()
h = df0.loc[df1['chr'] == 'chr1', 'SM-2'].to_list()
a = df0.loc[df1['chr'] == 'chr1', 'SM-3'].to_list()
b_num = pd.Series([*map(chr_mapping.get, b)])
h_num = pd.Series([*map(chr_mapping.get, h)])
a_num = pd.Series([*map(chr_mapping.get, a)])

# before correction seq
before_b = df1.loc[df1['chr'] == 'chr1', 'SM-1'].to_list()
before_h = df1.loc[df1['chr'] == 'chr1', 'SM-2'].to_list()
before_a = df1.loc[df1['chr'] == 'chr1', 'SM-3'].to_list()
before_b_num = pd.Series([*map(chr_mapping.get, before_b)])
before_h_num = pd.Series([*map(chr_mapping.get, before_h)])
before_a_num = pd.Series([*map(chr_mapping.get, before_a)])


def run_test(before_seq, groundtruth, round_n=3):
    before_acc = (before_seq != groundtruth).sum()/len(groundtruth)
    print(f'before correction error: {before_acc:.2f}')
    seq = before_seq
    round = 0
    while round < round_n:
        correct = Correction(win_size, error_a, error_h, error_b, seq)
        seq = correct.get_corrected()
        after_acc = (seq != groundtruth).sum()/len(groundtruth)
        print(f'correction error after round {round}: {after_acc:.2f}')
        round += 1
    return seq, before_acc, after_acc


def make_plot(gt_seq, before_seq, after_seq,
              before_acc, after_acc, fn,
              ):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 3))
    ax1.scatter(range(len(gt_seq)), gt_seq, s=0.5)
    ax1.set_title('ground truth')
    ax1.set_ylim(-0.5, 9.5)

    ax2.scatter(range(len(before_seq)), before_seq, s=0.5)
    ax2.set_title(f'before correction ({before_acc:.2f})')
    ax2.set_ylim(-0.5, 9.5)

    ax3.scatter(range(len(after_seq)), after_seq, s=0.5)
    ax3.set_title(f'after correction ({after_acc:.2f})')
    ax3.set_ylim(-0.5, 9.5)
    plt.tight_layout()
    plt.savefig(f'debug_{fn}.png', dpi=200)


def main():
    for seq_type, before_seq, gt_seq in zip(
            ['homo_b', 'homo_a', 'hete_h'],
            [before_b_num, before_h_num, before_a_num],
            [b_num, h_num, a_num]):
        print(f'test sequence type: {seq_type}')
        after_seq, before_acc, after_acc = run_test(before_seq, gt_seq)
        if after_acc > 0.05:
            make_plot(gt_seq, before_seq, after_seq, before_acc, after_acc,
                      seq_type)
            raise Exception(f'correction accuracy for {seq_type} task > 5%, check artifacts for details!')


if __name__ == "__main__":
    main()
