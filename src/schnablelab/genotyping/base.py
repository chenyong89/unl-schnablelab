"""
base class and functions to handle with vcf file
"""

import re
import sys
import pandas as pd
from collections import Counter


def find_sm(target_str, re_pattern):
    sms = re_pattern.findall(target_str)
    if len(sms) == 1:
        sm = sms[0][1:-1]
        return "-".join(re.split("[_-]", sm))
    else:
        sys.exit(f"bad file name '{target_str}'!")


class ParseVCF:
    """
    parse vcf file
    """

    def __init__(self, filename):
        self.fn = filename
        with open(filename) as f:
            n = 0
            hash_chunk = []
            hash_chunk2 = []
            num_SNPs = 0
            for i in f:
                if i.startswith("##"):
                    n += 1
                    hash_chunk.append(i)
                    hash_chunk2.append(i)
                    continue
                if i.startswith("#"):
                    SMs_header = i.split()[:9]
                    SMs = i.split()[9:]
                    numSMs = len(SMs)
                    n += 1
                    hash_chunk.append(i)
                else:
                    num_SNPs += 1
        self.num_SNPs = num_SNPs
        self.SMs_header = SMs_header
        self.SMs = SMs
        self.numSMs = numSMs
        self.numHash = n
        self.HashChunk = hash_chunk
        self.HashChunk2 = hash_chunk2
        self.numHeaderLines = len(self.HashChunk)
        self.hmpfield11 = "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode"
        self.hmpheader = self.hmpfield11 + "\t" + "\t".join(self.SMs) + "\n"

    def to_df(self):
        df = pd.read_csv(
            self.fn, skiprows=range(self.numHash - 1), delim_whitespace=True
        )
        return df

    @property
    def to_hmp(self):
        """
        yield line in hmp format
        """
        cen_NA = "+\tNA\tNA\tNA\tNA\tNA\tNA"
        with open(self.fn) as f:
            for _ in range(self.numHash):
                next(f)
            for i in f:
                j = i.split()
                a1, a2 = j[3], j[4]
                if len(a1) == len(a2) == 1:
                    a1a2 = "".join([a1, a2])
                    a2a1 = a1a2[::-1]
                    a1a1, a2a2 = a1 * 2, a2 * 2
                elif len(a1) == 1 and len(a2) > 1:
                    a1a2 = "".join(["-", a2[-1]])
                    a2a1 = a1a2[::-1]
                    a1a1, a2a2 = "--", a2[-1] * 2
                elif len(a1) > 1 and len(a2) == 1:
                    a1a2 = "".join([a1[-1], "-"])
                    a2a1 = a1a2[::-1]
                    a1a1, a2a2 = a1[-1] * 2, "--"
                else:
                    print("bad format line:\n", j[2])
                    continue
                geno_dict = {
                    "0/0": a1a1,
                    "0|0": a1a1,
                    "0/1": a1a2,
                    "0|1": a1a2,
                    "1/0": a2a1,
                    "1|0": a2a1,
                    "1/1": a2a2,
                    "1|1": a2a2,
                    "./.": "NN",
                    ".|.": "NN",
                }
                genos = list(map(geno_dict.get, j[9:]))
                if None in genos:
                    print(i)
                    sys.exit("unknow genotype detected!")
                genos = "\t".join(genos)
                rs, chr, pos = j[2], j[0], j[1]
                alleles = "/".join([a1a2[0], a1a2[1]])
                new_line = (
                    "\t".join([rs, alleles, chr, pos, cen_NA, genos]) + "\n"
                )
                yield new_line

    @property
    def missing_rate(self):
        """
        calculate missing rate
        """
        with open(self.fn) as f:
            for _ in range(self.numHash):
                next(f)
            for i in f:
                num_miss = i.count("./.") + i.count(".|.")
                yield i, num_miss / self.numSMs

    @staticmethod
    def count_geno(geno_ls):
        """
        count each genotypes
        """
        c = Counter(geno_ls)
        num_a = c["0/0"] + c["0|0"]
        num_b = c["1/1"] + c["1|1"]
        num_h = c["0/1"] + c["1/0"] + c["0|1"] + c["1|0"]
        return num_a, num_b, num_h

    @property
    def hetero_rate(self):
        """
        yield (line, heterozygous rate) for each line
        """
        with open(self.fn) as f:
            for _ in range(self.numHash):
                next(f)
            for i in f:
                num_a, num_b, num_h = self.count_geno(i.split()[9:])
                if num_h > max(num_a, num_b):
                    yield i, 1
                else:
                    yield i, num_h / float(num_a + num_b + num_h)

    @property
    def maf(self):
        """
        yield minor allele frequence for each line
        """
        with open(self.fn) as f:
            for _ in range(self.numHash):
                next(f)
            for i in f:
                num_a, num_b, num_h = self.count_geno(i.split()[9:])
                a1, a2 = num_a * 2 + num_h, num_b * 2 + num_h
                yield i, min(a1, a2) / (a1 + a2)
