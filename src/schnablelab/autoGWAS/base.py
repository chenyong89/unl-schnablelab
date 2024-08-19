import re
import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import Counter

geno_one2two = {
    "A": "AA",
    "C": "CC",
    "G": "GG",
    "T": "TT",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "N": "NN",
    "-": "--",
}

geno_two2one = {
    "AA": "A",
    "CC": "C",
    "GG": "G",
    "TT": "T",
    "GA": "R",
    "AG": "R",
    "TC": "Y",
    "CT": "Y",
    "CG": "S",
    "GC": "S",
    "TA": "W",
    "AT": "W",
    "TG": "K",
    "GT": "K",
    "CA": "M",
    "AC": "M",
    "NN": "N",
    "--": "-",
}


def geno_to_value(ref, alt, genosList):
    """
    convert ref homo, alt homo, and hetero genotype to 0, 2, and 1 respectively
    """
    newlist = []
    for k in genosList:
        if len(set(k)) == 1 and k[0] == ref:
            newlist.append("0")
        elif len(set(k)) == 1 and k[0] == alt:
            newlist.append("2")
        elif len(set(k)) == 2:
            newlist.append("1")
        else:
            raise Exception("genotype error !")
    return newlist


def num_from_chr(chr_name):
    """
    extract number from a corn/sorghum chromosome name
    """
    chr_num = re.findall(r"\d+$", chr_name)
    if len(chr_num) == 1:
        return int(chr_num[0])
    else:
        raise Exception(f"couldn't identify chromosome name: '{chr_name}' ")


class ParseHapmap:
    """
    parse file in hapmap format
    """

    def __init__(self, filename):
        """
        args:
            filename: filename in hapmap format
        """
        self.filename = filename
        with open(filename) as f:
            header = f.readline()
            attributes = header.split()[:11]
            samples = header.split()[11:]
            num_sample = len(samples)
            firstgeno = f.readline().split()[11]
            code_type = "single" if len(firstgeno) == 1 else "double"
            num_snp = sum(1 for _ in f)
            dtype_dict = {i: "str" for i in header.split()}
            dtype_dict["pos"] = np.int64
        self.header = header
        self.attributes = attributes
        self.samples = samples
        self.num_sample = num_sample
        self.num_snp = num_snp + 1
        self.code_type = code_type
        self.dtype_dict = dtype_dict

    def to_df(self, sorting=False):
        """
        convert hapmap file to pandas dataframe
        args:
            sort: if sorting the dataframe
        """
        df = pd.read_csv(
            self.filename, delim_whitespace=True, dtype=self.dtype_dict
        )
        if sorting:
            chrs = list(df["chrom"].unique())
            chrs_ordered = sorted(chrs, key=num_from_chr)
            df["chrom"] = pd.Categorical(
                df["chrom"], chrs_ordered, ordered=True
            )
            df = df.sort_values(["chrom", "pos"]).reset_index(drop=True)
        if self.code_type == "single":
            print("converting the single type to double type...")
            df_part2 = df.iloc[:, 11:].applymap(geno_one2two.get).fillna("NN")
            df = pd.concat([df.iloc[:, :11], df_part2], axis=1)
            self.type = "double"
        return df

    def to_map_ped(self, df, missing=True):
        """
        convert pandas dataframe to PLINK map and ped files
        args:
            df: pandas dataframe
            missing: if any missing gneotype data in the table
        """
        df_map = df[["rs#", "chrom", "pos"]]
        df_map["centi"] = 0
        map_cols = ["chrom", "rs#", "centi", "pos"]
        df_map = df_map[map_cols]

        if missing:
            df = df.replace("NN", "00")
        df_ped = np.zeros((self.num_sample, 6 + self.num_snp * 2), dtype="str")
        pbar = tqdm(self.samples)
        for idx, col in enumerate(pbar):
            col_a1 = df[col].str[0].to_numpy()
            col_a2 = df[col].str[1].to_numpy()
            df_ped[idx, 6:] = np.column_stack((col_a1, col_a2)).ravel()
            pbar.set_description(f"converting {col}...")
        df_ped = np.where(df_ped == "", 0, df_ped)
        df_ped = pd.DataFrame(df_ped, dtype="str")
        df_ped[1] = self.samples
        return df_map, df_ped

    def count_genotype(self, geno_list, a1, a2):
        """
        count homo_ref, homo_alt, and heteroygous genotypes in a POS across
        all samples

        args:
            geno_list (list): list of all genotypes
            a1 (str): character of allele 1
            a2 (str): character of allele2
        """
        if self.code_type == "single":
            h = geno_two2one[a1 + a2]
            counts = Counter(geno_list)
            num_homo1, num_homo2, num_hete = counts[a1], counts[a2], counts[h]
        if self.code_type == "double":
            a, b = a1 * 2, a2 * 2
            h1, h2 = a1 + a2, a2 + a1
            counts = Counter(geno_list)
            num_homo1, num_homo2 = counts[a], counts[b]
            num_hete = counts[h1] + counts[h2]
        return num_homo1, num_homo2, num_hete

    @property
    def missing_rate(self):
        """
        yield (line, missing rate) for each line
        """
        with open(self.filename) as f:
            next(f)
            for i in f:
                c = Counter(i.split()[11:])
                num_missing = c["NN"] + c["N"]
                yield i, num_missing / self.num_sample

    @property
    def maf(self):
        """
        calculate Minor Allele Frequency for each line the table

        yield (line, maf) for each line
        """
        with open(self.filename) as f:
            next(f)
            for i in f:
                j = i.split()
                alleles = j[1].split("/")
                try:
                    a1, a2 = alleles
                except ValueError:
                    yield i, 0
                else:
                    num_homo1, num_homo2, num_hete = self.count_genotype(
                        j[11:], a1, a2
                    )
                    num_a1, num_a2 = (
                        num_homo1 * 2 + num_hete,
                        num_homo2 * 2 + num_hete,
                    )
                    try:
                        maf = min(num_a1, num_a2) / (num_a1 + num_a2)
                    except ZeroDivisionError:
                        yield i, 0
                    else:
                        yield i, maf

    @property
    def heteros(self):
        """
        calculate heterozygous reate for each line the table

        yield (line, heterozgous rate) for each line
        """
        with open(self.filename) as f:
            next(f)
            for i in f:
                j = i.split()
                alleles = j[1].split("/")
                try:
                    a1, a2 = alleles
                except ValueError:
                    yield i, 1
                else:
                    num_homo1, num_homo2, num_hete = self.count_genotype(
                        j[11:], a1, a2
                    )
                    if num_hete > max(num_homo1, num_homo2):
                        yield i, 1
                    else:
                        hr = num_hete / float(num_homo1 + num_homo2 + num_hete)
                        yield i, hr


class ParseGWASresults:
    def __init__(self, filename, software, sorting=False, idx_cols=None):
        """
        Parse output tables generated from various GWAS tools

        Args:
            filename: filename of GWAS output table
            software (str): GWAS tools (gemma, farmcpu, mvp, gapit)
            sorting (bool): if sorting gwas output table
            idx_cols (list): specify the indices of which columns to use if
                using other GWAS tools
        """
        self.filename = filename
        self.software = software
        self.sorting = sorting
        self.idx_cols = idx_cols

        if self.software == "gemma":
            dtype_dict = {
                "chr": "str",
                "rs": "str",
                "ps": np.int64,
                "p_lrt": np.float64,
            }
            df = pd.read_csv(
                self.filename,
                delim_whitespace=True,
                usecols=["chr", "rs", "ps", "p_lrt"],
                dtype=dtype_dict,
            )
            df = df[["rs", "chr", "ps", "p_lrt"]]
        if self.software == "farmcpu":
            dtype_dict = {
                "Chromosome": "str",
                "SNP": "str",
                "Position": np.int64,
                "P.value": np.float64,
            }
            df = pd.read_csv(
                self.filename,
                usecols=["SNP", "Chromosome", "Position", "P.value"],
                dtype=dtype_dict,
            )
        if self.software == "gapit":
            dtype_dict = {
                "Chromosome": "str",
                "SNP": "str",
                "Position": np.int64,
                "P.value": np.float64,
            }
            df = pd.read_csv(
                self.filename,
                usecols=["SNP", "Chromosome", "Position ", "P.value"],
            )
        if self.software == "other":
            if self.idx_cols is None:
                raise Exception(
                    "specify which columns to use if choosing other gwas tools!"
                )
            if len(self.idx_cols) != 4:
                raise Exception("usecols must have the lenght of 4")
            with open(self.filename) as _:
                j = list(pd.read_csv(self.filename).columns)
                snp_idx, chr_idx, pos_idx, pv_idx = self.idx_cols
                snp_col = j[snp_idx]
                chr_col = j[chr_idx]
                pos_col = j[pos_idx]
                pv_col = j[pv_idx]
                dtype_dict = {
                    snp_col: "str",
                    chr_col: "str",
                    pos_col: np.int64,
                    pv_col: np.float64,
                }
            cols = [snp_col, chr_col, pos_col, pv_col]
            df = pd.read_csv(self.filename, usecols=cols, dtype=dtype_dict)[
                cols
            ]
        df.columns = ["snp", "chr", "pos", "pvalue"]
        df["pvalue"] = -np.log10(df["pvalue"])
        df.columns = ["snp", "chr", "pos", "-log10Pvalue"]
        if sorting:
            chrs = list(df["chr"].unique())
            chrs_ordered = sorted(chrs, key=num_from_chr)
            df["chr"] = pd.Categorical(df["chr"], chrs_ordered, ordered=True)
            df = df.sort_values(["chr", "pos"]).reset_index(drop=True)
        self.df = df
        self.num_snp = df.shape[0]

    def SignificantSNPs(self, p_level=0.05, MeRatio=1, sig_cutoff=None):
        """
        extract Significant SNPs

        sig_cutoff: the threshold of log10 transformed p values
        """
        cutoff = (
            -np.log10(p_level / (MeRatio * self.num_snp))
            if sig_cutoff is None
            else sig_cutoff
        )
        df_sigs = self.df[self.df["-log10Pvalue"] >= cutoff].reset_index(
            drop=True
        )
        return df_sigs


def cal_MAF_row(hmp_row):
    j = hmp_row.split()
    allele1, allele2 = j[1].split("/")
    genos = "".join(j[11:])
    a1, a2 = genos.count(allele1), genos.count(allele2)
    maf = min(a1, a2) / (a1 + a2)
    count = len(genos) * maf

    (minor_allele, major_allele) = (
        (allele1, allele2) if a1 <= a2 else (allele2, allele1)
    )
    minor_idx, major_idx, hetero_idx = [], [], []
    for m, n in enumerate(j[11:]):
        k = list(set(n))
        if len(k) == 1:
            if k[0] == minor_allele:
                minor_idx.append(m + 11)
            elif k[0] == major_allele:
                major_idx.append(m + 11)
            else:
                print(n)
                print("bad allele!!!")
        else:
            hetero_idx.append(m + 11)

    return j[0], maf, count, minor_idx, major_idx, hetero_idx


FarmCPU_header = """
library("bigmemory")
library("biganalytics")
library("compiler")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")
setwd(".")
myY <- read.table("{pheno}", head = TRUE)
myGM <- read.table("{geno_prefix}.GM", head = TRUE)
myGD <- read.big.matrix("{geno_prefix}.GD", type="char", sep="\t", head = TRUE)
myCV <- read.table("{pca}", head = TRUE)
#Step 2: Run FarmCPU
myFarmCPU <- FarmCPU(Y=myY, GD=myGD, GM=myGM, CV=myCV,
    method.bin="optimum", bin.size=c(5e5, 5e6, 5e7),
    bin.selection=seq(10, 100, 10),
    threshold.output=1,
    memo='{memo}')
"""

GAPIT_header = """library(multtest)
library(gplots)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")
setwd(".")
myY <- read.table("{pheno}", head=TRUE) # sep is the space related separator
myGM <- read.table("{geno_prefix}.GM", head=TRUE)
myGD <- read.table("{geno_prefix}.GD", head=TRUE)
myCV <- read.table("{pca}", head=TRUE)
myKI <- read.table("{kinship}", head=FALSE)
#Step 2: Run GAPIT
myGAPIT <- GAPIT(Y=myY, GD=myGD, GM=myGM, CV=myCV, KI=myKI, memo='{memo}')
"""

MVP_Data_header = """
library(MVP)
MVP.Data(fileHMP="{hapmap}", sep.hmp="\t", sep.phe="\t", SNP.effect="Add",
         fileKin=TRUE, filePC=TRUE, out="{output}", priority="speed",)
"""

MVP_Run_header = """
library(MVP)
phenotype <- read.table("{pheno}", head=TRUE)
genotype <- attach.big.matrix("{output_prefix}.geno.desc")
map <- read.table("{output_prefix}.map", head = TRUE)
Kinship <- attach.big.matrix("{output_prefix}.kin.desc")
Covariates <- attach.big.matrix("{output_prefix}.pc.desc")
imMVP <- MVP(phe=phenotype, geno=genotype, map=map, K=Kinship,
             CV.MLM=Covariates, CV.FarmCPU=Covariates, maxLoop=10,
             method.bin="FaST-LMM", priority="speed", threshold=0.05,
             method=c("MLM", "FarmCPU"))
"""

MSTMap_header = """
population_type RIL%s
population_name original
distance_function kosambi
cut_off_p_value 0.000001
no_map_dist 15.0
no_map_size 0
missing_threshold 0.3
estimation_before_clustering yes
detect_bad_data yes
objective_function COUNT
number_of_loci %s
number_of_individual %s

"""
