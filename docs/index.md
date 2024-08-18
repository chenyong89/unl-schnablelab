# Welcome to the `schnablelab` Package

The schnablelab package comprises code developed by Chenyong Miao during his PhD at the Schnable Lab, University of Nebraska-Lincoln, spanning from 2016 to 2020. This code primarily supported Chenyong's research activities in high-throughput phenotyping/genotyping, and genome-wide association studies (GWAS). The functionality has been encapsulated into command-line tools to facilitate analyses on both local computers and the HCC (Holland Computing Center) campus server. The package is structured into distinct modules, each containing several scripts that encompass multiple functions tailored for specific tasks. The general command-line usage adheres to the following pattern:

```
$ python -m schnablelab.[MODULE].[SCRIPT] [ACTION]
```

To tuilize this pakcage, please follow the [installation](installation.md) page for initial setup instructions.

## autoGWAS module

This module facilitates data preparation, conducts GWAS using various statistical models, and analyzes GWAS results.

```
$ python -m schnablelab.autoGWAS
Usage:
    python -m schnablelab.autoGWAS.SCRIPT


Available SCRIPTs:
    data_pipeline | Data QA/QC before conducting GWAS
             gwas | Conduct GWAS using various algorithms
        post_gwas | Post-GWAS analyses
```

## genotyping module

This module contains SNP calling pipelines

```
$ python -m schnablelab.genotyping
Usage:
    python -m schnablelab.genotyping.SCRIPT


Available SCRIPTs:
     pre_snp_calling | Get files ready for SNP calling
         snp_calling | Call SNPs on high throughput sequencing data using GATK, Freebayes, SAMtools
    post_snp_calling | Post SNP calling functions
```

## phenotyping module

This module processes high-throughput phenotyping data generated from the UNL Greenhouse Innovation Center.

```
$ python -m schnablelab.phenotyping
Usage:
    python -m schnablelab.phenotyping.SCRIPT


Available SCRIPTs:
       extract_traits | Extract sorghum inflorescence width and height, stem height, and plant height
             htp_data | Utils to deal with high throughput phenotyping data
    plot_growth_curve | Plot time-series traits extracted from sorghum RGB images
```

## zooniverse module

This module features a single script, `zookeeper`, designed to upload and download data to/from the Zooniverse crowdsourcing platform.

```
$ python -m schnablelab.zooniverse.zookeeper
Usage:
    python -m schnablelab.zooniverse.zookeeper ACTION


Available ACTIONs:
             toy | Random pick up some images for testing purporse
          divide | Divide a large number of images to several subsets
          upload | Load images to zooniverse
    batch_upload | Upload multiple dirs on HCC to zooniverse
          export | Get annotation and other exports
        manifest | Generate a manifest for zooniverse subject set upload
```

## GC (Genotype-Corrector) module

GC is a bioinformatics tool used to impute missing data and correct erroneous genotype calls in bi-parental populations. For more details, refer to our [publication](https://doi.org/10.1038/s41598-018-28294-0). Detailed [tutorials](gc_tutorial.md) on conducting GC with your own data can be found here.

```
$ python -m schnablelab.GC
Usage:
    python -m schnablelab.GC.SCRIPT


Available SCRIPTs:
      data_qc | Perform data QC before running Genotype-Corrector
    corrector | Correct Genotype Calls in biparental populations
```

## hcc module

The hcc module provides convenient tools for analyses conducted on the HCC.

```
$ python -m schnablelab.hcc
Usage:
    python -m schnablelab.hcc.SCRIPT


Available SCRIPTs:
         job | Create, submit, canceal jobs
     request | Request a node
    traverse | Traverse files to avoid purge policy on hcc
```

---

Please note that most of the code in this repository was developed over four years ago. Some scripts, particularly those related to HCC, may no longer be compatible with the updates to the campus server since my graduation in 2020. I apologize for any outdated or inefficient code that may still be present in the package.

- Chenyong Miao
- 06/23/2024
