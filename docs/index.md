# Welcome to the `schnablelab` package 

This package consists of codes developed by Chenyong Miao when he was a PhD student in the Schnable Lab at the University of Nebraska-Lincoln from 2016-2020. The majority of code was used in Chenyong's daily research on high throughput phenotyping, genotyping, as well as GWAS (genome-wide association studies). Code has been wrapped up in command line tools to faciliate analyses on local computers and HCC (Holland Computing Center) campus server. The package is organized into different modules and each module is composed by couple of scripts which contains multiple actions/functions to conduct specific tasks. The general usage of command line follows this pattern:

```
$ python -m schnablelab.[MODULE].[SCRIPT] [ACTION]
```

To use this pakcage, please follow the [installation](../docs/installation.md) page to install the package first. 

### 1. autoGWAS module

This module helps preparing data, conduting GWAS using various statistical models, and analyzing GWAS results

```
$ python -m schnablelab.autoGWAS
Usage:
    python -m schnablelab.autoGWAS.SCRIPT


Available SCRIPTs:
    data_pipeline | Data QA/QC before conducting GWAS
             gwas | Conduct GWAS using various algorithms
        post_gwas | Post-GWAS analyses
```

### 2. genotyping module

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

### 3. phenotyping module

This module is developed to process high throughtput phenotyping data generated from UNL Greenhouse Innovation Center. 

```
$ python -m schnablelab.phenotyping
Usage:
    python -m schnablelab.phenotyping.SCRIPT


Available SCRIPTs:
       extract_traits | Extract sorghum inflorescence width and height, stem height, and plant height
             htp_data | Utils to deal with high throughput phenotyping data
    plot_growth_curve | Plot time-series traits extracted from sorghum RGB images
```

### 4. zooniverse module

This module only has one script `zookeeper` which helps uploading and donwloading data to/from the zooniverse croud sourcing platform. 

```
$ python -m schnablelab.zooniverse.zookeeper
Usage:
    python -m schnablelab.zooniverse.zookeeper ACTION


Available ACTIONs:
             toy | Random pick up some images for testing purporse
          divide | Divide a large number of images to sevearl subsets
          upload | Load images to zooniverse
    batch_upload | Upload multiple dirs on HCC
          export | Get annotation and other exports
        manifest | Generate a manifest for zooniverse subject set upload
```

### 5. GC (Genotype-Corrector) module

GC is a bioinformatics tool to impuate missing data and correct wrong genotype calls in bi-parental population. You can find more details about this tool from our [publication](https://doi.org/10.1038/s41598-018-28294-0). Please follow more detailed [tutorials](../docs/gc_tutorial.md) to learn how to conduct GC on your own data. 

```
$ python -m schnablelab.GC
Usage:
    python -m schnablelab.GC.SCRIPT


Available SCRIPTs:
      data_qc | Perform data QC before running Genotype-Corrector
    corrector | Correct Genotype Calls in biparental populations
```

### 6. hcc module

hcc module provides convienet tools when conducting analyses on HCC 

```
$ python -m schnablelab.hcc
Usage:
    python -m schnablelab.hcc.SCRIPT


Available SCRIPTs:
         job | Create, submit, canceal jobs
     request | Request a node
    traverse | Traverse files to avoid purge policy on hcc
```