# OutcrossSeq

OutcrossSeq is a highly scalable genotyping software for identifying genomic variation,genotyping, haplotype phasing and mapping agronomically important loci in the populations of outcrossing crops, based on whole-genome low-coverage resequencing data. 

## Availability
You can download the source code and the Demo data here:  
[https://github.com/xhhuanglab/OutcrossSeq](https://github.com/xhhuanglab/OutcrossSeq) 
You can download the Docker container here:
[https://hub.docker.com/repository/docker/mjchenjojo/outcrossseq](https://hub.docker.com/repository/docker/mjchenjojo/outcrossseq)
The ditels of Docker container can be found in [supplementary note](https://github.com/xhhuanglab/OutcrossSeq/blob/master/Supplemental%20note.md).
``` 
docker pull mjchenjojo/outcrossseq:v1.0 
docker load --input mjchenjojo/outcrossseq:v1.0
```
## Dependencies
The code is implemented in perl (v5.16.3) and can be run from linux shell. To run OutcrossSeq the following libraries and tools should be available:
1. [R](https://stat.ethz.ch/pipermail/r-announce/2017/000616.html) version 3.4.1
2. [Statistics::Basic](https://metacpan.org/pod/distribution/Statistics-Basic/lib/Statistics/Basic.pod) module
3. [Number::Format](https://metacpan.org/pod/Number::Format) modular


## Software structure
There are three modules in OutcrossSeq.
- "Diploid-Outcrossing" module designed for diploid outcrossing F1 population
- "Double-Cross" module designed for double cross population
- "Autopolyploid Plant" module designed for autopolyploid F1 population. 

## Data preparation

The input file of "Diploid-Outcrossing" module is variant calling file of two heterozygous parents and all F1 individuals.Normally in VCF format.
See Demo data: `demo_diploid.vcf`

The input file of "Double-Cross" module is variant calling file of four inbred parents and all offspring individuals.
See Demo data:`demo_double.vcf`

The input files of "Autopolyploid Plant" module are two separated variant calling files of parents and offsprings, since for parents we call variants in polyploid mode, meanwhile, for offsprings we call variants in diploid mode.
See Demo data for two hexaploid parents: `demo-Parents.vcf` and their F1 offsprings `demo_autopolyploid.vcf`



- **For variant calling，we recommend to generate VCF files for each chromosome which is more convenient to run the code and to tune parameters in following steps.**


## Getting Started
You can choose to run in steps or run together.
The detail of run in steps, please see [Use cases](#cases).
Parameter information of run together, please see [supplementary note](https://github.com/xhhuanglab/OutcrossSeq/blob/master/Supplemental%20note.md).
```	
perl system_diploid_cross.pl 5 5 100 0.6 0.3 0.35 60 1000 100000 parent1 parent2 demo_diploid demo
#"Diploid-Outcrossing" module
perl system_double_cross.pl 16 17 18 19 5 5 100 100000 1 6 demo_double
#"Double-Cross" module
perl system_autopolyploid.pl 5 5 0.08 1 50 0.9 3 100 Parent1 Parent2 demo_autopolyploid-Parents demo_autopolyploid
#"Autopolyploid Plant" module
```
## Table of Contents
- [Users' Guide](#uguide)
  - [Installation](#Installation)
  - [Introduction](#Introduction)
  - [Use cases](#cases)
    - [Diploid-Outcrossing](#Diploid)
    - [Double-Cross](#Double)
    - [Autopolyploid Plant](#Autopolyploid)
  - [Getting help](#help)
  - [Citing OutcrossSeq](#citing)


## <a name="uguide"></a>Users' Guide
### <a name="Installation"></a>Installation
```
git clone https://github.com/xhhuanglab/OutcrossSeq.git
gunzip OutcrossSeq.gz
cd OutcrossSeq
#list three modules
ls 
#autopolyploid_crop  diploid_outcrossing_species  double_cross_populations  readme.md
```
### <a name="Introduction"></a>Introduction
OutcrossSeq is a method and pipeline for mapping agronomically important loci in outcrossing crop based on whole-genome low-coverage resequencing (which normally means low-cost) of large genetic population. It has three modules: "Diploid-Outcrossing","Double-Cross" and "Autopolyploid Plant".Typical OutcrossSeq pipeline includes three steps – sequencing, genotyping and imputation.

Step1 sequencing: Using Tn5 enzymes for low-cost library construction and performed population-scale multiplex sequencing. Hundreds of individuals are pooled within the same library (each individual with its unique barcode). Hundreds of individuals pooled within the same library (each individual with its unique barcode) were sequenced simultaneously and further demultiplexed according to the index tag, making the sequencing step rapidly and cost-efficiently, comparing with sequencing individuals one by one.

Step2 genotyping: Variant calling sofware, such as GATK, can be used for genotype calling from sequencing data. Due to the low coverage feature of individuals, we recommend to rule out problematic variants prior to genotype imputation. Here we can rely on the parental data if they holding higer coverage or well determind genotype.Following suggestions could be applied depend on which module will be practiced.

For "Diploid-Outcrossing" module: Only the SNP sites that are heterozygous in at least one parental line are selected for next step.

For "Double-Cross" module: SNP sites that segregated in four haplotypes satisfied 1:3 ratio, e.g., A, T, A and A for H1, H2, H3 and H4, respectively; or satisfied 2:2 ratio, e.g., A, T, A and T for H1, H2, H3 and H4, respectively. In other word, we only count in these biallelic SNP sites.

For “Autopolyploid Plant” module: Only simplex and double-simplex SNP sites in parents are counted in.


Step3 imputation: possibly the most important step of OutcrossSeq, we have different strategies for our three modules.

For "Diploid-Outcrossing" module: Base on the kinship calculation in 100kb window, we fill up the missing SNPs rely on known four genotypes to complete the bin matrix, we presume few to no recombination has happened within such window in most individuals.

For "Double-Cross" module: Base on the proportion of SNP in 100kb window to fill up the missing SNPs rely on known genotypes to get the complete bin matrix.

For “Autopolyploid Plant” module: Firstly, we divide bins according to the number of SNPs, let's say, 500 SNPs. Secondly, we cluster the SNPs according to the correlation coefficient between them. Lastly, fill up unknown sites based on the known information of each cluster to obtain a SNP matrix.

### <a name="cases"></a>Use cases
OutcrossSeq works with VCF format files as input. After imputation, the output file is a genotype file.

"Diploid-Outcrossing" and "Double-Cross" imputation results are in bin matrix, take the bin matrix came from OutcrossSeq and plus a formated phenotype table, one can run GACD, MapQTL etc.

“Autopolyploid Plant” imputation result is a SNP matrix, one can run fastGWA, EMMAX etc. using SNP matrix from OutcrossSeq and formated phenotype table.

For more details, please check Demo run below.

#### <a name="Diploid"></a>Diploid-Outcrossing module
This module is designed for all diploid outcrossing species such as trees,fruits etc.
F1 population with about 200 individauls is good enough for a well study using this module.
It can also be applied to other types of populations when the founders are no more than two heterozygous lines.

*1:Filter out low-quality SNP sites*
```
cd /OutcrossSeq/diploid_outcrossing_species
perl diploid_cross_Genotype_filter.pl 5 parent1 parent2 100 0.6 demo_diploid 0.3

```

| Parameter |  Description|
| ------------- |------------- | 
| `5` |  chromosome |
| `parent1` | the id of parent1 in vcf file |
| `parent2` |the id of parent2 in vcf file  |
|`100`|quality threshold we set here as an example,quality below threshold will be filter|
| `0.6` | miss rate threshold we set here as an example,missing values greater than the threshold will be filtered  |
| `demo_diploid` | the input file of the vcf file |
| `0.3` | heterozygosity between 0.3 and (1-0.3)|

*2:Calculate kinship of each sample in every windows*
```
perl diploid_cross_step1.pl 100000 1000 F1-chr5.geno 5 0.35
```
| Parameter |  Description|
| ------------- |------------- | 
| `100000` | window size (<1cM) |
| `1000` | min marker number in each window, at least 1000|
| `F1-chr5.geno` |  the input file |
| `5` |  chromosome |
| `0.35` | start cuttree threshold |

*3:Clustering in each window*
```
perl diploid_cross_step2.pl 5 60 0.35
```
| Parameter |  Description|
| ------------- |------------- | 
| `5` | chromosome |
| `60` | sample number |
|`0.35`|start cuttree threshold |

*4:Haplotying across the chromosome*
In this step, a sample.list is demanded which records the order of individauls in genome VCF file.
```
#less sample.list
Sim1
Sim2
Sim3
Sim4
Sim5
Sim6
Sim7
Sim8
Sim9
Sim10
...
```
```
perl diploid_cross_step3.pl 5
```

| Parameter |  Description|
| ------------- |------------- | 
| `5` | chromosome |

*5:Generate the genotype file*
```
perl diploid_cross_step4.pl 5 5 demo
```
| Parameter |  Description|
| ------------- |------------- | 
| `5` | chromosome |
| `demo` | output file's name |

*Optional steps: One may want to check the haplotypes of parents using following command line*
```
perl diploid_cross_Haplotype.pl demo.geno demo_diploid parent1 parent2
```
| Parameter |  Description|
| ------------- |------------- | 
| `demo.geno` | the genotype of bin map |
| `demo_diploid` | the file of vcf file such as demo_diploid.vcf |
| `parent1` | the id of parent1 in vcf file |
| `parent2` |the id of parent2 in vcf file  |

#### <a name="Double"></a>Double-Cross module
This module is designed for typical Double-Cross population in crop breeding programs, such as rice and mazie etc.
Population size with 500 individauls is sufficient to carry out a nice study.
It can also be applied to other types of populations when the founders are no more than four inbred lines.

*1.filter out low-quality low-rate variants and select homozygous variants without missing data and filter the variants that two parents of hybridization have the same mutation and then transform vcf file to genotype file*
```
cd /OutcrossSeq/double_cross_populations
perl double_cross_Genotype_filter.pl  demo_double 100 16 17 18 19
```
| Parameter |  Description|
| ------------- |------------- | 
| `demo_double` | the input file of the vcf file |
| `100` | quality threshold we set here as an example,quality below threshold will be filter |
| `16 17 18 19` |column number of four parents;16—parent1,17—parent2,18—parent3,19—parent4(Caution: 16,17 should be inbred parent1 x parent2 ,the column 18,19 should be inbred parent3 x parent4. If mix up, the result will be wrong) |

*2:Determine the genotype of the first allele*
```
perl double_cross_step1.pl demo_double_parent.geno demo_double_offspring.geno  100000 5 1 6
```
| Parameter |  Description|
| ------------- |------------- | 
| `demo_double_parent.geno,demo_double_offspring.geno` | input file |
| `100000` | window size (<1cM) |
| `5` |  chromosome |
| `1` | the first sample order |
| `6` |the last sample order|

*3:Determine the genotype of the second allele*
```
perl double_cross_step2.pl demo_double_parent.geno demo_double_offspring.geno  100000 5 1 6
```
| Parameter |  Description|
| ------------- |------------- | 
| `demo_double_parent.geno,demo_double_offspring.geno` | input file |
| `100000` | window size (<1cM) |
| `5` |  chromosome |
| `1` | the first sample order |
| `6` |the last sample order|

*4:Generate the genotype file*
In this step, a sample.list is demanded which records the order of individauls in genome VCF file.
```
#less sample.list
Sim1
Sim2
Sim3
Sim4
Sim5
Sim6

```
```
perl double_cross_step3.pl 5 5 6 100000
```
| Parameter |  Description|
| ------------- |------------- | 
| `5` | star of chromosome |
| `5` | end of chromosome |
| `6` |  Size of the population |
| `100000` |window size (<1cM)|

#### <a name="Autopolyploid"></a>Autopolyploid Plant module
This module is designed dedicated for F1 population of autopolyploid or near-autopolyploid plants such as potato and sweet potato etc.
A F1 Population with about 1000 or more individauls is recommended,  for a reasonable practice.
It might be applied to other types of large populations but which is not typical in autopolyploid plants.

*1:Selecet simplex or double-simplex SNPs*
```
perl autopolyploid_Genotype_filter_step1.pl demo_autopolyploid-Parents 100 0.08 1
```
| Parameter |  Description|
| ------------- |------------- | 
| `demo_autopolyploid-Parents` |the input file of the parents vcf file |
| `100` | quality threshold we set here as an example,quality below threshold will be filter |
| `0.08,1` | depth range from 0.08 to 1 |

*2:Transform .vcf to .geno file*
```
perl autopolyploid_Genotype_filter_step2.pl 5 demo_autopolyplo5id Parent1 Parent2
```
| Parameter |  Description|
| ------------- |------------- |
| `5` |  chromosome |
| `demo_autopolyploid` |the input file of the offsprings vcf file |
| `Parent1` | the id of parent1 in vcf file |
| `Parent2` |the id of parent2 in vcf file  |

*3:Filter low-quality and error ratio sites*
**This step can be skipped**
```
perl autopolyploid_step1.pl 5 demo_autopolyploid
```
| Parameter |  Description|
| ------------- |------------- |
| `5` |  chromosome |
| `demo_autopolyploid` |the input file of the .geno file |

*4:Select relevant sites*
```
perl autopolyploid_step2.pl select-demo_autopolyploid 5 50 0.9 3 100
```
| Parameter |  Description|
| ------------- |------------- |
|select-demo_autopolyploid|the input file of the select*.geno file |
| `5` |  chromosome |
|`50`|min SNP marker number of each window|
|`0.9`|cuttree|
|`3`|select min cuttree marker num|
|`100`|select max cuttree marker num|

*Optimized steps: If the input file is too big, one may want to parallelize the computation via split the file into multiple runs*
```
perl split-select-geno.pl 5 110 50 select-demo_autopolyploid
```
| Parameter |  Description|
| ------------- |------------- |
| `5` |  chromosome |
| `110` |  File lines (Can be obtained by command:wc -l file_name) of select-demo_autopolyploid.geno |
|`50`|line number of plit file,be careful this parameter must be a multiple of the min snp number for each window|
|`select-demo_autopolyploid`|the split file name |

*Optimized steps: Then one has to combine all the results together when all runs finished*
```
perl combine_split.pl 5 select-demo_autopolyploid 2

```
| Parameter |  Description|
| ------------- |------------- |
| `5` |  chromosome |
|`select-demo_autopolyploid`|the split file name |
|`3`|Number of split files|

*5:Select mark readsNum *
In this step, a sample.list is demanded which records the order of individauls in genome VCF file.
```
#less sample.list
Sim1
Sim2
Sim3
Sim4
Sim5
Sim6
Sim7
Sim8
Sim9
Sim10
...
```
```
perl autopolyploid_step3.pl demo_autopolyploid-readsNum.txt R2-select-demo_autopolyploid
```
| Parameter |  Description|
| ------------- |------------- |
|`demo_autopolyploid-readsNum.txt`|input file1 |
|`R2-select-demo_autopolyploid`|input file2 named R2-select-demo_autopolyploid.txt|

*6:Imputation missing SNPs*
```
perl autopolyploid_step4.pl 5 R2-select-demo_autopolyploid-readsNum.txt R2-select-demo_autopolyploid.txt
```
| Parameter |  Description|
| ------------- |------------- |
|`demo_autopolyploid-readsNum.txt`|input file1 |
|`R2-select-demo_autopolyploid`|input file2|

*7:Merge genotype files*
```
perl autopolyploid_step5.pl 5 5
#because demo data only contains the chromosome 5
```
| Parameter |  Description|
| ------------- |------------- |
|`first 5`|start chromosome |
|`second 5`|end chromosome |
### <a name="help"></a>Getting help
Manpage [supplementary note](https://github.com/xhhuanglab/OutcrossSeq/blob/master/Supplemental%20note.md) provides more detailed descriptions and tips of OutcrossSeq.
### <a name="citing"></a>Citing OutcrossSeq
The paper describing OutcrossSeq can be found here:


