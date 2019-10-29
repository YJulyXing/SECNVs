# SECNVs 2.7 (SimulateCNVs 2.7)

**Maintainer: Yue "July" Xing**<br>
**Contact: yue.july.xing@gmail.com**<br>
**Version: 2.7**<br>
**Date: 09/14/2019**


## Description
A tool for simulating CNVs for WES data. It simulates rearranged genomes, short reads (fastq) and bam files automatically in one single command as desired by the user. There are several ways and distributions to choose from to generate desired CNVs.
Custom codes and algorithms were used to simulate rearranged genomes.
Short read simulation is based on the modified script of [Wessim1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3624799/).
The error model is constructed using a modified script of [GimSIM](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-74).
BAM file generation uses [BWA](http://bio-bwa.sourceforge.net/), [samtools](http://samtools.sourceforge.net/), [picard](https://broadinstitute.github.io/picard/) and [GATK](https://software.broadinstitute.org/gatk/).

## Version 2.7 Update (09/14/2019):
* Updated the default error profile. Now it is Illumina HiSeq 2500 WXS paired end sequencing. The statistical reports for this error profile is also included.
* Included the modified scripts from [GimSIM](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-74) and instructions on how to make user-specific error profiles and generate statistical reports for the error profiles.

## Version 2.6 Update (08/31/2019):
* Users now can choose to only randomly replace gap regions (minimum length for defining the gap can be determined by the user), replace all "N"s or doesn't do any replacement.
* Users now can choose to only avoid gap regions (minimum length for defining the gap can be determined by the user), avoid all "N"s or doesn't avoid any "N"s.
* The control genome file randomly replaced "N"s/gap regions is now in the output directory with name "control.fa", instead of overwritting the original reference genome file.
* Modified the functions of SNP simulation to increase speed.
* Bug fixes.

## Version 2.5 Update (08/21/2019):
* Added an option to add slack regions up and down stream of target regions to simulate SNPs.
* The random simulation for SNPs is now weighted.
* Modified the functions of SNP and indel simulation to increase speed.
* Bug fixes.

## Version 2.4 Update (07/30/2019):
* Added indel simulation.

## Version 2.3 Update (07/16/2019):
* Changed SNP simulation from simulating on the whole genome to only in the target regions, which highly increase the speed.

## Version 2.2 Update (06/28/2019):
* Incorporated "replaceNs.py" to the main program as an option "-rN". No longer need to use "replaceNs.py" separately.
* Bug fixes.

## Version 2.1 Update (06/07/2019):
* Added Beta distribution for CNV length. User can define alpha and beta.
* User can now define mu and sigma for Gaussian distribution.
* Added feature: SNP simulation.
* Bug fixes.

## Version 2.0 Update (05/23/2019):
* Use a modified script of Wessim1 instead of ART_illumina to simulate short reads. This enables GC filtration in the step of short read generation. The modification of Wessim's script fixed the following issues in Wessim1:<br>
&#160;1) Only generate short reads at the start and end of each region, no reads generated in the middle of the regions.<br>
&#160;2) Can't read in the first line of the bed file with no header.<br>
&#160;3) The length of each region was 1 bp less than it should be.<br>
* Also, after modification Wessim1 doesn't require the genome to be indexed by "faidx", doesn't need python package "pysam", and the intermediate files can now be deleted automatically.

## For the old version SimulateCNVs 1.0, please click [here](https://yjulyxing.github.io/).

## Installation
No installation required. Simply type the following and it will be ready to use:
``` bash
git clone https://github.com/YJulyXing/SECNVs.git
```
and unpack it using "tar".
Or manually download the source codes [here](https://github.com/YJulyXing/SECNVs).

## Requirements
* General use: [Python 2.7](https://www.python.org/download/releases/2.7/). Required python packages: argparse, random, os, subprocess, math, sys, time, copy, numpy
* If to use Beta distribution for CNV length, and not using default parameters for Beta distribution, [R](https://www.r-project.org/) can be used to estimate the parameters for Beta distribution using a user-specified set of CNV lengths. (need packages MASS or fitdistrplus)
* To generate short reads (fastq) outputs (see requirements for [Wessim](http://sak042.github.io/Wessim/): <br>
&#160;1. [Python 2.7](https://www.python.org/download/releases/2.7/). Required python packages: bisect, gzip, cPickle, numpy, multiprocessing <br>
&#160;2. [GemSIM](https://sourceforge.net/projects/gemsim/) error models. The default error model is a Illumina HiSeq 2500 WXS paired end sequencing model made using a modified script of GemSIM. The statistical reports for this error profile is included. Users can also make their own error models using real data by the modified GimSIM script. Scripts for doing so is included. See below for instructions. [Python 2.7](https://www.python.org/download/releases/2.7/) is required for making user-specific error models. Required python packages: sys, getopt, cPickle, gzip, logging, numpy
* To generate BAM outputs: <br>
&#160;1. [Samtools](http://samtools.sourceforge.net/) <br>
&#160;2. [BWA](http://bio-bwa.sourceforge.net/) <br>
&#160;3. [JAVA](https://www.java.com/en/) <br>
&#160;4. [picard 2.15.0](https://broadinstitute.github.io/picard/) <br>
&#160;5. [GATK](https://software.broadinstitute.org/gatk/)

## Usage
``` bash
usage: SECNVs.py [-h] -G GENOME_FILE -T TARGET_REGION_FILE
                 [-rN {none,gap,all}] [-eN {none,gap,all}]
                 [-n_gap NUM_N_FOR_GAPS] [-e_cnv TARGET_CNV_LIST]
                 [-e_chr TARGET_CNV_CHR] [-e_tol TARGET_CNV_TOL]
                 [-e_cl TARGET_CNV_LEN_FILE] [-o_cnv OUT_CNV_LIST]
                 [-o_chr OUT_CNV_CHR] [-o_tol OUT_CNV_TOL]
                 [-o_cl OUT_CNV_LEN_FILE] [-ol OVERLAP_BPS]
                 [-min_len CNV_MIN_LENGTH] [-max_len CNV_MAX_LENGTH]
                 [-min_cn MIN_COPY_NUMBER] [-max_cn MAX_COPY_NUMBER]
                 [-p PROPORTION_INS] [-f MIN_FLANKING_LEN]
                 [-ms {random,uniform,gauss}]
                 [-ml {random,uniform,gauss,beta,user}] [-as AS1] [-bs BS]
                 [-al AL] [-bl BL] [-s_r S_RATE] [-s_s S_SLACK] [-i_r I_RATE]
                 [-i_mlen I_MAX_LEN] [-nr NREADS] [-fs FRAG_SIZE] [-s STDEV]
                 [-l READ_LENGTH] [-tf TARGET_REGION_FLANK] [-pr]
                 [-q QUALITY_SCORE_OFFSET] [-clr CONNECT_LEN_BETWEEN_REGIONS]
                 [-m MODEL] [-o OUTPUT_DIR] [-rn REARRANGED_OUTPUT_NAME]
                 [-n NUM_SAMPLES] [-sc] [-ssr] [-sb] [-picard PATH_TO_PICARD]
                 [-GATK PATH_TO_GATK]
```

## Arguments
#### Optional arguments:

|   Parameter                    |     Default value     |    Explanation                             | Restrictions |
| :----------------------------: | :-------------------: | :----------------------------------------- | :----------- |
| -h, --help | - | show this help message and exit | - |

#### Mandatory arguments:

|   Parameter                    |     Default value     |    Explanation                             | Restrictions |
| :----------------------------: | :-------------------: | :----------------------------------------- | :----------- |
| -G GENOME_FILE | - | Reference genome FASTA file  | - |
| -T TARGET_REGION_FILE | - | Target region file | - |

#### Arguments for simulating rearranged genomes:

|   Parameter                    |     Default value     |    Explanation                             | Restrictions |
| :----------------------------: | :-------------------: | :----------------------------------------- | :----------- |
| -rN {none,gap,all} | none | Replace sequences of "N"s by ATGC randomly? | Output is the "control.fa" file in output directory for any of the choices. This file is used as control for all steps. None: no change; gap: replace gap regions only; all: replace all "N"s. |
| -eN {none,gap,all} | none | Exclude sequences of "N"s for CNV simulation? | None: does not exclude anything; gap: exclude gap regions only; all: exclude all "N"s. Takes place after -rN option. |
| -n_gap NUM_N_FOR_GAPS | 50 |  Number of consecutive "N"s to be considered a gap region | - |
| -e_cnv TARGET_CNV_LIST | - | A user-defined list of CNVs overlapping with target regions | One and only one of -e_cnv, -e_chr, -e_tol and -e_cl can be used with WES simulation to generate CNVs overlapping with target regions.<br> If -e_cnv is provided, -em, -f, -ms, -ml, -ol, -min_cn, -max_cn, -min_len and -max_len will be ignored for CNVs overlapping with target regions. |
| -e_chr TARGET_CNV_CHR | - | Number of CNVs overlapping with target regions to be generated on each chromosome | Same as above. |
| -e_tol TARGET_CNV_TOL | - | Total number of CNVs overlapping with target regions to be generated across the genome (an estimate) |  Same as above. |
| -e_cl TARGET_CNV_LEN_FILE | - | User supplied file of CNV length for CNVs overlapping with target regions | Must be used with -ml user. Can’t be used with -o_cnv, -o_chr and -o_tol. Otherwise same as above. |
| -o_cnv OUT_CNV_LIST | - | A user-defined list of CNVs outside of target regions | One and only one of -o_cnv, -o_chr, -o_tol and -o_cl can be used with WES simulation to generate CNVs outside of target regions.<br> If -o_cnv is provided, -em, -f, -ms, -ml, -ol, -min_cn, -max_cn, -min_len and -max_len will be ignored for CNVs outside of target regions. |
| -o_chr OUT_CNV_CHR | - | Number of CNVs outside of target regions to be generated on each chromosome |  Same as above. |
| -o_tol OUT_CNV_TOL | - | Total number of CNVs outside of target regions to be generated across the genome (an estimate) |  Same as above. |
| -o_cl OUT_CNV_LEN_FILE | - | User supplied file of CNV length for CNVs outside of target regions | Must be used with -ml user. Can’t be used with -e_cnv, -e_chr and -e_tol. Otherwise same as above. |
| -ol OVERLAP_BPS | 100 | For each CNV overlapping with target regions, number of minimum overlapping bps | - |
| -min_len CNV_MIN_LENGTH | 1000 |  Minimum CNV length in bps | - |
| -max_len CNV_MAX_LENGTH | 100000 | Maximum CNV length in bps | - |
| -min_cn MIN_COPY_NUMBER | 2 | Minimum copy number for insertions | - |
| -max_cn MAX_COPY_NUMBER | 10 | Maximum copy number for insertions | - |
| -p PROPORTION_INS | 0.5 | Proportion of insertions | Must be between 0 and 1 |
| -f MIN_FLANKING_LEN | 50 |  Minimum length between each CNV | - |
| -ms {random,uniform,gauss} | random | Distribution of CNVs | - |
| -ml {random,uniform,gauss,beta,user} | random | Distribution of CNV length | -ml user must be used with -e_cl and/or -o_cl. If -ml user is used, -min_len and -max_len will be ignored. |
| -as AS1 | 0 | Mu for gauss CNV distribution | For other choices of -ms and -ml, this parameter will be ignored. |
| -bs BS | 1 | Sigma for Gaussian CNV distribution | For other choices of -ms and -ml, this parameter will be ignored. |
| -al AL | 0 for Gaussian distribution, and 0.5 for Beta distribution | Mu (Gaussian distribution) or alpha (Beta distribution) for CNV length distribution | If user has a set of CNV lengths, he/she can use "fitdistr" or "fitdist"in [R](https://www.r-project.org/) to estimate the value of the parameters. See below for instructions.<br> For other choices of -ms and -ml, this parameter will be ignored. |
| -bl BL | 1 for Gaussian distribution, and 2.3 for Beta distribution | Sigma (Gaussian distribution) or beta (Beta distribution) for CNV length distribution | If user has a set of CNV lengths, he/she can use "fitdistr" or "fitdist"in [R](https://www.r-project.org/) to estimate the value of the parameters. See below for instructions.<br> For other choices of -ms and -ml, this parameter will be ignored. |
| -s_r S_RATE | 0 | Rate of SNPs in target regions | Must be between 0 and 1 |
| -s_s S_SLACK | 0 | Slack region up and down stream of target regions to simulate SNPs | Must >= 0 |
| -i_r I_RATE | 0 | Rate of indels in target regions | Must be between 0 and 1 |
| -i_mlen I_MAX_LEN | 50 | The Maximum length of indels in target regions | If a deletion is equal or larger than the length of the target region it is in, the length of the deletion will be changed to (length of the target region it is in) - 1.|

#### Arguments for simulating short reads (fastq):

|   Parameter                    |     Default value     |    Explanation                             | Restrictions |
| :----------------------------: | :-------------------: | :----------------------------------------- | :----------- |
| -nr NREADS | 10000 | Number of reads / read pairs on target regions to be generated for each genome | - |
| -fs FRAG_SIZE | 200 | Mean fragment size to be generated | - |
| -s STDEV | 20 | Standard deviation of fragment sizes | - |
| -l READ_LENGTH | 100 | Read length of each short read | - |
| -tf TARGET_REGION_FLANK | 0 | Length of flanking region up and down stream of target regions to be sequenced | -tf takes place after -clr (target regions are first connected per request of the user, and then flanking regions up and down stream of target regions are included for sequencing). It is recommended to use -tf and -clr if target regions are very small and fragment size and read length are relatively large. |
| -pr | - | Select if paired-end sequencing | For paired-end sequencing, must use the option -pr. If using single-end sequencing, mean fragment size (-fs) and standard deviation of fragment size (-s) will be ignored. |
| -q QUALITY_SCORE_OFFSET | 33 | Quality score offset for short reads simulation | - |
| -clr CONNECT_LEN_BETWEEN_REGIONS | - | Maximum length bwtween target regions to connect the target regions | -tf takes place after -clr (target regions are first connected per request of the user, and then flanking regions up and down stream of target regions are included for sequencing). It is recommended to use -tf and -clr if target regions are very small and fragment size and read length are relatively large. |
| -m MODEL | Illumina_HiSeq_2500_p | GemSIM error model file (.gzip, need absolute path) | - |

#### Arguments for general settings:

|   Parameter                    |     Default value     |    Explanation                             | Restrictions |
| :----------------------------: | :-------------------: | :----------------------------------------- | :----------- |
| -o OUTPUT_DIR | simulation_output | Output directory | Will be generated if not exist. |
| -rn REARRANGED_OUTPUT_NAME | test | Prefix of the rearranged outputs | Do not include directory name. |
| -n NUM_SAMPLES | 1 | Number of test samples to be generated | Must be >= 1. |
| -sc | - | Simulation for control genome | - |
| -ssr | - | Simulate short reads (fastq) files | If the final output is BAM file(s), must first simulate short reads (-ssr) and then BAM (-sb). |
| -sb | - | Simulate BAM files | Same as above. |
| -picard PATH_TO_PICARD | - | Absolute path to picard | If the final output is BAM file(s) (using -ssr and -sb), must provide absolute paths to picard and GATK. |
| -GATK PATH_TO_GATK | - | Absolute path to GATK | Same as above. |

## Inputs
1. Sequence of a reference genome in fasta format. Anything after "_" or " " at header line will be ignored.
2. A tab delimited file of target regions in the order of chromosome, start and end. Header should not be included. Chromosome name should strictly match chromosome name in (1). Target regions should be sorted in ascending order.
3. For -e_cnv and -o_cnv, a tab delimited file generated by SimulateCNVs, in the order of chromosome, start, end, length and copy number. Header should be included. Chromosome name should strictly match chromosome name in (1).
4. For -e_cl and -o_cl, a tab delimited file of desired CNV lengths in the order of chromosome, CNV length and number of CNV of that length. Header should not be included. Chromosome name should strictly match chromosome name in (1).

## Outputs
1. Rearranged genome(s) (fasta)<br>
Control genome (fasta)<br>
2. List(s) of CNVs overlapping with target regions (bed)<br>
List(s) of CNVs outside of target regions (bed)<br>
3. Target regions for generating short reads for test(s) and control (\*.target_regions_for_gen_short_reads.bed). These rearranged target regions are for read generation ONLY. If don't want to use modified Wessim1, these can be used to extract target sequences and generate short reads by other read simulation tools like Wessim2 or ART_illumina (custom codes will be required for this, see [SimulateCNVs v1.0](https://github.com/YJulyXing/SimulateCNVs) which used a custom script to utilize ART_illumina). <strong>For CNV calling, the original target region file in input (2) should be used.</strong><br>
4. Short reads for rearranged genome(s) (fastq, if -ssr is chosen)<br>
Short reads for control genome (fastq, if -sc and -ssr is chosen)
5. Indexes for the control genome (.dict, .fai, .sa, etc., if -ssr and -sb is chosen and no indexes exist in the output directory)
6. BAM file(s) and index(es) for rearranged genome(s) (bam and bai, if -ssr and -sb is chosen)<br>
BAM file(s) and index(es) for control genome (bam and bai, if -sc, -ssr and -sb is chosen)

## Find parameters for Beta distribution from a set of known CNV lengths using R
* Use "fitdistr" in MASS library:
```
# x is an vector of CNV lengths, used to estimate parameters for Beta distribution.
x = c(10,20,30,40,35,25,11,37)
xs = (x-min(x))/(max(x)-min(x))
xs = xs[-which.max(x)]
xs = xs[-which.min(x)]

library(MASS)
fit=fitdistr(xs, "beta", list(shape1 = 0.5, shape2= 0.5))
print(fit)
```

* Or use "fitdist" in fitdistrplus library:
```
# x is an vector of CNV lengths, used to estimate parameters for Beta distribution.
x = c(10,20,30,40,35,25,11,37)
xs = (x-min(x))/(max(x)-min(x))
xs = xs[-which.max(x)]
xs = xs[-which.min(x)]

library(fitdistrplus)
fit <- fitdist(xs, "beta")
print(fit)
# To check fitness:
plot(fit, las = 1)
```

* Shape1 is alpha an shape2 is beta, which can be used in Beta distribution in SECNVs.

## Make user-specific error profiles
### To generate error model:
``` bash
python ~/SECNVs/GemSIM/GemErrxy.py \
      -r <maximium read length> \
      -f <reference genome (fasta)> \
      -s <SAM file for making error model with> \
      -n <prefix of error model> \
      -p
```
* Use "-p" only when data is paired!
* For extra options, see GimSIM manual

### To generate statistical reports:
``` bash
python ~/SECNVs/GemSIM/GemStats.py \
      -m <error model file name> \
      -n <prefix of the reports> \
      -p
```
* Use "-p" only when data is paired!

## Change SNP mutation rate
1. Open "snp_rate.py".
2. On line 14-18, find:
``` bash
"A": list(choices(["C","T","G"],1,p=[0.14, 0.04, 0.82])),
"T": list(choices(["C","A","G"],1,p=[0.84, 0.03, 0.13])),
"G": list(choices(["C","A","T"],1,p=[0.19, 0.7, 0.11])),
"C": list(choices(["G","A","T"],1,p=[0.17, 0.12, 0.71])),
```
The mutation rate from "A" to "C", "T" or "G" are the first, second and third numbers in p=[0.14, 0.04, 0.82] on "A" line, respectively. Change these numbers according to your need.
The mutation rate from "T" to "C", "A" or "G" are the first, second and third numbers in p=[0.84, 0.03, 0.13] on "T" line, respectively. Change these numbers according to your need.
The mutation rate from "G" to "C", "A" or "T" are the first, second and third numbers in p=[0.19, 0.7, 0.11] on "G" line, respectively. Change these numbers according to your need.
The mutation rate from "C" to "G", "A" or "T" are the first, second and third numbers in p=[0.17, 0.12, 0.71] on "C" line, respectively. Change these numbers according to your need.
3. Save "snp_rate.py".

## Examples
1. Simulate 10 CNVs overlapping with target regions, and 1 CNV outside of target regions randomly on each chromosome using default lengths, copy numbers, minimum distance between each of the 2 CNVs and proportion of insertions. For each CNV overlapping with target regions, the overlapping length is no less than 90 bps. CNV break points follow a Gaussian(1, 2) distribution, and CNV lengths follow a Beta(2, 5) distribution. CNVs are not generated in gaps. A total of 5 test and control samples are built. Short reads (fastq) files are generated using default settings, paired-end sequencing.
``` bash
SECNVs/SECNVs.py -G <input_fasta> -T <target_region> -o <output_dir> \
                  -e_chr 10 -o_chr 1 -ol 90 -ms gauss -as 1 -bs 2 -ml beta -al 2 -bl 5 -eN gap -n 5 -sc -pr -ssr
```

2. Simulate CNVs overlapping with target regions from a provided CNV list. Twenty CNVs are to be simulated outside of target regions randomly on the whole genome with default settings. CNVs are not to be generated on any stretches of "N"s. A pair of test and control genome are built.
``` bash
SECNVs/SECNVs.py -G <input_fasta> -T <target_region> -o <output_dir> \
                  -e_cnv <list_of_CNV_overlapping_with_target_regions> -o_tol 20 -eN all -sc
```

3.	Simulate 20 CNVs overlapping with target regions on the whole genome, and have at least 100 bps between any 2 CNVs. CNVs are not generated outside of target regions. Gaps (50 or more consecutive "N"s) are replaced by random nucleotides. SNP rate is 0.001 and indel rate is 0.00001, and the maximum indel length is 100 bps. Paired-end sequencing reads with quality offset 35 are then produced. For a pair of test and control genomes BAM files are generated.
``` bash
SECNVs/SECNVs.py -G <input_fasta> -T <target_region> -o <output_dir> \
                  -e_tol 20 -f 100 -rN gap -sc -pr -q 35 -ssr -sb \
                  -s_r 0.001 -i_r 0.00001 -i_mlen 100 \ 
                  -picard <absolute_path_to_picard> -GATK <absolute_path_to_GATK>
```

4.	Simulate CNVs overlapping with target regions and outside of target regions from provided files of CNV lengths. Combined single regions are formed from two or more regions originally separated by less than 100 bps. CNVs are not generated on gaps (60 or more consecutive "N"s). A total of 10 test and control samples are built. The paired-end sequencing must include sequences 50 bp upstream and downstream of the target regions. The final output consists of short reads (fastq) files with 100,000 reads.
``` bash
SECNVs/SECNVs.py -G <input_fasta> -T <target_region> -o <output_dir> \
                  -ml user -e_cl <length_file_1> -o_cl <length_file_2> \
                  -clr 100 -eN gap -n_gap 60 -n 10 -sc -pr -tf 50 -nr 100000 -ssr 
```
