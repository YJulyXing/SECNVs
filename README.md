# SECNVs 2.2 (SimulateCNVs 2.2)

**Maintainer: Yue "July" Xing**<br>
**Author: Yue "July" Xing**<br>
**Contact: yue.july.xing@gmail.com**<br>
**Version: 2.2**<br>
**Date: 06/028/2019**


## Description
A tool for simulating CNVs for WES data. It simulates rearranged genomes, short reads (fastq) and bam files automatically in one single command as desired by the user. There are several ways and distributions to choose from to generate desired CNVs.
Short read simulation is based on the modified script of Wessim.
Bam file generation uses BWA, samtools, picard and GATK.

## Version 2.2 Update (06/28/2019):
* Incorporated "replaceNs.py" to the main program as an option "-rN". No longer need to use "replaceNs.py" separately.
* Bug fixes.

## Version 2.1 Update (06/07/2019):
* Added beta distribution for CNV length. User can define alpha and beta.
* User can now define mu and sigma for gauss distribution.
* Added feature: SNP simulation.
* Bug fixes.


## Version 2.0 Update (05/23/2019):
* Use a modified script of Wessim1 instead of ART_illumina to simulate short reads. This enables GC filtration in the step of short read generation. The modification of Wessim's script fixed the following issues in Wessim1:<br>
&#160;1) Only generate short reads at the start and end of each region, no reads generated in the middle of the regions.<br>
&#160;2) Can't read in the first line of the bed file with no header.<br>
&#160;3) The length of each region was 1 bp less than it should be.<br>
* Also, after modification Wessim1 doesn't require the genome to be indexed by "faidx", doesn't need python package "pysam", and the intermediate files can now be deleted automatically.

## Installation
No installation required. Simply type the following and it will be ready to use:
``` bash
git clone https://github.com/YJulyXing/SECNVs.git
```
and unpack it using "tar".
Or manually download the source codes [here](https://github.com/YJulyXing/SECNVs).

## Requirements
* General use: [Python 2.7](https://www.python.org/download/releases/2.7/). Required python packages: argparse, random, os, subprocess, math, sys, time
* To generate short reads (fastq) outputs (see requirements for [Wessim](http://sak042.github.io/Wessim/): <br>
&#160;1. [Python 2.7](https://www.python.org/download/releases/2.7/). Required python packages: bisect, gzip, cPickle, numpy, multiprocessing <br>
&#160;2. [GemSim](https://sourceforge.net/projects/gemsim/) error models. This tool by default uses the GemSim model "Illumina  GA  IIx  with  TrueSeq  SBS  Kit  v5‐GA,  paired reads", and has "Illumina  GA  IIx  with  TrueSeq  SBS  Kit  v5‐GA,  single reads" bulit in as well. Users can make their own error models using real data by GemErr.py  in GemSim.
* To generate bam outputs: <br>
&#160;1. [Samtools](http://samtools.sourceforge.net/) <br>
&#160;2. [BWA](http://bio-bwa.sourceforge.net/) <br>
&#160;3. [JAVA](https://www.java.com/en/) <br>
&#160;4. [picard 2.15.0](https://broadinstitute.github.io/picard/) <br>
&#160;5. [GATK](https://software.broadinstitute.org/gatk/)

## Usage
``` bash
usage: usage: SECNVs.py [-h] -G GENOME_FILE -T TARGET_REGION_FILE [-rN] [-em]
                        [-e_cnv EXON_CNV_LIST] [-e_chr EXON_CNV_CHR]
                        [-e_tol EXON_CNV_TOL] [-e_cl EXON_CNV_LEN_FILE]
                        [-o_cnv OUT_CNV_LIST] [-o_chr OUT_CNV_CHR]
                        [-o_tol OUT_CNV_TOL] [-o_cl OUT_CNV_LEN_FILE]
                        [-ol OVERLAP_BPS] [-min_len CNV_MIN_LENGTH]
                        [-max_len CNV_MAX_LENGTH] [-min_cn MIN_COPY_NUMBER]
                        [-max_cn MAX_COPY_NUMBER] [-p PROPORTION_INS]
                        [-f MIN_FLANKING_LEN] [-ms {random,uniform,gauss}]
                        [-ml {random,uniform,gauss,beta,user}] [-as AS1] [-bs BS]
                        [-al AL] [-bl BL] [-r RATE] [-nr NREADS] [-fs FRAG_SIZE]
                        [-s STDEV] [-l READ_LENGTH] [-tf TARGET_REGION_FLANK] [-pr]
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
| -rN | - | Replace gap regions (Ns) by ATGC randomly | This option will create a copy of the the original genome file with name genome_file_copy.fa in the same directory, and replace the original one. |
| -em | - | Exclude gap sequences for CNV simulation | - |
| -e_cnv TARGET_CNV_LIST | - | A user-defined list of CNVs overlapping with target regions | One and only one of -e_cnv, -e_chr, -e_tol and -e_cl can be used with WES simulation to generate CNVs overlapping with target regions.<br> If -e_cnv is provided, -em, -f, -ms, -ml, -ol, -min_cn, -max_cn, -min_len and -max_len will be ignored for CNVs overlapping with target regions. |
| -e_chr TARGET_CNV_CHR | - | Number of CNVs overlapping with target regions to be generated on each chromosome | Same as above. |
| -e_tol TARGET_CNV_TOL | - | Total number of CNVs overlapping with target regions to be generated across the genome (an estimate) |  Same as above. |
| -e_cl TARGET_CNV_LEN_FILE | - | User supplied file of CNV length for CNVs overlapping with target regions | Must be used with -ml user. Can’t be used with -o_cnv, -o_chr and -o_tol. Otherwise same as above. |
| -o_cnv OUT_CNV_LIST | - | A user-defined list of CNVs outside of target regions | One and only one of -o_cnv, -o_chr, -o_tol and -o_cl can be used with WES simulation to generate CNVs outside of target regions.<br> . If -o_cnv is provided, -em, -f, -ms, -ml, -ol, -min_cn, -max_cn, -min_len and -max_len will be ignored for CNVs outside of target regions. |
| -o_chr OUT_CNV_CHR | - | Number of CNVs outside of target regions to be generated on each chromosome |  Same as above. |
| -o_tol OUT_CNV_TOL | - | Total number of CNVs outside of target regions to be generated across the genome (an estimate) |  Same as above. |
| -o_cl OUT_CNV_LEN_FILE | - | User supplied file of CNV length for CNVs outside of target regions | Must be used with -ml user. Can’t be used with -e_cnv, -e_chr and -e_tol. Otherwise same as above. |
| -ol OVERLAP_BPS | 100 | For each CNV overlapping with target regions, number of minimum overlapping bps | - |
| -min_len CNV_MIN_LENGTH | 1000 |  Minimum CNV length in bps | - |
| -max_len CNV_MAX_LENGTH | 100000 | Maximum CNV length in bps | - |
| -min_cn MIN_COPY_NUMBER | 2 | Minimum copy number for insertions | - |
| -max_cn MAX_COPY_NUMBER | 10 | Maximum copy number for insertions | - |
| -p PROPORTION_INS | 0.5 | Proportion of insertions | - |
| -f MIN_FLANKING_LEN | 50 |  Minimum length between each CNV | - |
| -ms {random,uniform,gauss} | random | Distribution of CNVs | - |
| -ml {random,uniform,gauss,beta,user} | random | Distribution of CNV length | -ml user must be used with -e_cl and/or -o_cl. If -ml user is used, -min_len and -max_len will be ignored. |
| -as AS1 | 0 | Mu for gauss CNV distribution | For other choices of -ms and -ml, this parameter will be ignored. |
| -bs BS | 1 | Sigma for gauss CNV distribution | For other choices of -ms and -ml, this parameter will be ignored. |
| -al AL | 0 for gauss distribution, and 2 for beta distribution | Mu (gauss distribution) or alpha (beta distribution) for CNV length distribution | For other choices of -ms and -ml, this parameter will be ignored. |
| -bl BL | 1 for gauss distribution, and 2 for beta distribution | Sigma (gauss distribution) or beta (beta distribution) for CNV length distribution | For other choices of -ms and -ml, this parameter will be ignored. |
| -r RATE | 0 | Rate of SNPs | - |

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
| -m MODEL | ill100v5_p | GemSim error model file (.gzip, need absolute path) | - |

#### Arguments for general settings:

|   Parameter                    |     Default value     |    Explanation                             | Restrictions |
| :----------------------------: | :-------------------: | :----------------------------------------- | :----------- |
| -o OUTPUT_DIR | simulation_output | Output directory | Will be generated if not exist. |
| -rn REARRANGED_OUTPUT_NAME | test | Prefix of the rearranged outputs | Do not include directory name. |
| -n NUM_SAMPLES | 1 | Number of test samples to be generated | Must be >= 1. |
| -sc | - | Simulation for control genome | - |
| -ssr | - | Simulate short reads (fastq) files | If the final output is bam file(s), must first simulate short reads (-ssr) and then bam (-sb). |
| -sb | - | Simulate bam files | Same as above. |
| -picard PATH_TO_PICARD | - | Absolute path to picard | If the final output is bam file(s) (using -ssr and -sb), must provide absolute paths to picard and GATK. |
| -GATK PATH_TO_GATK | - | Absolute path to GATK | Same as above. |

## Inputs
1. Sequence of a reference genome in fasta format. Anything after "_" or " " at header line will be ignored.
2. A tab delimited file of target regions in the order of chromosome, start and end. Header should not be included. Chromosome name should strictly match chromosome name in (1). Target regions should be sorted in ascending order.
3. For -e_cnv and -o_cnv, a tab delimited file generated by SimulateCNVs, in the order of chromosome, start, end, length and copy number. Header should be included. Chromosome name should strictly match chromosome name in (1).
4. For -e_cl and -o_cl, a tab delimited file of desired CNV lengths in the order of chromosome, CNV length and number of CNV of that length. Header should not be included. Chromosome name should strictly match chromosome name in (1).

## Outputs
1. Rearranged genome(s) (fasta)<br>
Control genome (fasta, if -sc is chosen)<br>
2. List(s) of CNVs overlapping with target regions (bed)<br>
List(s) of CNVs outside of target regions (bed)<br>
3. Target regions for generating short reads for test(s) and control (\*.target_regions_for_gen_short_reads.bed). These rearranged target regions are for read generation ONLY. If don't want to use Wessim, these can be used to extract target sequences and generate short reads by other read simulation tools like ART_illumina (custom codes will be required for this, see [SimulateCNVs v1.0](https://github.com/YJulyXing/SimulateCNVs) which used a custom script to utilize ART_illumina). <strong>For CNV calling, the original target region file in input (2) should be used.</strong><br>
4. Short reads for rearranged genome(s) (fastq, if -ssr is chosen)<br>
Short reads for control genome (fastq, if -sc and -ssr is chosen)
5. Indexes for the control genome (.dict, .fai, .sa, etc., if -ssr and -sb is chosen and no indexes exist in the output directory)
6. Bam file(s) and index(es) for rearranged genome(s) (bam and bai, if -ssr and -sb is chosen)<br>
Bam file(s) and index(es) for control genome (bam and bai, if -sc, -ssr and -sb is chosen)

## Examples
1. Simulate 10 CNVs overlapping with target regions, and 1 CNV outside of target regions randomly on each chromosome using default lengths, copy numbers, minimum distance between each of the 2 CNVs and proportion of insertions. For each CNV overlapping with target regions, the overlapping length is not less than 90 bps. CNV start points following Gauss(1, 2) distribution, and lengths follow Beta(2, 5) distribution. Don’t generate CNVs on gap regions. Make 5 test samples and control. Generate short reads (fastq) files by default settings, using paired-end sequencing.
``` bash
SECNVs/SECNVs.py -G <input_fasta> -T <target_region> -o <output_dir> \
                  -e_chr 10 -o_chr 1 -ol 90 -ms gauss -as 1 -bs 2 -ml beta -al 2 -bl 5 -em -n 5 -sc -pr -ssr
```

2. Simulate CNVs overlapping with target regions from a provided CNV list. Simulate approximately 20 CNVs outside of target regions randomly on the whole genome with default settings. Don’t generate CNVs on gap regions. Make a pair of test and control genome.
``` bash
SECNVs/SECNVs.py -G <input_fasta> -T <target_region> -o <output_dir> \
                  -e_cnv <list_of_CNV_overlapping_with_exons> -o_tol 20 -em -sc 
```

3.	Simulate approximately 20 CNVs overlapping with target regions on the whole genome, and at least 100 bps between any 2 CNVs. Don’t generate CNVs outside of target regions. Replace gap regions randomly by nucleotides. Paired-end sequencing, with quality offset 35. Make a pair of test and control. The final outputs are BAM files.
``` bash
SECNVs/SECNVs.py -G <input_fasta> -T <target_region> -o <output_dir> \
                  -e_tol 20 -f 100 -rN -sc -pr -q 35 -ssr -sb \
                  -picard <absolute_path_to_picard> -GATK <absolute_path_to_GATK>
```

4.	Simulate CNVs overlapping with target regions and outside of target regions from provided files of CNV lengths. If the length between 2 target regions are smaller than 100 bps, connect them as 1 target region. Don’t generate CNVs on gap regions. Make 10 test samples and control. Use paired-end sequencing; sequence 50 bp up and down stream of the target regions (after connecting the target regions) as well. The final output is short reads (fastq) files with 100000 reads.
``` bash
SECNVs/SECNVs.py -G <input_fasta> -T <target_region> -o <output_dir> \
                  -ml user -e_cl <length_file_1> -o_cl <length_file_2> \
                  -clr 100 -em -n 10 -sc -pr -tf 50 -nr 100000 -ssr 
```


***


## ReplaceNs.py (Not needed for SECNVs version 2.2)

**Maintainer: Yue "July" Xing**<br>
**Author: Yue "July" Xing**<br>
**Version: 1.0**<br>
**Date: 06/27/2018**

### Description
A small program to fix genomes which have too many ‘N’s to generate desired CNVs. It replaces all ‘N’s in the genome sequence to ‘A’s, ‘T’s, ‘G’s, or ‘C’s randomly.

### Installation
It is included in the package of SECNVs.

### Requirements
[Python 2.7](https://www.python.org/download/releases/2.7/)

### Usage
``` bash
ReplaceNs.py [-h] -i INPUT_FASTA_FILE -o OUTPUT_FASTA_FILE
```
### Arguments
-h: help<br>
-i: input genome sequence in fasta format<br>
-o output genome sequence in fasta format
