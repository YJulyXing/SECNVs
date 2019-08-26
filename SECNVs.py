#!/usr/bin/python

import argparse
import random
import os
import subprocess
import math
import sys
import time
import copy
from numpy.random import choice as choices

from WES_simulator import *

def main():
	parser = argparse.ArgumentParser(description='Simulator for WES or WGS data', \
		formatter_class=argparse.RawTextHelpFormatter)
	group1 = parser.add_argument_group('Mandatory inputs')
	group1.add_argument('-G', type=str, dest='genome_file', required=True, \
		help='Reference genome FASTA file')
	group1.add_argument('-T', type=str, dest='target_region_file', required=True, \
		help='Target region file')

	group2 = parser.add_argument_group('Arguments for simulating rearranged genomes')
	group2.add_argument('-rN', dest='replaceNs', action='store_true', \
		help='Replace gap regions (Ns) by ATGC randomly. Caution: This option will create a copy of the the original genome file with name genome_file_copy.fa in the working directory, and replace the original one.')
	group2.add_argument('-em', dest='exclude_missing', action='store_true', \
		help='Exclude gap sequences for CNV simulation')
	group2.add_argument('-e_cnv', dest='exon_cnv_list', type=str, default=None, \
		help='A user-defined list of CNVs overlapping with exons')
	group2.add_argument('-e_chr', dest='exon_cnv_chr', type=int, default = None, \
		help='Number of CNVs overlapping with exons to be generated on each chromosome')
	group2.add_argument('-e_tol', dest='exon_cnv_tol', type=int, default = None, \
		help='Total number of CNVs overlapping with exons to be generated across the genome (estimate)')
	group2.add_argument('-e_cl', dest='exon_cnv_len_file', type=str, default=None, \
		help='User supplied file of CNV length for CNVs overlapping with exons')
	group2.add_argument('-o_cnv', dest='out_cnv_list', type=str, default=None, \
		help='A user-defined list of CNVs outside of exons')
	group2.add_argument('-o_chr', dest='out_cnv_chr', type=int, default = None, \
		help='Number of CNVs outside of exons to be generated on each chromosome')
	group2.add_argument('-o_tol', dest='out_cnv_tol', type=int, default = None, \
		help='Total number of CNVs outside of exons to be generated across the genome (estimate)')
	group2.add_argument('-o_cl', dest='out_cnv_len_file', type=str, default=None, \
		help='User supplied file of CNV length for CNVs outside of exons')
	group2.add_argument('-ol', dest='overlap_bps', type=int, default = 100, \
		help='For each CNV overlapping with exons, number of minimum overlapping bps [100]')
	group2.add_argument('-min_len', dest='cnv_min_length', type=int, default=1000, \
		help='Minimum CNV length [1000]')
	group2.add_argument('-max_len', dest='cnv_max_length', type=int, default=100000, \
		help='Maximum CNV length [100000]')
	group2.add_argument('-min_cn', dest='min_copy_number', type=int, default=2, \
		help='Minimum copy number for insertions [2]')
	group2.add_argument('-max_cn', dest='max_copy_number', type=int, default=10, \
		help='Maximum copy number for insertions [10]')
	group2.add_argument('-p', dest='proportion_ins', type=float, default=0.5, \
		help='Proportion of insertions [0.5]')
	group2.add_argument('-f', dest='min_flanking_len', type=int, default=50, \
		help='Minimum length between each CNV [50]')	
	group2.add_argument('-ms', dest='method_s', choices=['random','uniform','gauss'], default="random", \
		help='Distribution of CNVs [random]')
	group2.add_argument('-ml', dest='method_l', choices=['random','uniform','gauss','beta','user'], default="random", \
		help='Distribution of CNV length [random]')
	group2.add_argument('-as', dest='as1', type=float, default=None, \
		help='Mu for gauss CNV distribution [0]')
	group2.add_argument('-bs', dest='bs', type=float, default=None, \
		help='Sigma for gauss CNV distribution [1]')
	group2.add_argument('-al', dest='al', type=float, default=None, \
		help='Mu (gauss distribution) or alpha (beta distribution) for CNV length distribution [0 for gauss distribution, and 2 for beta distribution]')
	group2.add_argument('-bl', dest='bl', type=float, default=None, \
		help='Sigma (gauss distribution) or beta (beta distribution) for CNV length distribution [1 for gauss distribution, and 2 for beta distribution]')
	group2.add_argument('-s_r', dest='s_rate', type=float, default=0, \
		help='Rate of SNPs in target regions [0]')
	group2.add_argument('-s_s', dest='s_slack', type=int, default=0, \
		help='Slack region up and down stream of target regions to simulate SNPs [0]')
	group2.add_argument('-i_r', dest='i_rate', type=float, default=0, \
		help='Rate of indels in target regions [0]')
	group2.add_argument('-i_mlen', dest='i_max_len', type=float, default=50, \
		help='The Maximum length of indels in target regions [50]. (If a deletion is equal or larger than the length of the target region it is in, the length of the deletion will be changed to (length of the target region it is in) - 1.)')
	
	group3 = parser.add_argument_group('Arguments for simulating short reads (fastq)')
	group3.add_argument('-nr', dest='nreads', type=int, default=10000, \
		help='Number of reads / read pairs on target regions to be generated for each genome [10000]')
	group3.add_argument('-fs', dest='frag_size', type=int, default=200, \
		help='Mean fragment size to be generated [200]')
	group3.add_argument('-s', dest='stdev', type=int, default=20, \
		help='Standard deviation of fragment sizes [20]')
	group3.add_argument('-l', dest='read_length', type=int, default=100, \
		help='Read length of each short read [100]')
	group3.add_argument('-tf', dest='target_region_flank', type=int, default=0, \
		help='Length of flanking region up and down stream of target regions to be sequenced (this step take place after -clr). [0]')
	group3.add_argument('-pr', dest='paired_end', action='store_true', \
		help='Select if paired-end sequencing')
	group3.add_argument('-q', dest='quality_score_offset', type=int, default=33, \
		help='Quality score offset for short reads simulation [33]')
	group3.add_argument('-clr', dest='connect_len_between_regions', type=int, default=None, \
		help='Maximum length bwtween target regions to connect the target regions.')
	group3.add_argument('-m', dest='model', type=str, \
		default=os.path.join(os.path.dirname(os.path.realpath(__file__)) + "/ill100v5_p.gzip"), \
		help='GemSim error model file (.gzip, need absolute path) [ill100v5_p]')

	group4 = parser.add_argument_group('Arguments for general settings')
	group4.add_argument('-o', dest='output_dir', type=str, default="simulation_output", \
		help='Output directory [simulator_output]')
	group4.add_argument('-rn', dest='rearranged_output_name', type=str, default="test", \
		help='Prefix of the rearranged outputs (do not include directory name) [test]')
	group4.add_argument('-n', dest='num_samples', type=int, default=1, \
		help='Number of test samples to be generated [1]')
	group4.add_argument('-sc', dest='sim_control', action='store_true', \
		help='Simulation for control genome')
	group4.add_argument('-ssr', dest='sim_short_reads', action='store_true', \
		help='Simulate short reads (fastq) files')
	group4.add_argument('-sb', dest='sim_bam', action='store_true', \
		help='Simulate bam files')
	group4.add_argument('-picard', dest='path_to_picard', type=str, default=None, \
		help='Absolute path to picard')
	group4.add_argument('-GATK', dest='path_to_GATK', type=str, default=None, \
		help='Absolute path to GATK')

	args = parser.parse_args()
	
	if not os.path.exists(args.genome_file):
		log_print('Error: The reference genome file does not exist!')
		exit(1)
	if not os.path.exists(args.target_region_file):
		log_print('Error: The target region file does not exist!')
		exit(1)

	param = {}
	param['type'] = 'e'
	param['genome_file'] = os.path.join(os.getcwd(), args.genome_file)
	if args.target_region_file:
		param['target_region_file'] = os.path.join(os.getcwd(), args.target_region_file)
	param['replaceN'] = args.replaceNs
	param['cnv_min_len'] = args.cnv_min_length
	param['cnv_max_len'] = args.cnv_max_length
	param['min_cn'] = args.min_copy_number
	param['max_cn'] = args.max_copy_number
	param['p_ins'] = args.proportion_ins
	param['e_cnv_list'] = args.exon_cnv_list
	param['o_cnv_list'] = args.out_cnv_list
	param['out_dir'] = os.path.join(os.getcwd(), args.output_dir)
	param['e_cnv_chr'] = args.exon_cnv_chr
	param['e_cnv_tol'] = args.exon_cnv_tol
	param['o_cnv_chr'] = args.out_cnv_chr
	param['o_cnv_tol'] = args.out_cnv_tol
	param['overlap_bp'] = args.overlap_bps
	param['tmp_dir'] = os.path.join(param['out_dir'], "tmp")
	#param['rearranged_out'] = args.rearranged_output_name
	param['nreads'] = args.nreads
	param['frag_size'] = args.frag_size
	param['stdev'] = args.stdev
	param['read_length'] = args.read_length
	param['paired_end'] = args.paired_end
	param['qual'] = args.quality_score_offset
	param['model'] = args.model
	#param['sim_control'] = args.sim_control
	param['sim_short_reads'] = args.sim_short_reads
	param['sim_bam'] = args.sim_bam
	param['path_to_picard'] = args.path_to_picard
	param['path_to_GATK'] = args.path_to_GATK
	param['method_s'] = args.method_s
	param['method_l'] = args.method_l
	param['e_cnv_len_file'] = args.exon_cnv_len_file
	param['o_cnv_len_file'] = args.out_cnv_len_file
	param['opt'] = args.exclude_missing
	param['flank'] = args.min_flanking_len
	param['fl'] = args.target_region_flank
	param['inter'] = args.connect_len_between_regions
	param['s_rate'] = args.s_rate
	param['i_rate'] = args.i_rate
	param['i_mlen'] = args.i_max_len
	param['as'] = args.as1
	param['bs'] = args.bs
	param['al'] = args.al
	param['bl'] = args.bl
	param['snp_slack'] = args.s_slack

	t = args.num_samples
	if t < 1:
		log_print("Error: The number of test samples (-n) must be at least 1!")
		exit(1)

	if (param['p_ins'] < 0) or (param['p_ins'] > 1):
		log_print("Error: Insertion rate must be between 0 and 1.")
		exit(1)

	if (param['s_rate'] < 0) or (param['s_rate'] > 1):
		log_print("Error: SNP rate must be between 0 and 1.")
		exit(1)

	if (param['i_rate'] < 0) or (param['i_rate'] > 1):
		log_print("Error: indel rate must be between 0 and 1.")
		exit(1)

	if (param['i_mlen'] < 0):
		log_print("Error: the maximium length of indels must > 0.")
		exit(1)

	if (param['snp_slack'] < 0):
		log_print("Error: the slack region for making SNPS must > 0.")
		exit(1)

	if param['method_s'] == 'gauss':
		if not param['as']:
			param['as'] = 0
		if not param['bs']:
			param['bs'] = 1

	if param['method_l'] == 'gauss':
		if not param['al']:
			param['al'] = 0
		if not param['bl']:
			param['bl'] = 1

	if (param['method_l'] == 'beta'):
		if not param['al']:
			param['al'] = 2
		if not param['bl']:
			param['bl'] = 2
		if (param['al']) and (param['al'] <= 0):
			log_print("Error: alpha must > 0 for beta distribution.")
			exit(1)
		if (param['bl']) and (param['bl'] <= 0):
			log_print("Error: beta must > 0 for beta distribution.")
			exit(1)

	if (param['method_s'] != 'gauss'):
		if param['as'] or param['bs']:
			log_print("Warning: parameters mu and sigma for CNV distribution are not used! (Only used in gauss distribution.)")

	if (param['method_l'] != 'gauss') and (param['method_l'] != 'beta'):
		if param['al'] or param['bl']:
			log_print("Warning: parameters mu and sigma / alpha and beta for CNV length distribution are not used! (Only used in gauss and beta distribution.)")

	if param['sim_bam']:
		if (not param['path_to_picard']) or (not param['path_to_GATK']):
			log_print('Error: Must provide path to picard (-picard) and path to GATK (-GATK)!')
			exit(1)

	if param['sim_short_reads'] and not param['paired_end']:
		log_print("Warning: Chose single-end sequencing. Mean fragment size (-fs) and standard deviation of fragment size (-s) will be ignored.")

	if param['type'] == 'e':	
		e_ct = 0
		if param['e_cnv_list']:
			e_ct += 1
		if param['e_cnv_chr']:
			e_ct += 1
		if param['e_cnv_tol']:
			e_ct += 1
		if param['e_cnv_len_file']:
			e_ct += 1
		if e_ct != 1:
			log_print('Error: One and only one of -e_cnv, -e_chr, -e_tol and -e_cl must be present!')
			exit(1)

		o_ct = 0
		if param['o_cnv_list']:
			o_ct += 1
		if param['o_cnv_chr']:
			o_ct += 1
		if param['o_cnv_tol']:
			o_ct += 1
		if param['o_cnv_len_file']:
			o_ct += 1
		if not (o_ct == 0 or o_ct ==1):
			log_print('Error: Only one of -o_cnv, -o_chr, -o_tol and -o_cl can be present!')
			exit(1)

		if param['e_cnv_list']:
			log_print('Warning: A list of CNVs overlapping with exons are provided. -em, -f, -ms, -ml, -ol, -min_cn, -max_cn, -min_len and -max_len will be ignored for CNVs on this list!')
		if param['o_cnv_list']:
			log_print('Warning: A list of CNVs outside of exons are provided. -em, -f, -ms, -ml, -ol, -min_cn, -max_cn, -min_len and -max_len will be ignored for CNVs on this list!')


		if param['method_l'] == 'user':
			log_print('Warning: -min_len and -max_len will be ignored since "-ml user" is chosen!')
			if not param['e_cnv_len_file']:
				log_print('Error: "-ml user" must be used with -e_cl!')
				exit(1)
			if o_ct == 1 and not param['o_cnv_len_file']:
				log_print('Error: If CNVs outside of exons are to be generated, "-ml user" must be used with -o_cl!')
				exit(1)
		else:
			if param['e_cnv_len_file']:
				log_print('Error: Only "-ml user" could be used with -e_cl!')
				exit(1)
			if o_ct == 1 and param['o_cnv_len_file']:
				log_print('Error: Only "-ml user" could be used with -o_cl!')
				exit(1)

	if param['sim_bam']:
		if not param['sim_short_reads']:
			log_print('Error: Must simulate short reads (-ssr) to simulate bam files!')
			exit(1)

	if os.path.exists(param['tmp_dir']):
		subprocess.call(['rm', '-rf', param['tmp_dir']], stderr=None)
		#shutil.rmtree(param['tmp_dir'])
		os.makedirs(param['tmp_dir'])
	else:      
		os.makedirs(param['tmp_dir'])	

	print '    ==================== SECNVs ====================    '
	sys.stdout.flush()
	print '                      SECNVs (2019)                     '
	sys.stdout.flush()
	print '                     Version 2.4  (July 2019)                    '
	sys.stdout.flush()
	print '        Bug report: Yue Xing <yue.july.xing@gmail.com>        '
	sys.stdout.flush()
	print '    ------------------------------------------------------    '
	sys.stdout.flush()
    
	log_print('Reading genome file...')
	if param['replaceN'] == True:
		subprocess.call(['cp', param['genome_file'], 'genome_file_copy.fa'], stderr=None)
		listi=list("ATGC")
		log_print("Replace gap regions...")
		log_print("Opening the original fasta file...")
		file = open(param['genome_file'])
		log_print("Writing the imputed fasta file...")
		file_write = open("tmp.fa", 'w')
		while True:
			line = file.readline().rstrip('\n')
			if not line:
				break
			if line[0]=='>':
				file_write.writelines(line)
				file_write.write('\n')
			else:	
				seq=list(line)
				for i in range(len(seq)):
					if (seq[i] == "N") or (seq[i] == "n"):
						seq[i] = random.choice(listi)
				seq2 = ''.join(seq)
				file_write.writelines(seq2)
				file_write.write('\n')
		file_write.close()
		file.close()
		subprocess.call(['rm', '-f', param['genome_file']], stderr=None)
		subprocess.call(['mv', 'tmp.fa', param['genome_file']], stderr=None)

	iin_seqs, iin_chrs = read_fasta(param['genome_file'])

	if param['type'] == 'e':
		log_print('Reading target region file...')
		iin_st, iin_ed = read_target(param['target_region_file'], iin_chrs)

	if t == 1:
		param['rearranged_out'] = args.rearranged_output_name
	else:
		log_print('Processing the 1st sample and control (if required)...')
		param['rearranged_out'] = args.rearranged_output_name + "1"

	simulate_WES(param, iin_seqs, iin_chrs, iin_st, iin_ed, args.sim_control, 0)

	if t > 1:
	 	for i in range(1,t):
	 		mess = 'Processing the ' + str(i+1) + 'th sample...'
	 		log_print(mess)
	 		param['rearranged_out'] = args.rearranged_output_name + str(i+1)
			simulate_WES(param, iin_seqs, iin_chrs, iin_st, iin_ed, None, 1)
	
	#shutil.rmtree(param['tmp_dir'])
	subprocess.call(['rm', '-rf', param['tmp_dir']], stderr=None)

	log_print('Done')

if __name__ == '__main__':
	main()