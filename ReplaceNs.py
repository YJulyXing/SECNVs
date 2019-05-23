#!/usr/bin/python

import argparse
import sys
import os
import random
import time

def log_print(message):
	print '[SimulateCNVs ' + time.asctime( time.localtime(time.time()) ) + '] ' + message
	sys.stdout.flush()

parser = argparse.ArgumentParser(description='Impute fasta files', \
	formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i", type=str, dest='input_fasta_file', required=True, help="The original fasta file name")
parser.add_argument("-o", type=str, dest='output_fasta_file', required=True, help="The imputed fasta file name")

args = parser.parse_args()

file_name = args.input_fasta_file
file2_name = args.output_fasta_file

listi=list("ATGC")

if not os.path.exists(file_name):
	log_print("The provided file does not exist!")
	exit(1)
if os.path.exists(file2_name):
	os.remove(file2_name)

log_print("Opening the original fasta file...")
file = open(file_name)

log_print("Writing the imputed fasta file...")
file_write = open(file2_name, 'w')
while True:
	line = file.readline().rstrip('\n')
	if not line:
		break
	seq=list(line)
	for i in range(len(seq)):
		if seq[i] == "N":
			seq[i] = random.choice(listi)
	seq2 = ''.join(seq)
	file_write.writelines(seq2)
	file_write.write('\n')

file_write.close()
file.close()
