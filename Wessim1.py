import sys
import random
import bisect
#import pysam
import gzip
import cPickle
import numpy
from time import time, localtime, strftime
import argparse
from multiprocessing import Process
import os
import math
#import pysam

inds={'A':0,'T':1,'G':2,'C':3,'N':4,'a':0,'t':1,'g':2,'c':3,'n':4}

def subprogram(command, name):
	os.system(command)
	print "exiting subprocess " + str(name)

def main(argv):
	t0 = time()
	arguline = " ".join(argv)
	parser = argparse.ArgumentParser(description='Wessim1: Whole Exome Sequencing SIMulator 1 (Ideal target region-based version)', prog='Wessim1', formatter_class=argparse.RawTextHelpFormatter)
	group1 = parser.add_argument_group('Mandatory input files')
	group1.add_argument('-R', metavar = 'FILE', dest='reference', required=True, help='faidx-indexed (R)eference genome FASTA file')
	group1.add_argument('-B', metavar = 'FILE', dest='region', required=True, help='Target region .(B)ED file')

	group2 = parser.add_argument_group('Parameters for exome capture')
	group2.add_argument('-f', metavar = 'INT', type=int, dest='fragsize', required=False, help='mean (f)ragment size. this corresponds to insert size when sequencing in paired-end mode. [200]', default=200)
	group2.add_argument('-d', metavar = 'INT', type=int, dest='fragsd', required=False, help='standard (d)eviation of fragment size [50]', default=50)
	group2.add_argument('-m', metavar = 'INT', type=int, dest='fragmin', required=False, help='(m)inimum fragment length [read_length + 20]')
	group2.add_argument('-x', metavar = 'INT',type=int, dest='slack', required=False, help='slack margin of the given boundaries [0]', default=0) 
	
	group3 = parser.add_argument_group('Parameters for sequencing')
	group3.add_argument('-p', action='store_true', help='generate paired-end reads [single]')
	group3.add_argument('-n', metavar = 'INT', type=int, dest='readnumber', required=True, help='total (n)umber of reads')	
	group3.add_argument('-l', metavar = 'INT', type=int, dest='readlength', required=True, help='read (l)ength (bp)')
	group3.add_argument('-M', metavar = 'FILE', dest='model', required=True, help='GemSim (M)odel file (.gzip)')
	group3.add_argument('-t', metavar = 'INT', type=int, dest='threadnumber', required=False, help='number of (t)hreaded subprocesses [1]', default=1) 

	group4 = parser.add_argument_group('Output options')
	group4.add_argument('-o', metavar = 'FILE', dest='outfile', help='(o)utput file header. ".fastq.gz" or ".fastq" will be attached automatically. Output will be splitted into two files in paired-end mode', required=True)
	group4.add_argument('-z', action='store_true', help='compress output with g(z)ip [false]')
	group4.add_argument('-q', metavar = 'INT', type=int, dest='qualbase', required=False, help='(q)uality score offset [33]', default=33)
	group4.add_argument('-v', action='store_true', help='(v)erbose; print out intermediate messages.')

	args = parser.parse_args()
	reffile = args.reference
	regionfile = args.region
	 
	isize = args.fragsize
	isd = args.fragsd
	imin = args.fragmin
	slack = args.slack

	getRegionVector(reffile, regionfile, slack)
	paired = args.p
	readlength = args.readlength
	readnumber = args.readnumber		
	threadnumber = args.threadnumber
	
	if imin==None:
		if paired:
			imin = readlength + 20
		else:
			imin = readlength + 20
	if isize < imin:
		print "too small mean fragment size (" + str(isize) + ") compared to minimum length (" + str(imin) + "). Increase it and try again."
		sys.exit(0) 
	model = args.model
		
	outfile = args.outfile
	compress = args.z
	qualbase = args.qualbase
	verbose = args.v

	print 
	print "-------------------------------------------"
	print "Reference:", reffile
	print "Region file:", regionfile
	print "Fragment:",isize, "+-", isd, ">", imin
	print "Paired-end mode?", paired
	print "Sequencing model:", model
	print "Read length:", readlength, "Read number:", readnumber
	print "Output File:", outfile
	print "Gzip compress?", compress
	print "Quality base:", qualbase
	print "Thread number:", threadnumber
	print "Job started at:", strftime("%Y-%m-%d %H:%M:%S", localtime())
	print "-------------------------------------------"
	print

	processes = []
	for t in range(0, threadnumber):
		readstart = int(float(readnumber) / float(threadnumber) * t) + 1
		readend = int(float(readnumber) / float(threadnumber) * (t+1))
		command = "python __sub_wessim1.py " + arguline + " -1 " + str(readstart) + " -2 " + str(readend) + " -i " + str(t+1)
		p = Process(target=subprogram, args=(command, t+1))
		p.start()
		processes.append(p)
	for p in processes:
		p.join()
	t1 = time()
	print "Done generating " + str(readnumber) + " reads in %f secs" % (t1 - t0)
	print "Merging subresults..."
	wread = None
	wread2 = None
	if paired and compress:
		wread = gzip.open(outfile + "_1.fastq.gz", 'wb')
		wread2 = gzip.open(outfile + "_2.fastq.gz", 'wb')
	elif paired and not compress:
		wread = open(outfile + "_1.fastq", 'w')
		wread2 = open(outfile + "_2.fastq", 'w')
	elif not paired and compress:
		wread = gzip.open(outfile + ".fastq.gz", 'wb')
	else:
		wread = open(outfile + ".fastq", 'w')
	if not paired:
		for t in range(0, threadnumber):
			suboutfile = outfile + "-" + str(t+1)
			fread = None
			if compress:
				suboutfile += ".fastq.gz"
				fread = gzip.open(suboutfile, 'rb')
			else:
				suboutfile += ".fastq"
				fread = open(suboutfile, 'r')
			line = fread.readline()
			while line:
				wread.write(line)
				line = fread.readline()
			fread.close()
			os.remove(suboutfile)
		wread.close()
	else:
		for t in range(0, threadnumber):
			suboutfile1 = outfile + "-" + str(t+1) + "_1"
			suboutfile2 = outfile + "-" + str(t+1) + "_2"
			fread1 = None
			fread2 = None
			if compress:
				suboutfile1 += ".fastq.gz"
				suboutfile2 += ".fastq.gz"
				fread1 = gzip.open(suboutfile1, "rb")
				fread2 = gzip.open(suboutfile2, "rb")
			else:
				suboutfile1 += ".fastq"
				suboutfile2 += ".fastq"
				fread1 = open(suboutfile1, "r")
				fread2 = open(suboutfile2, "r")
			line1 = fread1.readline()
			line2 = fread2.readline()
			while line1 and line2:
				wread.write(line1)
				wread2.write(line2)
				line1 = fread1.readline()
				line2 = fread2.readline()
			fread1.close()
			fread2.close()
			os.remove(suboutfile1)
			os.remove(suboutfile2)
		wread.close()
		wread2.close()
	sys.exit(0)

def read_fasta(fname):
	with open(fname, "r") as fh:
		name = None
		seqs = {}
		seqs2 = {}
		for line in fh.readlines():
			if line[0]=='>':
				line = line.rstrip()
				name = line[1:]
				name = name.split()[0]
				name = name.split("_")[0]
				seqs[name] = []
			else:
				seqs[name].append(line.rstrip())
	#chrs = []
	for key in seqs:
		#chrs.append(key)
		seqs2[key] = ''.join(seqs[key])
	fh.close()
	return(seqs2)

def getRegionVector(fastafile, regionfile, slack):
	print("Generating fasta file for given regions...")
	faoutfile = regionfile + ".fa"
	abdoutfile = regionfile + ".abd"
	ref=read_fasta(fastafile)
	f = open(regionfile)
	wfa = open(faoutfile, 'w')
	wabd = open(abdoutfile, 'w')
	abd = 0
	for i in f.readlines():
		values = i.split("\t")
		if i.startswith("#") or len(values)<3:
			continue
		chrom = values[0]
		start = max(int(values[1]) - slack, 1)
		end = int(values[2]) + slack
		header = ">" + chrom + "_" + str(start) + "_" + str(end)
		x0 = ref[chrom]
		x = x0[(start-1):end]
		length = len(x)
		abd += length
		wfa.write(header + "\n")
		wfa.write(x + "\n")
		wabd.write(str(abd) + "\n")
	f.close()
	wfa.close()
	wabd.close()
	
if __name__=="__main__":
	main(sys.argv[1:])

# fewer dependencies
# can read in first line of bed
# fix "length is one less" issue: started from (start+1)
