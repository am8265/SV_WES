#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
import sys
import os
VCF_HEADER = 'header.vcf.txt'
inputBed = sys.argv[1]
sample = '.'.join(inputBed.split('.')[:2])

def removeDuplicates(infilename, outfilename):
   lines_seen = set() # holds lines already seen
   outfile = open(outfilename, "w")
   for line in open(infilename, "r"):
      if line not in lines_seen: # not a duplicate
         outfile.write(line)
         lines_seen.add(line)
   outfile.close()


with open(VCF_HEADER, 'r') as vcf_fh:
	vcf_header=vcf_fh.readlines()

vcf_header = ''.join(vcf_header)

delVCF_fh = open(sample + '.del.vcf.tmp', 'w')
dupVCF_fh = open(sample + '.dup.vcf.tmp', 'w')

# write the template header to the del and dup vcf, also write ^#CHROM header

chrom_header = '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample])
delVCF_fh.write(vcf_header + '\n' + chrom_header + '\n')
dupVCF_fh.write(vcf_header + '\n' + chrom_header + '\n')



with open(inputBed,'r') as bed:
	next(bed)
	for line in bed:
		bed_fields = line.strip('\n').split('\t')
		chrom = bed_fields[0]
		start = bed_fields[1]
		end = bed_fields[2]
		svlen = int(end) - int(start)
		svtype = bed_fields[3]

		if svtype == "deletion":
			ID = '_'.join(['ED','DEL', chrom, start, end])
			record_start = '\t'.join([chrom, start, ID, 'N', '<DEL>', '.', 'PASS'])
			record_end = ';'.join(['END=' + end, 'SVLEN=-' + str(svlen), 'SVTYPE=' + 'DEL', 'CIEND=-1000,1000', 'CIPOS=-1000,1000'])
			record_final = record_start + '\t' + record_end + '\t' + "GT" + '\t' + "0/1"
			delVCF_fh.write(record_final + '\n')
			

		else:
			ID = '_'.join(['ED','DUP', chrom, start, end])
			record_start = '\t'.join([chrom, start, ID, 'N', '<DUP>', '.', 'PASS'])
			record_end = ';'.join(['END=' + end, 'SVLEN=' + str(svlen), 'SVTYPE=' + 'DUP', 'CIEND=-1000,1000', 'CIPOS=-1000,1000'])
			record_final = record_start + '\t' + record_end + '\t' + "GT" + '\t' + "0/1"
			dupVCF_fh.write(record_final + '\n')

delVCF_fh.close()
dupVCF_fh.close()


removeDuplicates(sample + '.del.vcf.tmp', sample + '.del.vcf')
removeDuplicates(sample + '.dup.vcf.tmp', sample + '.dup.vcf')

os.system('rm ' + sample + '.del.vcf.tmp')
os.system('rm ' + sample + '.dup.vcf.tmp')

