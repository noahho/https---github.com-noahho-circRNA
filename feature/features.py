import sys
import numpy as np
import gzip
import random
import os
import bz2
import cPickle
import pdb
import urllib
import argparse
import subprocess
import itertools

SYSTEM = "WIN"
POSITIVE_EXAMPLES_FILE = 'data/positive.fa';
NEGATIVE_EXAMPLES_FILE = 'data/negative.fa';

chromosome_list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12","chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chrX", "chrY", "chrUn_gl000220"]
	
def get_features(data):
	print 'extracting repeat feature...'
	#repeat_features = get_tandem_feature()
	#print repeat_features

	print 'extracting length and chromosome feature...'
	length_features = get_other_features(data)
	
	print 'extracting kmer feature...'
	kmer_features = get_kmer_features(data, 4)
	
	i = 0
	features = []
	for gene in data:
		features.append(length_features[i] + kmer_features[i])
		i = i+1
		
	return features
def get_other_features(data):
	other_features = []
	
	for gene in data:
		chrom_features = []
		for chr in chromosome_list:
			if chr == gene[1]:
				chrom_features.append(1);
			else:
				chrom_features.append(0);
		if chr in chromosome_list:
			chrom_features.append(0);
		else:
			chrom_features.append(1);
		other_features.append([len(gene[0])] + chrom_features);
		
	
	return other_features

def get_kmer_features(data, k):
	key_kmers = []
	for i in range(1, k+1):
		key_kmers += [ "".join(c) for c in itertools.combinations_with_replacement('ACGT', i) ]
	key_kmers += ["GTAG"];

	kmer_features = []
	i = 0;
	p1 = (len(data)/10);
	for gene in data:
		gene = gene[0]
		if (i % p1 == 0):
			print "## 1% step..."
		kmer_features.append([ float(gene.count(kmer)) / len(gene)
			for kmer in key_kmers ])
		i += 1;
	return kmer_features
	
def get_tandem_feature():
	feature = get_tandem_feature_from_file(POSITIVE_EXAMPLES_FILE)
	feature += get_tandem_feature_from_file(NEGATIVE_EXAMPLES_FILE)
	
def get_tandem_feature_from_file(file):
    run_trftandem_cmd(file)
    repeat_fea = []
    repeat_file = file + '_repeat'
    fp1 = open(repeat_file, 'r')
    for line in fp1:
        values = line.split()
        repeat_fea.append(float(values[0]))
        
    fp1.close()
    os.remove(repeat_file)
	
    return repeat_fea
def parse_tandem_feature(tandem_file, out_file):
    fp = open(tandem_file, 'r')

    parameter_flag = False
    fw = open(out_file + '_tmp', 'w')
    for line in fp:
        if line.find('Sequence:') != -1:
            fw.write('>%s' %line)
            parameter_flag = False
            continue

        if line.find('Parameters:') != -1:
            parameter_flag = True
            continue

        if parameter_flag and len(line) > 5 :
            fw.write('%s' %line)
    fp.close()
    fw.close()
    fw1 = open(out_file, 'w')
    fp1 = open(out_file + '_tmp', 'r')
    first_flag = True
    tmp = [0] * 12
    for line in fp1:
        values = line.split()
        if line[0] == '>':
            if not first_flag:
                fw1.write('%d\t' %tmp[0])
                if tmp[0] > 0:
                    for val in tmp[1:]:
                        fw1.write('%0.3f\t' %(float(val)/tmp[0]))
                else:
                    for val in tmp[1:]:
                        fw1.write('%0.3f\t' %val)
                fw1.write('\n')        
            seq_name = values[1]
            tmp = [0] * 12
            #fw1.write('%s:\t' %seq_name)
            first_flag = False
        else:
            #print line
            ext_vals = values[2:13]
            #print len(ext_vals)
            tmp[0] = tmp[0] + 1
            for index in range(11):
                tmp[index + 1] = tmp[index + 1] + float(ext_vals[index])
    
    fw1.write('%d\t' %tmp[0])
    if tmp[0] > 0:
        for val in tmp[1:]:
            fw1.write('%0.3f\t' %(float(val)/tmp[0]))
    else:
        for val in tmp[1:]:
            fw1.write('%0.3f\t' %val)
    fw1.write('\n')     
    fp1.close()  
    fw1.close()
	
    os.remove(out_file + '_tmp')
                
	
def run_trftandem_cmd(fasta_file):
	if (SYSTEM == "LIN"):
		cli_str = 'trf404.linux64 ' + fasta_file + ' 2 7 7 80 10 50 500 -f -d -h'
	else :
		cli_str = 'trf409.dos64 ' + fasta_file + ' 2 7 7 80 10 50 500 -f -d -h'
	print cli_str;
	FNULL = open(os.devnull, 'w')
	subprocess.call(cli_str, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
	#filelist = [ f for f in os.listdir("../data/") if f.endswith(".html") ]
	#for f in filelist:
	#	os.remove(f)     
    
	trf_out_file = fasta_file + '.2.7.7.80.10.50.500.dat'
	out_file = fasta_file + '_repeat'
	features = parse_tandem_feature(trf_out_file, out_file)

	os.remove(trf_out_file)
	
	return features