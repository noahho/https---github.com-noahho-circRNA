## This module finds positive and negative data for the task
## Positive training data: a single gene that is converted to circRNA
## hsa_hg19_Rybak2015.bed stores the loci of circRNA genes, we will read the sequence of those genes and use as positive data
## Negative training data: any sequence of bases 
## we will use parts of the exon, that do not contain a whole circRNA gene as negative test data

import re
from intervaltree import Interval, IntervalTree
from collections import namedtuple
from sklearn import svm
from parse import *
import random
import itertools
import time
import json
from bintrees import AVLTree
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

GENOME_FILE = '../data/hg19.fa';
CIRC_FILE = '../data/hsa_hg19_Rybak2015.bed';
EXON_FILE = '../data/all_exons.bed';
TESTMODE = True;

Gene = namedtuple('Gene', ['start', 'end', 'dir'])

# Reads the Data in the genome file and returns a dictionary with an entry for each chromosome
def read_genome_sequence(file, testmode):
    fasta = SeqIO.parse(file, "fasta")

    seqs = {}

    for record in fasta:
        # Make all upper case
		seqs[record.id] = record.seq.upper()
		print("%s sequence loaded" % record.id)
		if (testmode): break;
    fasta.close()

    return seqs
	
# Reads a bed file and returns an array of genes per chromosome
def read_bed(file, testmode):
	bedFile = open(file, 'r')
	genes = {}

	for line in bedFile:
		vals = line.split()
		chr_n, start, end = vals[0:3]

        # In Testmode only the first chromosome is used
		if testmode and chr_n != 'chr1':
			continue

        # Short BED
		if vals[3] in ['+','-']:
			dir = vals[3]
			gene_name = ""
        # Long BED
		else:
			dir = vals[5]
			gene_name = vals[3]

		start, end = int(start), int(end)

		# Omit weird data
		if end <= start:
			continue

		if not chr_n in genes:
			genes[chr_n] = []
		
		genes[chr_n].append(Gene(int(start), int(end), dir))

	bedFile.close()
	return genes

# Returns the sequence of the given gene
def get_gene_data(gene, chr, seqs):
    data = seqs[chr][gene.start:gene.end]

    if gene.dir == '-':
        data = data.reverse_complement();

    return SeqRecord(data, description=chr)

# hsa_hg19_Rybak2015.bed stores the loci of circRNA genes, we will read the sequence of those genes and use as positive data
# we will use parts of the exon, that do not contain a whole circRNA gene as negative test data
def generate_examples(testmode):
	print("Loading files...")
	seqs = read_genome_sequence(GENOME_FILE, TESTMODE) # A dictionary with an entry for each chromosome
	crnas = read_bed(CIRC_FILE, TESTMODE) # A dictionary with an array of genes per chromosome
	exons = read_bed(EXON_FILE, TESTMODE) # A dictionary with an array of genes per chromosome
	notcrnas = [] # A dictionary with an array of genes per chromosome, that will be filled later
	print("## Loaded "+str(len(crnas))+" CircRNA Chromosomes")	
	print("## Loaded "+str(len(exons))+" Exon Chromosomes")
	
	print("Loading data structures...")
	
	# Setting up a tree to check if a point is contained in a crna
	per_chr_crna_tree = {}
	per_chr_crna_dict = {}
	
	# Building an AVL Tree of starts and ends of the circRNAs for every chromosome
	for nchr in crnas:
		if not nchr in per_chr_crna_dict:
			per_chr_crna_dict[nchr] = {}
		for gene in crnas[nchr]:
			per_chr_crna_dict[nchr][gene.start] = gene.end
		per_chr_crna_tree[nchr] = AVLTree(per_chr_crna_dict[nchr])
	
	print("Generating Negative Examples...")
	# test if an exon contains a circRNA gene or is contained in a circRNA gene
	per_chr_not_circ_gene_array = {};
	ngene = 0;
	for nchr in exons:
		if nchr in per_chr_crna_tree:
			crna_tree = per_chr_crna_tree[nchr]
		else:
			per_chr_not_circ_gene_array[nchr] = exons[nchr];
			continue
		maxlocus = max(crna_tree)
		minlocus = min(crna_tree)
		per_chr_not_circ_gene_array[nchr] = [];
		
		for gene in exons[nchr]:
			ngene = ngene + 1
			# if gene contains circRNA it is not a negative example
			# that is if there is a circRNA that starts after the gene begins and before the gene ends
			circRNA_before_end = crna_tree.floor_item(gene.start)[1] if (gene.start >= minlocus) else 0
			circRNA_after_start = crna_tree.ceiling_item(gene.start)[0] if (gene.start < maxlocus) else maxlocus+1
			if not (circRNA_before_end > gene.start or circRNA_after_start < gene.end):
				per_chr_not_circ_gene_array[nchr].append(gene)
	
	print("Loading the gene information...")
	# Read the dna sequeneces of the positive and negative examples
	positive_examples, negative_examples = [], []
	for nchr in crnas:
		print "## Positive Chromosome "+str(nchr)+" has "+str(len(crnas[nchr]))+ " genes"
		for gene in crnas[nchr]:
			positive_examples.append(get_gene_data(gene, nchr, seqs))
			
	for nchr in per_chr_not_circ_gene_array:
		print "## Negative Chromosome "+str(nchr)+" has "+str(len(per_chr_not_circ_gene_array[nchr]))+ " genes"
		for gene in per_chr_not_circ_gene_array[nchr]:
			negative_examples.append(get_gene_data(gene, nchr, seqs))
			
	print("Saving test data...")
	SeqIO.write(positive_examples, "../data/positive.fa", "fasta")
	SeqIO.write(negative_examples, "../data/negative.fa", "fasta")
generate_examples(False);