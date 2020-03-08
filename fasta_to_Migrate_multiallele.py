# The script will convert fasta-format sequences to Migrate format based on a ID nums and pop IDs in another file

## Args in order:
## 1. File with list of input fasta file paths. One per line
## 2. File with ID numbers, populations, and allele IDs. Comma-separated with ID followed by population, followed by allele ID.
##### Assumes fasta ID line contains the ID number and allele ID.
##### Example of .csv file:
##### I29224,ALSYM,seq1
##### I29224,ALSYM,seq2
##### I12715,ALALLO,seq1
##### I12715,ALALLO,seq2
## 3. Desired order of populations. One population ID per line. Must match population IDs in file with ID numbers.
## 4. Filepath with filename of output file
from typing import Any

import sys
import os
import re
import pandas as pd
from Bio import SeqIO

indir = sys.argv[1]
inidnums = sys.argv[2]
poporder = sys.argv[3]
outdir = sys.argv[4]

df = pd.read_csv( inidnums , header=None )
all_ids = list( df.iloc[ : , 0 ] )
#pops = set( list( df.iloc[ : , 1 ] ) )
pops = list(pd.read_csv( poporder , header=None )[0])

fasta_files = open( indir, 'r' )
fasta_files = fasta_files.read().splitlines()
fasta_files = list(fasta_files)

with open('seqs.tmp', "w", newline='') as f:
	for pop in pops:
		singlepop_df = df[df[1]==pop]
		f.write(str(len(singlepop_df[1])) + ' ' + pop + '\n')
		for id in set(list(singlepop_df[0])):
			id_df = singlepop_df[singlepop_df[0] == id].iloc[ : , [0,2] ]
			for allele in list(id_df[2]):
				f.write(str(id) + '\t')
#				f.write(str(id) + '_' + str(allele) '\t') ## Optional if different id per allele are wanted. Turn off previous line.
				for fastafp in fasta_files:
					fasta_seqs = list(SeqIO.parse(fastafp, 'fasta'))
					idx_list = []
					for idx in fasta_seqs:
						idx_list.append(idx.id)
					select_idx: Any = [i for i, s in enumerate(idx_list) if id in s and allele in s]
					if select_idx == []:
						print("Warning: Sequence " + id + " " + allele + " could not be found in alignment " + fastafp)
						seq_len = len(fasta_seqs[0])
						seq_str = '-' * seq_len
						f.write(str(seq_str) + ' ')
					else:
						seq_len = len(fasta_seqs[select_idx[0]].seq)
						seq_gaps = fasta_seqs[select_idx[0]].seq.count('-')
						seq_gaps = int(float(seq_gaps))
						percent_gaps = (seq_gaps/seq_len)*100
						if percent_gaps > 90:
							print("Warning: Sequence " + id + " " + allele + " from file " + fastafp + " has 90% or more gaps")
						f.write(str(fasta_seqs[select_idx[0]].seq) + ' ')
				f.write('\n')

with open('len_header.tmp', "w", newline='') as header:
	for fastafp in fasta_files:
		fasta_seqs = SeqIO.parse(fastafp, 'fasta')
		locus_len = len(next(fasta_seqs))
		header.write('(s' + str(locus_len) + ') ')
	header.write('\n')

filenames = ['len_header.tmp', 'seqs.tmp']
with open(outdir, "w") as outfile:
	outfile.write(str(len(pops)) + ' ' + str(len(fasta_files)) + '\n')
	for fname in filenames:
		with open(fname) as infile:
			outfile.write(infile.read())

for fname in filenames:
	os.remove(fname)   