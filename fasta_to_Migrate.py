# The script will convert fasta-format sequences to Migrate format based on a ID nums and pop IDs in another file

# Args
## File with list of input fasta paths
## File for ID numbers and populations file. Comma-separated with ID followed by population
## Filepath output file

import sys
import os
import re
import pandas as pd
from Bio import SeqIO

indir = sys.argv[1]
inidnums = sys.argv[2]
outdir = sys.argv[3]

df = pd.read_csv( inidnums , header=None )
all_ids = list( df.iloc[ : , 0 ] )
pops = set( list( df.iloc[ : , 1 ] ) )

fasta_files = open( indir, 'r' )
fasta_files = fasta_files.read().splitlines()
fasta_files = list(fasta_files)

with open('seqs.tmp', "w", newline='') as f:
	for pop in pops:
		pop_df = df[df[1]==pop]		
		idnums = list( pop_df.iloc[ : , 0 ] )
		f.write(str(len(pop_df[1])) + ' ' + pop + '\n')
		for id in idnums:
			f.write(str(id) + '\t')
			for fastafp in fasta_files:
				fasta_seqs = SeqIO.parse(fastafp, 'fasta')
				fasta_handle = open(fastafp, "r")
				if id in fasta_handle.read():
					for seq in fasta_seqs:
						if id in seq.id:
							seq_len = len(seq.seq)
							seq_gaps = seq.seq.count('-')
							seq_gaps = int(float(seq_gaps))
							percent_gaps = (seq_gaps/seq_len)*100
							if percent_gaps > 90:
								print("Warning: Sequence " + id + " from file " + fastafp + " has 90% or more gaps")
							f.write(str(seq.seq) + ' ')
				else:
					fasta_seqs_len = SeqIO.parse(fastafp, 'fasta')
					seq_len = len(next(fasta_seqs_len).seq)	
					seq_str = '-' * seq_len
					f.write(str(seq_str) + ' ')
					print("Warning: Sequence " + id + " could not be found in alignment " + fastafp)
			f.write('\n')

with open('len_header.tmp', "w", newline='') as header:
	for fastafp in fasta_files:
		fasta_seqs = SeqIO.parse(fastafp, 'fasta')
		locus_len = len(next(fasta_seqs).seq)
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
