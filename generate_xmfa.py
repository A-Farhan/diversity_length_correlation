## require
import os, sys
from Bio import AlignIO

## arguments
indr = sys.argv[1] # input directory of codon MSAs
fout = sys.argv[2] # output XMFA alignment file

# list of input files
fls = os.listdir(indr)
nf = len(fls)
print("# MSAs = %i"%nf)

# initialize variables
alns = []
genome_names = []
gene_end = 0
gene_start = 1

for x in range(nf):
	# gene name
	gene = fls[x].replace('.fasta','')
	# input msa file path
	fpath = os.path.join( indr, fls[x])
	
	# biopython alignment
	try:
		aln = AlignIO.read( fpath, 'fasta')
	except ValueError:
		print("Empty file: %s"%fpath)
		continue

	# list of genomes names present in the alignment
	aln_genome_names = [ i.id.split('|')[-1] for i in aln]
	# list of all genome names seen so far
	genome_names = genome_names + [ i for i in aln_genome_names if i not in genome_names]
	# indices of alignment's genomes in the above list
	aln_genome_ixs = [ genome_names.index(i) for i in aln_genome_names]
	# number of sites in the alignment
	l = aln.get_alignment_length()
	# update gene start and end position
	gene_start = gene_end + 1
	gene_end =  gene_start + l - 1
	
	for y,i in enumerate(aln):
		# modified header for a sequence in the alignment
		header = ''.join( [ '>', str( aln_genome_ixs[y]), '__', aln_genome_names[y], ':', 
			str(gene_start), '-', str(gene_end), ' + ', gene])
		seq = str(i.seq)
		alns.extend([ header, seq])

	# end of a gene alignment
	alns.extend('=')

print("# number of genomes = %i"%len(genome_names))

## write output
with open(fout,'w') as flob:
	for i in alns:
		flob.write(i+'\n')
			
print("Program finished!")
