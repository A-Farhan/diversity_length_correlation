### calculate open state probability of mRNAs of a given genome 
## require
import os, sys, RNA, pandas, multiprocessing
from Bio import SeqIO, AlignIO

## arguments
fft = sys.argv[1] # full path for feature table
fgen = sys.argv[2] # full path for genome fasta
fout = sys.argv[3] # full path to output table
ncpu = int(sys.argv[4]) # no. of cpu threads
up = int(sys.argv[5]) # length of region upstream of translation start
down = int(sys.argv[6]) # length of region downstream to translation star

## functions
def pub(bppm):
	"""Return list of per-base probabilities of being unpaired in a secondary rna structure."""
	L = len(bppm)
	out = [ round( 1 - ( sum(bppm[i][j] for j in range(i+1,L)) + sum(bppm[k][i] for k in range(i))),3) 
	for i in range(L)]
	return out

def extract_seq_tls(gl,gr,ot,chrom,up,down):
	"""Extract sequence of a gene from a chromosome with upstream region."""
	# check that strand is either + or -
	assert ot in ['+','-']

	# find left and right slicing index
	if ot == '+':
		xl,xr = gl-up, gl+down
	else:
		xl,xr = gr-down, gr+up
	
	# check that both ends are within genome
	if xl < 0 or xr >= len(chrom):
		print("Boundary gene! Returning None..")
		return

	seq = chrom[xl:xr].seq
	out = seq if ot == '+' else seq.reverse_complement()
	return out

def mainfun(gl,gr,ot,chrom):
	# extract sequence from the chromosome, alongwith upstream region
	gene_seq = str( extract_seq_tls(gl=gl,gr=gr,ot=ot,chrom=chrom,up=up,down=down))
	# create object to evaluate rna secondary structures
	fc = RNA.fold_compound(gene_seq)
	# predict Minimum Free Energy and corresponding secondary structure
	(mfe_struct, mfe) = fc.mfe()
	# rescale boltzmann factors for partition function computation
	fc.exp_params_rescale(mfe)
	# compute partition function
	(pp, pf) = fc.pf()
	# base-pair probability matrix
	mat = fc.bpp()
	# prob for each site to be unpaired
	probs = pub(mat)
	return probs

##############################
# feature table
ft = pandas.read_csv(fft,sep='\t')
print("# features = %i"%(len(ft)))
# subset rows to CDS excluding pseudo-genes and genes from plasmid(s) 
sft = ft[( ft.iloc[:,0]=='CDS') & (ft['class']=='with_protein')]
print("# CDS, with_protein, chromosomal features = %i"%(len(sft)))
# subset to genes longer than the size of downstream part
sft = sft[ sft['feature_interval_length'] > down]
print("# genes satisfying length threshold = %i"%(len(sft)))

# subset columns
sft = sft.loc[:,['symbol','start','end','strand','genomic_accession']]
sft.index = range(len(sft))
# assembly seqrecords
asm = SeqIO.to_dict( SeqIO.parse(fgen, 'fasta'))
## main
with multiprocessing.Pool(ncpu) as mp_pool:
	ptab = mp_pool.starmap( mainfun, [ ( *sft.iloc[x,1:4], asm[sft.iloc[x,4]]) 
		for x in range(sft.shape[0])])  

out = pandas.concat( [ sft.symbol, pandas.DataFrame(ptab)], axis=1)
out.to_csv(fout,index=False)
print("Program finished!")
