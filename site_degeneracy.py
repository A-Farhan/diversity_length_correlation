### Get degeneracy of each site from the consensus of a given codon MSA

## require
import os, sys, multiprocessing
from Bio import AlignIO
from Bio.Data.CodonTable import unambiguous_dna_by_id as codon_table
from collections import Counter
from utils.bioutils import consensus, n2c

## arguments
din = sys.argv[1] # input directory of codon MSAs
dout = sys.argv[2] # output directory of site-degeneracy tables
ncpu = int(sys.argv[3]) # number of CPUs

# create output dir if needed
not os.path.exists(dout) and os.mkdir(dout)

## Function ##
## Function to get degeneracy of each site in a codon
def sitedegen( ct = codon_table[11], nucs = 'ACGT'):
    """
    Return Synonymous, Non-synonymous & Four-fold degenerate site count, in that order,
    for all codons according to CodonTable(Bacteria).

    Arguments:
    - ct 	- codon table [Bacteria]
    - nucs 	- string of allowed nucleotides ['ACGT']
    """
    ixs = range(3)
    table = ct.forward_table
    stops = ct.stop_codons
    codons = table.keys()
    out = { k:() for k in codons}

    # for each codon in the table
    for codon in codons:
        # corresponding amino acid
        aa = table[codon]
        # initialize list to hold fraction of synonymous changes at each site of a codon with zeros
        f = [ 0 for _ in ixs]
        # for every site in the codon
        for i in ixs:
            # produce all possible codons modifying this site
            altc = [ ''.join( [ n if x == i else codon[x] for x in ixs]) for n in nucs]
            # corresponding amino acids for above codons
            alta = [ table[c] for c in altc if c not in stops]
            # count no. of unique amino acid possible
            #nua = len(set(alta))
            # degeneracy
            f[i] = Counter(alta).most_common(1)[0][1]
        out[codon] = f
    return out

# dict of codons and their degeneracy for each site
degen_tab = sitedegen()

def get_dseq( fin, fout, dtab=degen_tab):
	# check if output already present
    if os.path.exists(fout):
        print("Output: %s exists! Moving to next.."%fout)

    # load alignment
    try:
        aln = AlignIO.read(fin, 'fasta')
    except ValueError:
        print("No alignment in %s!"%fin)
        return

    # consensus sequence
    conseq = consensus(aln, minc=0.6)
    # split into codons
    codons = n2c(conseq.seq)
    dlist = [ dtab[i] if i in dtab.keys() else list('000') for i in codons]
    dseq = ''.join( str(i) for x in dlist for i in x)
    with open(fout,'w') as flob:
    	flob.write(dseq+'\n')
    return

## Main
fnames = os.listdir(din)
fins = [ os.path.join( din, i) for i in fnames]
fouts = [ os.path.join( dout, i.replace('.fasta','.txt')) for i in fnames]
fls = list(zip(*[fins,fouts]))

with multiprocessing.Pool(ncpu) as mp_pool:
	mp_pool.starmap( get_dseq, fls)
