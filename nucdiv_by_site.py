### Get nucleotide diversity from a MSA give a list of sites

## require
import os, sys, numpy, pandas, multiprocessing
from Bio import AlignIO
from utils._divs import pwd_ar1d

## arguments
daln = sys.argv[1] # directory of codon MSAs
ddeg = sys.argv[2] # dir of site degeneracy
dout = sys.argv[3] # output dir of site diversity tables
alnsfx = sys.argv[4] # alignment file suffix
gapfreq = float(sys.argv[5]) # max allowed freq of gaps
ncpu = int(sys.argv[6]) # number of CPUs

# create output dir if needed
not os.path.exists(dout) and os.mkdir(dout)

## Function ##
def get_div( faln, fdsq, fout, gf = gapfreq, alpha = set('ACGT')):
    # load alignment
    try:
        aln = AlignIO.read(faln, 'fasta')
    except ValueError:
        print("No alignment in %s!"%faln)
        return
        
    # numpy array form of the MSA
    ar = numpy.array( [ list(rec) for rec in aln]) 
    # number of sequences and number of sites
    N,L = ar.shape
    # indices of columns where the element exceeds frequency threshold
    fc = gf*N
    gixs = [ x for x in range(L) if list(ar[:,x]).count('-') > fc ]
    # no. of pairs
    numpairs = N*(N-1)/2
    # div for each site
    sdiv = [ round( pwd_ar1d(ar[:,k],alpha)/numpairs, 4) if k not in gixs else 2 for k in range(L)]    
    
    # degeneracy
    with open(fdsq) as flob:
        dseq = flob.read().rstrip('\n')
    
    # write output table
    out = pandas.DataFrame([ (x,int(i),sdiv[x]) for x,i in enumerate(dseq)])
    out.to_csv(fout,index=False,header=['site','degeneracy','diversity'])

## Main
# names of input files
fnames = os.listdir(daln)
# tuples of all corresponding IO file paths
fls = [ ( os.path.join( daln, i), \
          os.path.join( ddeg, i.replace(alnsfx,'.txt')), \
          os.path.join( dout, i.replace(alnsfx,'.csv'))) for i in fnames]

with multiprocessing.Pool(ncpu) as mp_pool:
	mp_pool.starmap( get_div, fls)
