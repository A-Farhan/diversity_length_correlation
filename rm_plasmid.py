### Remove plasmids from genome sequences

import os, sys
from Bio import SeqIO

fin = sys.argv[1]
fout = sys.argv[2]

recs = SeqIO.index(fin,'fasta')

# non-plasmid keys
npk = [ i for i in recs.keys() if 'plasmid' not in recs[i].description]
# corresponding records
subrecs = [ v for k,v in recs.items() if k in npk]
SeqIO.write(subrecs,fout,'fasta')