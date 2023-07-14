## generate multi-FASTA files for cds and protein sequence records of a gene across genomes

# require
import os, sys, re, click, pandas
import multiprocessing as mp
from collections import Counter
from utils.basic import *
from utils.bioutils import extractgenename
from Bio import SeqIO, SeqRecord

# start clock
start = time()

# regular expressions
rcloc = re.compile('\[locus_tag=(.+?)\]')

## Functions ##
def list_seqrecs(strain,ot,cdr,pdr): 
    """
    Make a list of CDS and proteome sequence records of given protein ids for a given strain.

    Arguments:
    - strain - strain name in the ortholog table
    - ot 	 - ortholog table
    
	note:	this function doesn't make any checks; if the input files are not present, it throws an error
	else it continues and processes all records.
    """                                    
    print("%s: for %s\n"%(timer(start),strain))
    # extract corresponding column of protein ids
    ids = ot[strain]    
    
    ## NUCLEOTIDE
    # input FASTA file of coding sequences
    fin = os.path.join( cdr, strain)
    # dict of sequence records
    cds_recs = SeqIO.to_dict( SeqIO.parse(fin,'fasta'))
    
    # modify CDS records' id
    for k,v in cds_recs.items():
        # extract gene name from the record description
        gene = extractgenename(v.description)
        if gene is None:
            gene = ''
        # extract protein id from records' name
        pid = re.search('.*?_cds_?(.*)_[0-9]+$',v.name).group(1)
        # extract chromosome accession from record's name
        chrom = re.split('\||_cds_',v.name)[1]
        # generate new id using the above 
        v.id = '|'.join([pid,gene,chrom,strain])
    
    # change keys to proteins ids
    cds_recs = { re.search('.*?_cds_?(.*)_[0-9]+$',v.name).group(1):v for k,v in cds_recs.items()}
    # records corresponding to ids
    out1 = [ cds_recs[i] if i in cds_recs.keys() else None for i in ids ] ## CDS

    ## AMINO ACID
    # input FASTA file of protein sequences
    fin = os.path.join( pdr, strain)
    # dict of sequence records
    prot_recs = SeqIO.to_dict( SeqIO.parse(fin,'fasta'))
    # initial output list of protein records mapped from ortholog table, with modified ids
    out2 = []
    
    # for every protein id present in the ortholog table, from this genome 
    for i in ids:

        # if the id is not present among CDS records
        if i not in cds_recs.keys():
            out2.append(None)
            continue

        # if id is present in the proteome FASTA file
        if i in prot_recs.keys():
            # extract the record
            rec = prot_recs[i]
            # change protein records id to the corresponding CDS records id, that was modified above
            rec.id = cds_recs[i].id
            out2.append(rec)
        else:
            out2.append(None)

    return {'cds':out1,'prots':out2}

@click.command()
@click.option('--fin', help='full path to table of ortholog groups', type=click.Path(), required=True)
@click.option('--cidr', help='full path to input directory of CDS', type=click.Path(), required=True)
@click.option('--codr', help='full path to output directory of CDS', type=click.Path(), required=True)
@click.option('--pidr', help='full path to input directory of Proteins', type=click.Path(), required=True)
@click.option('--podr', help='full path to output directory of Proteins', type=click.Path(), required=True)
@click.option('--minp', help='min fraction of ortholog presence', default=.75, show_default=True, type=click.FLOAT)
@click.option('--ncpu', help='number of CPUs to use', default=4, show_default=True, type=int)

def main_fun(fin,cidr,codr,pidr,podr,minp,ncpu):
    ## arguments
    print("Ortholog file: %s"%fin)
    print("Input CDS directory: %s"%cidr)
    print("Input Protein directory: %s"%pidr)
    print("Output CDS directory: %s"%codr)
    print("Output Protein directory: %s"%podr)    
    print("Min. Fraction of Genomes with the protein: %f"%minp)
    print("# CPUs: %i"%ncpu)

    if not os.path.exists(codr):
        os.mkdir(codr)

    if not os.path.exists(podr):
        os.mkdir(podr)

    odf = pandas.read_csv(fin,sep='\t',index_col=0,header=0)
    print("%s: Ortholog table had %i entries and %i fields"%(timer(start),odf.shape[0],odf.shape[1]))
    strains = odf.columns
    # extract single-copy orthologs
    odf = odf[ odf.apply( func = lambda x: not any(',' in i for i in x), axis=1)]
    nog,ns = odf.shape
    print("# Single-copy orthologs = %i"%nog)
    # extract orthologs present in minp fraction of genomes
    odf = odf[ odf.apply( func = lambda x: sum( i!='*' for i in x)/ns >= minp, axis=1)]
    print("# orthologs present in at least %f fraction of genomes = %i"%(minp,odf.shape[0]))
    
    with mp.Pool(ncpu) as mp_pool:
        allrecs = mp_pool.starmap( list_seqrecs, [(i,odf,cidr,pidr) for i in strains])

    nucrecs = [ i['cds'] for i in allrecs]
    protrecs = [ i['prots'] for i in allrecs]
    
    # transpose above matrix, so that all records for a gene comes in the same list
    tmat_n = [ [ j for j in i if j is not None] for i in zip(*nucrecs)]
    tmat_p = [ [ j for j in i if j is not None ] for i in zip(*protrecs)]
    nr = len(tmat_n)
    print("%s: # Orthologs sequences extracted from across genomes = %i"%(timer(start),nr))
    assert len(tmat_p) == nr
    assert nr == odf.shape[0]
    
    seengenes = []
    # for every ortholog group
    for x in range(nr):
        ids = [ i.id for i in tmat_n[x]]
        # gene name
        gnames = [ i.split('|')[1] for i in ids]
        gnames = [ i for i in gnames if i != '']

        try:    
            gene = Counter(gnames).most_common(1)[0][0]
        except IndexError:
            gene = rcloc.search(tmat_n[x][0].description).group(1)

        # add an index to gene name if its a copy
        orig_gname = gene

        if gene in seengenes:
            gene = gene + str( seengenes.count(gene))

        seengenes.append(orig_gname)        
        # output file paths
        fout_n = os.path.join( codr, gene+'.fasta')
        fout_p = os.path.join( podr, gene+'.fasta')
        # write sequence records to a file
        SeqIO.write( tmat_n[x], fout_n, 'fasta')
        SeqIO.write( tmat_p[x], fout_p, 'fasta')

    print("%s:Program Finished!\n"%timer(start))
        
if __name__ == '__main__':
    main_fun()
