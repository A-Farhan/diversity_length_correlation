### generate protein MSAs and convert to codon MSAs

## require
import os, click, multiprocessing, subprocess
from time import time
from Bio.Align.Applications import ClustalOmegaCommandline as clustalo
from Bio import SeqIO
from utils.basic import timer

# start clock
start = time()

### Functions ###
def clustal(infile,outdr,prog='clustalo',outfmt='fa'):    
    """Run Clustal Omega to generate protein MSAs."""
    # basename of input file path
    fname = os.path.basename(infile)
    # extract gene name
    gene = fname.replace('.fasta', '')
    # generate output file path
    outfile = os.path.join( outdr, fname)
    
    # if output already presents, then skip the computation
    if os.path.exists(outfile):
        print("output- %s - already present!"%outfile)
        return

    # generate the command
    cmd = clustalo(prog,infile=infile,outfile=outfile,outfmt=outfmt)
    # run the command
    cmd() 

def prot2cds_msa(infile,pdr,cidr,codr,prog='pal2nal'):
    """
    Convert protein MSA to CODON MSA.

    Arguments:
    - infile    - input protein MSA file name
    - pdr       - directory path where these MSAs are stored
    - cidr      - directory path for CDS multi-FASTA files
    - codr      - output directory path for codon based MSAs
    - prog      - program to convert amino-acid alignments to codon alignments
    """

    # output file path
    f_out = os.path.join( codr, infile)
    
    # terminate, if output already exists
    if os.path.exists(f_out):
        print("Output- %s - already present!"%f_out)
        return

    # protein MSA file path
    fp = os.path.join( pdr, infile)
    precs = SeqIO.parse( fp, 'fasta')
    # cds multi-FASTA file path
    fn = os.path.join( cidr, infile)
    nrecs = SeqIO.parse( fn, 'fasta')

    # generate commands to convert amino acid alignments to codon alignments
    cmd = [ prog, fp, fn, '-output', 'fasta', '-codontable','11']
    with open(f_out,'w') as flob:
        s = subprocess.run(cmd, stdout=flob, stderr=subprocess.DEVNULL)    

#######################################################   

@click.command()
@click.option('--cidr', help='full path to input directory of CDS', type=click.Path(), required=True)
@click.option('--codr', help='full path to output directory of codon MSAs', type=click.Path(), required=True)
@click.option('--pidr', help='full path to input directory of Proteins', type=click.Path(), required=True)
@click.option('--podr', help='full path to output directory of protein MSAs', type=click.Path(), required=True)
@click.option('--aprog', help='program to generate protein MSAs, full path is required if it is not available\
        system-wide', default='clustalo', show_default=True)
@click.option('--cprog', help='''program to generate codon-based MSAs using protein MSAs and
        CDS multi-FASTA files, full path is required if it is not available system-wide''',
        default='pal2nal',show_default=True)
@click.option('--ncpu', help='number of CPUs to use', default=4, show_default=True, type=int)

def main_fun(cidr,pidr,codr,podr,aprog,cprog,ncpu):
    """This program generates multiple sequence alignments for both amino acids and codons."""

    if not os.path.exists(podr):
        click.echo("Creating protein MSA directory, %s"%podr)
        os.mkdir(podr)
    
    if not os.path.exists(codr):
        click.echo("Creating codon MSA directory, %s"%codr)
        os.mkdir(codr)

    # list of cds FASTA files
    fls_n = [ os.path.join( cidr, i) for i in os.listdir(cidr)]
    # list of protein FASTA files
    fls_p = [ os.path.join( pidr, i) for i in os.listdir(pidr)]

    ## Generate MSA ##
    # run the clustalo function parallely
    print('%s: Building MSAs of proteins\n'%timer(start))
    with multiprocessing.Pool(ncpu) as mp_pool:
        mp_pool.starmap( clustal, [(i,podr,aprog) for i in fls_p])
    
    
    ## Generate codon-alignments with protein MSA using PAL2NAL ###
    # list of files of protein MSA
    fls = os.listdir(podr)
    print("%s: Running PAL2NAL to generate codon alignments\n"%timer(start))
    
    with multiprocessing.Pool(ncpu) as mp_pool:
        mp_pool.starmap( prot2cds_msa, [(i,podr,cidr,codr,cprog) for i in fls]) 

    print("%s: Program finished!\n"%timer(start))

if __name__ == '__main__':
    main_fun()
