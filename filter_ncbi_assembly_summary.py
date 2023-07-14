### Filter NCBI assembly summary for Bacteria, and select species with a minimum assembly count

## require
import sys, pandas

## arguments

fin = sys.argv[1] # input NCBI assembly summary table
fout = sys.argv[2] # filtered output table
mina = int(sys.argv[3]) # min no. of assemblies for each species

df = pandas.read_csv( fin, sep='\t',skiprows=1) 
print("Original dataframe has %i rows"%df.shape[0])

# remove those excluded from refseq
df = df[df['excluded_from_refseq'].isna()]
print("# rows excluding those without refseq = %i"%df.shape[0])
# remove misclassified ones
df = df[-df['organism_name'].str.contains("[",regex=False)]
print("# rows excluding misclassified ones = %i"%df.shape[0])
# remove uncultured ones
df = df[-df['organism_name'].str.contains("uncultured",regex=False)]
print("# rows excluding uncultured ones = %i"%df.shape[0])
# remove those without species names
df = df[-df['organism_name'].str.contains("sp.",regex=False)]
print("# rows excluding those w/o species names = %i"%df.shape[0])
# remove those with candidatus in the name
df = df[-df['organism_name'].str.contains("Candidatus",regex=False)]
print("# rows excluding Candidatus = %i"%df.shape[0])
# remove those with bacterium
df = df[~df['organism_name'].str.contains('\s{1}bacterium$')]
print("# rows excluding bacterium = %i"%df.shape[0])
# reduce organism names to first two components
df['organism_name'] = df['organism_name'].str.split(' ').apply( lambda x: ' '.join(x[:2]))
# consider only chromosomal-level or higher assemblies
df = df[df['assembly_level'].isin(['Chromosome','Complete Genome'])]
print("# rows with assembly level chromosomal or higher = %i"%df.shape[0])
# get assembly count per species
spc = df['organism_name'].value_counts()
spc = spc[spc >= mina]
df = df[df['organism_name'].isin(spc.index)]
print("# rows for species with at least %i assemblies = %i"%(mina,df.shape[0]))
df.to_csv(fout,sep='\t',index=False)
