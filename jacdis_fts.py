### Compute pairwise jaccard distances on feature tables of strains of a given species
## require
import os, sys, pandas, multiprocessing
from utils.mathsfun import jacdis

## arguments
indr = sys.argv[1] # full path to input directory of feature tables
fdis = sys.argv[2] # full path to output distance table
ncpu = int( sys.argv[3]) # number of CPUs for parallel processing

## file checks
if os.path.exists(fdis):
	raise FileExistsError

# list of feature table files
fls = os.listdir(indr)

## function to get the set of protein features for a given strain 
def get_ft( strain, dr = indr):
	fpath = os.path.join( dr, strain)
	
	try:
		df = pandas.read_csv( fpath, sep='\t')
	except:
		print("Failed to read the file!")
		return

	df = df[ (df.iloc[:,0]=='CDS') & (df.iloc[:,1]=='with_protein')]
	out = set(df['non-redundant_refseq'])
	return out

# dict of strains and their feature set
fts = {i:get_ft(i) for i in fls}
fts = { k:v for k,v in fts.items() if v is not None}
# list of strains with data
strains = list(fts.keys())
# number of strains
N = len(strains)
print("Feature set dict created.")

# initialize dataframe
odf = pandas.DataFrame( columns=strains, index=range(N))

# main function to parallelize over
def fun(x,y,names=strains,ft=fts):
	return (x,y,jacdis( ft[names[x]], ft[names[y]]))

# parallely run the above function
with multiprocessing.Pool(ncpu) as mp_pool:
	dls = mp_pool.starmap( fun, [(x,y) for x in range(N-1) for y in range(x+1,N)])

# fill dataframe with these values
for x,y,d in dls:
	odf.iat[y,x] = d

odf.to_csv( fdis, index=False)
print("Program finished!")