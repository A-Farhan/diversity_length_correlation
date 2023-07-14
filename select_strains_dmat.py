### Select strains given a pairwise distance matrix

## require
import os, sys, pandas, numpy, random
#from utils.mathsfun import selbydist

## arguments
fin = sys.argv[1] # input table of pairwise distances
fout = sys.argv[2] # output list of selected elements
mind = float(sys.argv[3]) # minimum distance threshold

## Function 
def selbydist( df, cut, is_lower=True):
    """
    Select set of objects satisfying minimum distance threshold,
    given a dataframe of pairwise distances.

    Arguments:
    - df - input dataframe of distances
    - cut - minimum distance among all selected elements
    - is_lower - boolean specifying if the data is present in the lower triangle [True]
    """
    print("\tDistance threshold = %.2f"%mind)
    # strain names
    strains = df.columns
    df.index = strains
    # no. of objects
    R = len(strains)
    print("# Total strains = %i"%R)
    # convert to numpy 2-D array 
    dm = numpy.array(df.values)
    # transpose to shift data to upper triangle
    if is_lower:
        dm = numpy.transpose(dm)

    # initialize set of selected indices with the first one
    sel = {0}
    # empty set for the indices to be removed
    removed = set()
    # for every index, except the last
    for x in range(0,R-1):
        # if it is present in "removed", skip to next index
        if x in removed: continue
        # for every index following the present one
        for y in range(x+1,R):
            # if it is present in "removed", skip to next index
            if y in removed: continue
            # otherwise, get the distance between runs corresponding to these two indices
            d = dm[x,y]
            # if it is more than threshold, select the second index
            if d > cut:
                sel.add(y)
            # else, mask it from further consideration
            else:
                removed.add(y)
                # if the second index was previously selected, remove it
                if y in sel:
                    sel.remove(y)
    # get runs corresponding to selected indices
    out = [ i for x,i in enumerate(strains) if x in sel]
    return out

# load pairwise distance table
ddf = pandas.read_csv(fin)

sels = selbydist( ddf, mind)

with open( fout, 'w') as flob:
	for i in sels:
		flob.write( i + '\n')

print("# Selected strains = %i"%len(sels))
print("Program finished!")

