## tabulate results of Mcorr runs on different bacterial species
import os, sys, pandas
from glob import glob


indr = sys.argv[1] # input directory with sub-directories for each species
fout = sys.argv[2] # output table path
pat = sys.argv[3] #"*/mcorr/codon_lmfit_report.csv"

fls = glob( os.path.join( indr, pat))

def func(fin):
	with open(fin) as flob:
		contents = [ i.strip('\n') for i in flob]

	contents = [ i.split(',') for i in contents[-2:]]
	parm_val = { i[0]:float(i[1]) for i in zip(*contents)} 
	parm_val['rbym'] = parm_val['phi_p']/parm_val['theta_p']
	out = pandas.Series(parm_val)
	return out

mcr = [ func(f) for f in fls]
out = pandas.concat( mcr, axis = 1).transpose()
out.index = [ i.split('/')[-3] for i in fls]
out.to_csv(fout)
