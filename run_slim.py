import os, sys, multiprocessing, warnings, json, pandas
from contextlib import redirect_stdout
from subprocess import call

slim_script = sys.argv[1] # SLIM code
fprm = sys.argv[2] # table of parameters
fchr = sys.argv[3] # chromosome map
nrep = int(sys.argv[4]) # number of replicates
ncpu = int(sys.argv[5]) # number of CPUs

### Functions ####
def run_slim(slim_script,Ne,N_generations,Mu,Rho,tractlen,dsel,sellen,fchr,runID,logfile=None):
    
    slim_cmd = [ 'slim', '-t', '-m', '-l', '0', 
    "-d", "fchr='{}'".format(fchr),
    "-d", "runID='{}'".format(runID),
    "-d", "Ne={}".format(Ne),
    "-d", "N_generations={}".format(N_generations),
    "-d", "Mu={}".format(Mu),
    "-d", "Rho={}".format(Rho),
    "-d", "tractlen={}".format(tractlen),
    "-d", "dsel={}".format(dsel),
    "-d", "sellen={}".format(sellen),
    slim_script]

    if logfile is None:
        print(f"Simulation Parameters : \n"
        f"\tchromosomeMap = {fchr}\n"
        f"\tNe = {Ne}\n"
        f"\tt = {N_generations}\n"
        f"\tMu = {Mu}\n"
        f"\tRho = {Rho}\n"
        f"\ttractlen = {tractlen}\n"
        f"\ts = {dsel}\n"
        f"\tslen = {sellen}\n"
        )
        exit_code = call(slim_cmd)
    else:
        with open(logfile,'w') as lob:
            print(f"Simulation Parameters : \n"
            f"\tchromosomeMap = {fchr}\n"
            f"\tNe = {Ne}\n"
            f"\tt = {N_generations}\n"
            f"\tMu = {Mu}\n"
            f"\tRho = {Rho}\n"
            f"\ttractlen = {tractlen}\n"
            f"\ts = {dsel}\n"
            f"\tslen = {sellen}\n"
            ,file=lob)
            exit_code = call(slim_cmd,stdout=lob,stderr=lob)

    return exit_code

def func(params,slim_script=slim_script):
    ## simulation parameters
    runID = params["runID"]
    Ne = int(params["Ne"])
    N_generations = int(params["N_generations"])
    Mu = params["Mu"]
    Rho = params["Rho"]
    tractlen = int(params["tractlen"])
    dsel = params["dsel"]
    sellen = int(params["sellen"])
    flog = runID + ".log"
    ## run slim
    exit_code = run_slim(slim_script, Ne, N_generations, Mu, Rho, tractlen, dsel, sellen, fchr, runID, logfile=flog)
    return

## main
prms = pandas.read_csv(fprm)
# no. of parameter sets
np = prms.shape[0]
# add nrep copies of each line
prms = pandas.concat([ prms for _ in range(nrep)])
# add a column of name extensions
prms['ext'] = [ x+1 for x in range(nrep) for _ in range(np)]
# add the extension to runID
prms['runID'] = prms[['runID','ext']].apply( lambda x: '_'.join(x.values.astype(str)), axis=1)
# drop the extra columns
prms.drop('ext',axis=1,inplace=True)
# reset the index
prms.reset_index(inplace=True,drop=True)
parasets = prms.T.to_dict().values()

with multiprocessing.Pool(ncpu) as mp_pool:
	mp_pool.map( func, parasets)