import uproot as ur
import pandas as pd
import numpy as np
import glob
#help(pd.read_csv)
good_runs=[]
for f in glob.glob("good_run_lists/*"):
    df_run_period=pd.read_csv(f, names=("run","charge"))
    print(f, len(df_run_period))                              
    good_runs=np.concatenate([good_runs,df_run_period.run])

good_runs=sorted(good_runs)

def read_root(filename,treename=None,N=None):
    if type(filename) != str:
        return pd.concat([read_root(f,treename,N) for f in filename])
        
    with ur.open(filename) as f:
        if treename is None:
            if len(f.keys()) == 1:
                treename = f.keys()[0]
            elif len(f.keys()) == 0:
                raise Exception("error: no tree names in file " + filename)
            else:
                raise Exception("error: treename must be specified; more than one tree in " + filename)
            
        df = f[treename].arrays(library="pd")
    return df

#input_file="hadd.root"
#topo=0
#output_file="D0bar_topo.pkl"
input_file="lcp.root"
topo=2
output_file="lcp_topo.pkl"

tmp=read_root(input_file, "events")
tmp_cut=tmp.query(f"topo=={topo}")
tmp_goodruns=tmp_cut[tmp_cut["run_num"].isin(good_runs)]
tmp_goodruns.reset_index(drop=True).to_pickle(output_file)
