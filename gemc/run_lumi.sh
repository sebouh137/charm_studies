#!/bin/bash
#!/bin/bash

#SBATCH --job-name=mvtx_lumi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
####  SBATCH --partition=short
#SBATCH --time=02:00:00

##set this to the number of files in INPUT_DIR
#SBATCH --array=0-99
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sebouh.paul@gmail.com
i=$SLURM_ARRAY_TASK_ID

outdir=/work/clas12/spaul/charm_sim/lumi/
mkdir -p $outdir


outfile=${outdir}/test_bkgd_maxlumi_${i}.hipo
recfile=${outdir}/test_bkgd_maxlumi_rec_${i}.hipo
tuplefile=${outdir}/flux_tuples_bkgd_maxlumi_${i}.root
GCARD=mvtx_and_foils.gcard
#gemc $GCARD -USE_GUI=0 -OUTPUT="hipo, "${outfile} -N=1 -INTEGRATEDRAW='flux' -LUMI_EVENT="14000000, 5000*ns, 4.008*ns" -LUMI_P="e-, 10.6*GeV, 0*deg, 0*deg" -LUMI_V="(0.0, 0.0, -10.0)cm" -LUMI_SPREAD_V="(0.02, 0.02)cm"

gemc $GCARD -USE_GUI=0 -OUTPUT="hipo, "${outfile} -N=50 -INTEGRATEDRAW='flux' -LUMI_EVENT="88000, 100*ns, 4.008*ns" -LUMI_P="e-, 10.6*GeV, 0*deg, 0*deg" -LUMI_V="(0.0, 0.0, -10.0)cm" -LUMI_SPREAD_V="(0.02, 0.02)cm"

rm $recfile
recon-util -i ${outfile} -o ${recfile} -y /group/clas12/packages/coatjava/10.1.0/config/mc-cv.yaml
clas12root -l -b -q ../src/FluxTuples.C --in=${recfile} --out=${tuplefile}



