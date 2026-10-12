#!/bin/bash
#!/bin/bash
#
#SBATCH --job-name=clas12_mvtx
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
####  SBATCH --partition=short
#SBATCH --time=01:00:00

##set this to the number of files in INPUT_DIR
#SBATCH --array=0-99
#SBATCH --mem-per-cpu=2500
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sebouh.paul@gmail.com
i=$SLURM_ARRAY_TASK_ID

topo=D0bar

dir=/work/clas12/spaul/charm_sim/lund_${topo}

lundfile=${dir}/lund_${topo}_${i}.lund
simfile=${dir}/${topo}_sim_${i}.hipo
recfile=${dir}/${topo}_rec_${i}.hipo
tuplesfile=${dir}/${topo}_tuples_${i}.root

GCARD=clas12_mvtx_and_foils.gcard
RECON_YAML=/group/clas12/packages/coatjava/10.1.0/config/mc-cv.yaml

vtx=$(( -50 + 5*($i % 5) ))

if (( i % 5 == 2 )); then
    half_thickness=0.025
else
    half_thickness=0.05
fi

rm -rf $recfile

gemc $GCARD -USE_GUI=0 -INPUT_GEN_FILE="LUND, "${lundfile} -USE_GUI=0 -OUTPUT="hipo, "${simfile} -N=1000 -SHIFT_LUND_VERTEX="(0, 0, "${vtx}")mm" -BEAM_SPOT="0*mm, 0*mm, 0.2*mm, 0.2*mm, 0*deg" -INTEGRATEDRAW='flux'
recon-util -i $simfile -o $recfile -y $RECON_YAML
clas12root -l -b -q ~/charm_studies/src/MakeHistogramsMC_noPID.C --in=$recfile --out=$tuplesfile --isMC
