#!/bin/bash
#
#SBATCH --job-name=****
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
####  SBATCH --partition=short
#SBATCH --time=00:20:00

##set this to the number of files in INPUT_DIR
#SBATCH --array=0-185   
#SBATCH --mem-per-cpu=16000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sebouh.paul@gmail.com
srun hostname

echo $SLURM_ARRAY_TASK_ID

INPUT_DIR=/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass1/v1/dst/train/ElecFTKaon/
OUTPUT_DIR=/work/clas12/spaul/charm_hists/

FILE=`ls ${INPUT_DIR}/*.hipo | awk '{if (NR=='$SLURM_ARRAY_TASK_ID') print $0;}' `


echo ${FILE}
OUTFILE=${OUTPUT_DIR}/$(basename ${FILE} .hipo).root


srun clas12root -l -b -q /home/spaul/charm_studies/src/MakeHistograms.C+ --in=$FILE --out=$OUTFILE --FTonly


