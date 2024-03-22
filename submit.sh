#!/bin/bash
#
#SBATCH --job-name=****
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
####  SBATCH --partition=short
#SBATCH --time=02:00:00

##set this to the number of files in INPUT_DIR
#SBATCH --array=0-537   
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sebouh.paul@gmail.com
srun hostname

echo $SLURM_ARRAY_TASK_ID

TRAIN_NAME=ElecFTKaon
INPUT_DIR1=/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/${TRAIN_NAME}/
INPUT_DIR2=/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/${TRAIN_NAME}/
INPUT_DIR3=/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/supplemental/train/${TRAIN_NAME}/
INPUT_DIR4=/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/${TRAIN_NAME}/
OUTPUT_DIR=/work/clas12/spaul/charm_hists/

#include both inbending and outbending data here
FILE=`ls ${INPUT_DIR1}/*.hipo ${INPUT_DIR2}/*.hipo ${INPUT_DIR3}/*.hipo ${INPUT_DIR4}/*.hipo | awk '{if (NR=='$SLURM_ARRAY_TASK_ID'+1) print $0;}' `


echo ${FILE}
OUTFILE=${OUTPUT_DIR}/$(basename ${FILE} .hipo).root


srun clas12root -l -b -q /home/spaul/charm_studies/src/MakeHistograms.C --in=$FILE --out=$OUTFILE --FTonly
