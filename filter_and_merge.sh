#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH --time=02:00:00

#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sebouh.paul@gmail.com
input_dir=$1
tmp_dir=$2
output_file=$3
mkdir $tmp_dir
for f in `ls $input_dir/`; do
    ./LcpSkim $input_dir/$f $tmp_dir/$f
done
rm $output_file
hipo-utils -merge -o $output_file $tmp_dir/*.hipo
#rm -rf $tmp_dir
