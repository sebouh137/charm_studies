mkdir /work/clas12/spaul/tmp/

for dir in fall2018/torus+1/pass2/dst/recon/ fall2018/torus-1/pass2/main/dst/recon/ spring2019/torus-1/pass2/dst/recon/ ; do
    for run in `ls /cache/clas12/rg-a/production/recon/$dir/` ; do
	#skip the README.json files
	if [[ $run == *"json"* ]]; then
            continue  # Continue to the next iteration of the loop
	fi
	input_dir=/cache/clas12/rg-a/production/recon/${dir}/${run}/
	#echo sbatch ./filter_and_merge.sh $input_dir /work/clas12/spaul/tmp/${run} /work/clas12/spaul/lcp_skim/merged_${run}.hipo
	sbatch ./filter_and_merge.sh $input_dir /work/clas12/spaul/tmp/${run} /work/clas12/spaul/lcp_skim/merged_${run}.hipo
    done
done

