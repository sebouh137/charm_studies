To make histograms at JLab:
1) edit submit.sh to have the locations for the output files be set to where you want to put the output files
2) Also make sure the input files exist in the cache.  If not, use the "jchache" utility:  "jcache get /mss/.../*.hipo" which will copy the input files to /cache/mss/.../*.hipo
3) Run the "sbatch submit.sh" command.  This will process all of the files from the input directory in parallel
4) Use "squeue -u $USER" to check the status of the event processing (should take about 12 minutes)
5) Use hadd to merge the output histogram files.  "hadd -f merged.root outputdir/*.root"
