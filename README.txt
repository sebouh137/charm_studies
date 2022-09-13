The histograms are created in src/MakeHistograms.C.  The slurm script submit.sh calls this program on all of the files in a given directory.  

To make histograms at JLab:
1) make sure the input files you want to use are accessible.  If they are on the MSS storage, you may need to cache them with jcache.
   For instance, if the files you want are .hipo files in the directory /mss/xxx/, running  "jcache get /mss/xxx/*.hipo" will cache
   these files in /cache/mss/xxx/.   This may take some time depending on the size of the files, and how busy the system is.
2) edit the INPUT_DIR and OUTPUT_DIR locations in submit.sh to specify the locations of the input and output files.  
3) Run the "sbatch submit.sh" command.  This will process all of the files from the input directory in parallel.  
4) Use "squeue -u $USER" to check the status of the event processing (should take about 12 minutes)
5) Use hadd to merge the output histogram files.  "hadd -f merged.root outputdir/*.root"
