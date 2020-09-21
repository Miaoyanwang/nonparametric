# we need 3 files to request jobs on the server

1. transform.R: R script file that you want to run in the server

2. transform.sub: The file that you can specify the number of jobs, cpus, memory, and disk. 
                  You can specify the location where you want to save ``*.log, *.err, *.out'' files. 
                  You have to write down files needed for the jobs: 
                   ex) R361.tar.gz(Rfile), packages.tar.gz(R packages), transform.R(script file), *.Rdata(Dataset)
   
3. transform.sh: This file actually performs what your request. 
                 1) untar R. 2) untar packages. 3) run the trnasform.R


For bbnet68_spatial_orientation.RData, I used my_rscript.R, run_R.sub, run_R.sh

For brain_binary_IQ.RData, I used my_rscriptIQ.R, run_RIQ.sub, run_RIQ.sh
