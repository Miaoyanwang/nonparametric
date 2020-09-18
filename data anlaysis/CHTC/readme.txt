# we need 3 files to request jobs on the server

1. transform.R: R script file that you want to run in the server

2. transform.sub: The file that you can specify the number of jobs, cpus, memory, and disk. 
   You can specify the location where you want to save ``*.log, *.err, *.out'' files. 
   You have to write down files needed for the jobs: R361.tar.gz(Rfile), packages.tar.gz(packages you will use), transform.R(script file), *.Rdata
   
3. transform.sh: This file actually performs what your request. 
                 1) untar R. 2) untar packages. 3) run the trnasform.R
