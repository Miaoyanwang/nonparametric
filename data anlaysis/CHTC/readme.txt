transform.R: R script file that you want to run in the server
transform.sub: The file that you can specify the number of jobs, cpus, memory, and disk. you can specify the location where you want to save ``*.log, *.err, *.out'' files. You have to write down files you will use in the jobs such as R file, packages.tar.gz(packages you will use), transform.R(script file), *.Rdata
transform.sh: This file actually performs what your request. 1. untar R, 2. untar packages, 3. run the trnasform.R
