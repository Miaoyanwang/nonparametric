#!/bin/bash

# untar your R installation. Make sure you are using the right version!
tar -xzf R361.tar.gz
# (optional) if you have a set of packages (created in Part 1), untar them also
tar -xzf packages.tar.gz
#tar -xzf data.tar.gz
#mv $(pwd)/data/*.csv $(pwd)/

# make sure the script will use your R installation, 
# and the working directory as its home location
export PATH=$PWD/R/bin:$PATH
export RHOME=$PWD/R
export R_LIBS=$PWD/packages

# run your script
Rscript my_rscriptIQ.R $1
