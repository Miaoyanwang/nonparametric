There are around 4000 RData files which represent cross validation result on each r and sparsity.
CVIQ* files are based on IQ brain dataset while CV* files on VSPLOT brain dataset.

CVIQ_sp_r_sparse.RData: a vector ``cvresult''  is saved, which represent log-likelihood values on each test datset.
combinations (r,sparse) = (14,25),(19,9) are ommited (the jobs failed to get output because of lack of memory)

CV_sp_r_sparse.RData:  a matrix``cvresult'' is saved, 
                       which the first row represent the log-likelihood values on each test datset 
                       while the second row saved the results on each training datset.
combinations (r,sparse) = (45,22),(50,4),(51,13) are ommited (the jobs failed to get output because of lack of memory)
                       
