Rcpp::sourceCpp("src/vecchia_.cpp")

################################################################
# checking that obtained vecchia_ factor gives sensible samples #
################################################################

set.seed(100000)
n = 100000
# spatial locations 
locs = cbind(runif(n), runif(n))
NNarray = GpGp::find_ordered_nn(locs, 10)
NNarray[is.na(NNarray)] = 0



#isotropic, with smoothness = 1.5
set.seed(1)
tatato = array(0, dim = c(ncol(NNarray), n, 1 + ncol(NNarray)))
log_range = matrix(0, n, 1)
log_range[,1] = -8 + 3*locs[,1]
vecchia_(
  log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 8, 
  compute_derivative = T, smoothness = 1.5, result = tatato, start_idx = 1)
GeoNonStat::plotPointillistPainting(pch = ".", cex = 2,locs, GpGp::fast_Gp_sim_Linv (t(tatato[,,1]), NNarray))

# nominally anisotropic, with smoothness = 1.5
set.seed(1)
log_range = cbind(log_range[,1], log_range[,1], 0)
tatato = array(0, dim = c(ncol(NNarray), n, 1 + 3*ncol(NNarray)))
vecchia_(
  log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 8, compute_derivative = T, smoothness = 1.5, result = tatato, start_idx = 1)
GeoNonStat::plotPointillistPainting(pch = ".", cex = 2,locs, GpGp::fast_Gp_sim_Linv (t(tatato[,,1]), NNarray))

# actually anisotropic, with smoothness = 1.5
set.seed(1)
log_range[,2] = -10 + 3*locs[,2]
log_range[,3] = locs[,2] - 2*locs[,1]
vecchia_(
  log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 1, 
  compute_derivative = T, smoothness = 1.5, result = tatato, start_idx = 1)
GeoNonStat::plotPointillistPainting(pch = ".", cex = 2,locs, GpGp::fast_Gp_sim_Linv (t(tatato[,,1]), NNarray))

#isotropic, with smoothness = 0.5
set.seed(1)
tatato = array(0, dim = c(ncol(NNarray), n, 1 + ncol(NNarray)))
log_range = matrix(0, n, 1)
log_range[,1] = -8 + 3*locs[,1]
vecchia_(
  log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 8, 
  compute_derivative = T, smoothness = .5, result = tatato, start_idx = 1)
GeoNonStat::plotPointillistPainting(pch = ".", cex = 2,locs, GpGp::fast_Gp_sim_Linv (t(tatato[,,1]), NNarray))

# nominally anisotropic, with smoothness = .5
set.seed(1)
log_range = cbind(log_range[,1], log_range[,1], 0)
tatato = array(0, dim = c(ncol(NNarray), n, 1 + 3*ncol(NNarray)))
vecchia_(
  log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 8,
  compute_derivative = T, smoothness = .5, result = tatato, start_idx = 1)
GeoNonStat::plotPointillistPainting(pch = ".", cex = 2,locs, GpGp::fast_Gp_sim_Linv (t(tatato[,,1]), NNarray))

# actually anisotropic, with smoothness = 0.5
set.seed(1)
log_range[,2] = -10 + 3*locs[,2]
log_range[,3] = locs[,2] - 2*locs[,1]
vecchia_(
  log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 1, 
  compute_derivative = T, smoothness = .5, result = tatato, start_idx = 1)
GeoNonStat::plotPointillistPainting(pch = ".", cex = 2,locs, GpGp::fast_Gp_sim_Linv (t(tatato[,,1]), NNarray))


########################################################################
# Checking derivatives by plotting derivatives found with the formula  #
# against derivatives found with finite differences, ansisotropic case #
########################################################################

set.seed(1)
# spatial locations 
n = 1000
locs = cbind(runif(n), runif(n))
# range parameters 
log_range = matrix(0, n, 3)
log_range[,1] = -8 + 3*locs[,1]
log_range[,2] = -8 + 3*locs[,2]
log_range[,3] = -3*locs[,2]
# vecchia_ approx 
NNarray = GpGp::find_ordered_nn(locs, 10)
NNarray[is.na(NNarray)] = 0
smoothness = c(0.5, 1.5)
par(mfrow = c(2,2))

tatato = array(0, dim = c(ncol(NNarray), n, 1 + 3*ncol(NNarray)))
tatato_ = array(0, dim = c(ncol(NNarray), n, 1 + 3*ncol(NNarray)))

for(sm in smoothness){
  # computing vecchia_ 
  vecchia_(
    log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = T, smoothness = sm, tatato)
  # computing new vecchia_ approx with a new range parameter
  range_row_idx = 30
  range_col_idx = 2
  for(range_col_idx in seq(3)){
    # computing new vecchia_ approx and derivatrives
    log_range_ = log_range 
    log_range_[range_row_idx, range_col_idx] = log_range_[range_row_idx, range_col_idx] + .001
    vecchia_(
      log_range = t(log_range_), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = T, smoothness = sm, tatato_)
    
    # testing impact when impacted range is child in the DAG
    plot(
      (tatato_[,range_row_idx,1] - tatato[,range_row_idx,1])*1000, 
      tatato[,range_row_idx, 2 + 11 * (range_col_idx-1)], 
      xlab = "derivative using finite diff", ylab = "derivative computed using formula", 
      main = paste("test for the derivative of ", range_row_idx, 
                   "-th vecchia_ approx \n when its range is moved at column", range_col_idx, "\n smoothness = ", sm)
    )
    abline (a=0, b=1)
    
    
    # testing impact when impacted range is parent in the DAG
    
    idx_impacted_child = row(NNarray)[which(NNarray == range_row_idx)[6]]
    any(NNarray[idx_impacted_child,] == range_row_idx) # testing that the chosen index is present among parents of the impacted child
    idx_among_parents_of_impacted_child = col(NNarray)[which(NNarray == range_row_idx)[6]]-1
    (NNarray[idx_impacted_child,idx_among_parents_of_impacted_child + 1] == range_row_idx)# testing that the chosen index is at the right place among parents of the impacted child
    
    
    plot(
      (tatato_[,idx_impacted_child,1]-tatato[,idx_impacted_child,1])*1000,
      tatato[,idx_impacted_child, 2 + 11 * (range_col_idx-1) + idx_among_parents_of_impacted_child],
      xlab = "derivative using finite diff", ylab = "derivative computed using formula", 
      main = paste("test for the derivative of ", idx_impacted_child, 
                   "-th vecchia_ approx \n when the range of its", idx_among_parents_of_impacted_child, 
                   "-th parent is moved \n smoothness = ", sm )
    )
    abline(a = 0, b=1)
  }
}



########################################################################
# Checking derivatives by plotting derivatives found with the formula  #
# against derivatives found with finite differences, isotropic case    #
########################################################################

set.seed(1)
# spatial locations 
n = 1000
locs = cbind(runif(n), runif(n))
# range parameters 
log_range = matrix(0, n, 1)
log_range[,1] = -8 + 3*locs[,1]

# vecchia_ approx 
NNarray = GpGp::find_ordered_nn(locs, 10)
NNarray[is.na(NNarray)] = 0
smoothness = c(0.5, 1.5)

tatato = array(0, dim = c(ncol(NNarray), n,  1 + 1*ncol(NNarray)))
tatato_ = array(0, dim = c(ncol(NNarray), n, 1 + 1*ncol(NNarray)))

for(sm in smoothness){
  # computing vecchia_ 
  vecchia_(
    log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = T, smoothness = sm, tatato)
  # computing new vecchia_ approx with a new range parameter
  range_row_idx = 30
  range_col_idx = 1
  
  # computing new vecchia_ approx and derivatrives
  log_range_ = log_range 
  log_range_[range_row_idx, range_col_idx] = log_range_[range_row_idx, range_col_idx] + .001
  vecchia_(log_range = t(log_range_), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = T, smoothness = sm, tatato_)
  
  # testing impact when impacted range is child in the DAG
  plot(
    (tatato_[,range_row_idx,1] - tatato[,range_row_idx, 1])*1000, 
    tatato[,range_row_idx, 2 + 11 * (range_col_idx-1)], 
    xlab = "derivative using finite diff", ylab = "derivative computed using formula", 
    main = paste("test for the derivative of ", range_row_idx, 
                 "-th vecchia_ approx \n when its range is moved at column", range_col_idx, "\n smoothness = ", sm)
  )
  abline (a=0, b=1)
  
  
  # testing impact when impacted range is parent in the DAG
  
  idx_impacted_child = row(NNarray)[which(NNarray == range_row_idx)[6]]
  any(NNarray[idx_impacted_child,] == range_row_idx) # testing that the chosen index is present among parents of the impacted child
  idx_among_parents_of_impacted_child = col(NNarray)[which(NNarray == range_row_idx)[6]]-1
  (NNarray[idx_impacted_child,idx_among_parents_of_impacted_child + 1] == range_row_idx)# testing that the chosen index is at the right place among parents of the impacted child
  
  
  plot(
    (tatato_[,idx_impacted_child,1]-tatato[,idx_impacted_child,1])*1000,
    tatato[,idx_impacted_child, 2 + 11 * (range_col_idx-1) + idx_among_parents_of_impacted_child],
    xlab = "derivative using finite diff", ylab = "derivative computed using formula", 
    main = paste("test for the derivative of ", idx_impacted_child, 
                 "-th vecchia_ approx \n when the range of its", idx_among_parents_of_impacted_child, 
                 "-th parent is moved \n smoothness = ", sm )
  )
  abline(a = 0, b=1)
}



###############################
# testing derivative sandwich #
###############################

Rcpp::sourceCpp("src/vecchia_.cpp")

# spatial locations 
n = 10000
locs = cbind(runif(n), runif(n))
# range parameters 
log_range = matrix(0, n, 3)
log_range[,1] = -8 + 3*locs[,1]
log_range[,2] = -8 + 3*locs[,2]
log_range[,3] = -3*locs[,2]
# vecchia approx 
NNarray = GpGp::find_ordered_nn(locs, 10)
smoothness = c(0.5, 1.5)
# computing Vecchia 
tatato = array(0, dim = c(ncol(NNarray), n,  1 + 3*ncol(NNarray)))
tatato_ = array(0, dim = c(ncol(NNarray), n, 1 + 3*ncol(NNarray)))
vecchia_(
  log_range = t(log_range), locs = t(locs), 
  NNarray = t(NNarray), num_threads = 8, compute_derivative = T, 
  smoothness = 1.5, tatato)
# computing new Vecchia approx with a new range parameter
range_row_idx = 30
range_col_idx = 1


for(range_col_idx in seq(3)){
  print(paste("range_col_idx = ", range_col_idx))
  # computing new Vecchia approx and derivatrives
  log_range_ = log_range 
  log_range_[range_row_idx, range_col_idx] = log_range_[range_row_idx, range_col_idx] + .00001
  vecchia_(
    log_range = t(log_range_), locs = t(locs), NNarray = t(NNarray), num_threads = 8, compute_derivative = T, smoothness = 1.5, tatato_)
  
  
  v_left = rnorm(n)
  v_right = rnorm(n)
  
  t1 = Sys.time()
  tatata = derivativeSandwiches_(vecchia = tatato, left_vector = v_left, right_vector = v_right, NNarray = t(NNarray), num_threads = 1, sauce_determinant_chef = F)
  print(paste("computing derivative sandwich took", Sys.time()-t1, "s"))
  
  print("analytical sandwich")
  print(((
    v_left %*% Matrix::sparseMatrix(i = row(NNarray)[!is.na(NNarray)], j= NNarray[!is.na(NNarray)], x= t(tatato_[,,1])[!is.na(NNarray)]) %*% v_right -
      v_left %*% Matrix::sparseMatrix(i = row(NNarray)[!is.na(NNarray)], j= NNarray[!is.na(NNarray)], x= t(tatato [,,1])[!is.na(NNarray)]) %*% v_right
  )*100000)[1,1])
  print("finite diff sandwich")
  print(tatata[range_col_idx, range_row_idx])
  
  
  t1 = Sys.time()
  tatata = derivativeSandwiches_(vecchia = tatato, left_vector = v_left, right_vector = v_right, NNarray = t(NNarray), num_threads = 4, sauce_determinant_chef = T)
  print(paste("computing derivative sandwich with determinant sauce  took", Sys.time()-t1, "s"))
  print("analytical sandwich")
  print(((
    v_left %*% Matrix::sparseMatrix(i = row(NNarray)[!is.na(NNarray)], j= NNarray[!is.na(NNarray)], x= t(tatato_[,,1])[!is.na(NNarray)]) %*% v_right -
      v_left %*% Matrix::sparseMatrix(i = row(NNarray)[!is.na(NNarray)], j= NNarray[!is.na(NNarray)], x= t(tatato [,,1])[!is.na(NNarray)]) %*% v_right 
    + sum(log(tatato_[1,,1])) - sum(log(tatato[1,,1]))
  )*100000)[1,1])
  print("finite diff sandwich with determinant sauce")
  print(tatata[range_col_idx, range_row_idx])
  print("=========================")
  
  print("determinant only, analytical")
  ttt = derivativeSandwiches_(vecchia = tatato, left_vector = 0*v_left, right_vector = v_right, NNarray = t(NNarray), num_threads = 4, sauce_determinant_chef = T)
  print(ttt[range_col_idx, range_row_idx])
  print("determinant only, finite diff")
  print(sum(+log(tatato_[1,,1]) - log(tatato[1,,1]))*100000)
  print("=========================")
}


########################
# speed test
########################





Rcpp::sourceCpp("src/vecchia_.cpp")

# spatial locations 
n = 200000
locs = cbind(runif(n), runif(n))
# range parameters 
log_range = matrix(0, n, 3)
log_range[,1] = -8 + 3*locs[,1]
log_range[,2] = -8 + 3*locs[,2]
log_range[,3] = -3*locs[,2]
# vecchia approx 
NNarray = GpGp::find_ordered_nn(locs, 10)
smoothness = c(0.5, 1.5)
# computing Vecchia 

tatato = array(0, dim = c(ncol(NNarray), n,  1 + 3*ncol(NNarray)))



timediff1 = mean(sapply(seq(10), function(x){
  t1 = Sys.time()
vecchia_(
  log_range = t(log_range), locs = t(locs), 
  NNarray = t(NNarray), num_threads = 1, compute_derivative = T, 
  smoothness = 1.5, tatato)
  Sys.time()-t1
}))

timediff5 = mean(sapply(seq(10), function(x){
  t1 = Sys.time()
vecchia_(
  log_range = t(log_range), locs = t(locs), 
  NNarray = t(NNarray), num_threads = 5, compute_derivative = T, 
  smoothness = 1.5, tatato)
  Sys.time()-t1
}))

timediff8 = mean(sapply(seq(10), function(x){
  t1 = Sys.time()
vecchia_(
  log_range = t(log_range), locs = t(locs), 
  NNarray = t(NNarray), num_threads = 8, compute_derivative = T, 
  smoothness = 1.5, tatato)
  Sys.time()-t1
}))

timediff10 = mean(sapply(seq(50), function(x){
  t1 = Sys.time()
vecchia_(
  log_range = t(log_range), locs = t(locs), 
  NNarray = t(NNarray), num_threads = 10, compute_derivative = T, 
  smoothness = 1.5, tatato)
  Sys.time()-t1
}))

timediff20 = mean(sapply(seq(10), function(x){
  t1 = Sys.time()
vecchia_(
  log_range = t(log_range), locs = t(locs), 
  NNarray = t(NNarray), num_threads = 20, compute_derivative = T, 
  smoothness = 1.5, tatato)
  Sys.time()-t1
}))

plot(c(1,5,8,10,20), c(timediff1, timediff5, timediff8, timediff10, timediff20), ylim = c(0, timediff1))










########################################################################
# Checking derivatives with high range
########################################################################

set.seed(1)
# spatial locations 
n = 30000
locs = cbind(runif(n), runif(n))
# range parameters 
log_range = matrix(0, n, 1)
log_range[,1] = 2 + 3*locs[,1]

# vecchia_ approx 
NNarray = GpGp::find_ordered_nn(locs, 10)
NNarray[is.na(NNarray)] = 0
smoothness = c(0.5, 1.5)

tatato = array(0, dim = c(ncol(NNarray), n,  1 + 1*ncol(NNarray)))
tatato_ = array(0, dim = c(ncol(NNarray), n, 1 + 1*ncol(NNarray)))

for(sm in smoothness){
  # computing vecchia_ 
  vecchia_(
    log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = T, smoothness = sm, tatato)
  GeoNonStat::plotPointillistPainting(pch = 16, cex = .8,locs, GpGp::fast_Gp_sim_Linv (t(tatato[,,1]), NNarray))
  
  
  # computing new vecchia_ approx with a new range parameter
  range_row_idx = 15000
  range_col_idx = 1
  
  # computing new vecchia_ approx and derivatrives
  log_range_ = log_range 
  log_range_[range_row_idx, range_col_idx] = log_range_[range_row_idx, range_col_idx] + .00001
  vecchia_(log_range = t(log_range_), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = T, smoothness = sm, tatato_)
  
  # testing impact when impacted range is child in the DAG
  plot(
    (tatato_[,range_row_idx,1] - tatato[,range_row_idx, 1])*100000, 
    tatato[,range_row_idx, 2 + 11 * (range_col_idx-1)], 
    xlab = "derivative using finite diff", ylab = "derivative computed using formula", 
    main = paste("test for the derivative of ", range_row_idx, 
                 "-th vecchia_ approx \n when its range is moved at column", range_col_idx, "\n smoothness = ", sm)
  )
  abline (a=0, b=1)
  
  
  # testing impact when impacted range is parent in the DAG
  
  idx_impacted_child = row(NNarray)[which(NNarray == range_row_idx)[6]]
  any(NNarray[idx_impacted_child,] == range_row_idx) # testing that the chosen index is present among parents of the impacted child
  idx_among_parents_of_impacted_child = col(NNarray)[which(NNarray == range_row_idx)[6]]-1
  (NNarray[idx_impacted_child,idx_among_parents_of_impacted_child + 1] == range_row_idx)# testing that the chosen index is at the right place among parents of the impacted child
  
  
  plot(
    (tatato_[,idx_impacted_child,1]-tatato[,idx_impacted_child,1])*100000,
    tatato[,idx_impacted_child, 2 + 11 * (range_col_idx-1) + idx_among_parents_of_impacted_child],
    xlab = "derivative using finite diff", ylab = "derivative computed using formula", 
    main = paste("test for the derivative of ", idx_impacted_child, 
                 "-th vecchia_ approx \n when the range of its", idx_among_parents_of_impacted_child, 
                 "-th parent is moved \n smoothness = ", sm )
  )
  abline(a = 0, b=1)
}









