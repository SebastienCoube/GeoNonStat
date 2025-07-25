Rcpp::sourceCpp("src/vecchia.cpp")

################################################################
# checking that obtained Vecchia factor gives sensible samples #
################################################################

set.seed(100000)
n = 100000
# spatial locations 
locs = cbind(runif(n), runif(n))
log_range = matrix(0, n, 1)
log_range[,1] = -10 + 3*locs[,1]
NNarray = GpGp::find_ordered_nn(locs, 10)
NNarray[is.na(NNarray)] = 0

#isotropic, with smoothness = 1.5
set.seed(1)
tatato = vecchia(
  log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = F, smoothness=  1.5)
GeoNonStat::plot_pointillist_painting(locs, GpGp::fast_Gp_sim_Linv (t(tatato[,1,]), NNarray))

# nominally anisotropic, with smoothness = 1.5
set.seed(1)
log_range = cbind(log_range, log_range, 0)
tatato = vecchia(
  log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = F, smoothness = 1.5)
GeoNonStat::plot_pointillist_painting(locs, GpGp::fast_Gp_sim_Linv (t(tatato[,1,]), NNarray))

# actually anisotropic, with smoothness = 1.5
set.seed(1)
log_range[,2] = -10 + 3*locs[,2]
log_range[,3] = locs[,2] - 2*locs[,1]
tatato = vecchia(
  log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = F, smoothness = 1.5)
GeoNonStat::plot_pointillist_painting(locs, GpGp::fast_Gp_sim_Linv (t(tatato[,1,]), NNarray))

#isotropic, with smoothness = 0.5
set.seed(1)
log_range = matrix(0, n, 1)
log_range[,1] = -10 + 3*locs[,1]
tatato = vecchia(
  log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = F, smoothness=  0.5)
GeoNonStat::plot_pointillist_painting(locs, GpGp::fast_Gp_sim_Linv (t(tatato[,1,]), NNarray))

# nominally anisotropic, with smoothness = 0.5
set.seed(1)
log_range = cbind(log_range, log_range, 0)
tatato = vecchia(
  log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = F, smoothness = 0.5)
GeoNonStat::plot_pointillist_painting(locs, GpGp::fast_Gp_sim_Linv (t(tatato[,1,]), NNarray))

# actually anisotropic, with smoothness = 0.5
set.seed(1)
log_range[,2] = -10 + 3*locs[,2]
log_range[,3] = locs[,2] - 2*locs[,1]
tatato = vecchia(
  log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = F, smoothness = 0.5)
GeoNonStat::plot_pointillist_painting(locs, GpGp::fast_Gp_sim_Linv (t(tatato[,1,]), NNarray))


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
# vecchia approx 
NNarray = GpGp::find_ordered_nn(locs, 10)
NNarray[is.na(NNarray)] = 0
smoothness = c(0.5, 1.5)
for(sm in smoothness){
  # computing Vecchia 
  tatato = vecchia(
    log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = T, smoothness = sm)
  # computing new Vecchia approx with a new range parameter
  range_row_idx = 30
  range_col_idx = 2
  for(range_col_idx in seq(3)){
    # computing new Vecchia approx and derivatrives
    log_range_ = log_range 
    log_range_[range_row_idx, range_col_idx] = log_range_[range_row_idx, range_col_idx] + .001
    tatato_ = vecchia(
      log_range = t(log_range_), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = T, smoothness = sm)
    
    # testing impact when impacted range is child in the DAG
    plot(
      (tatato_[,1,range_row_idx] - tatato[,1,range_row_idx])*1000, 
      tatato[,2 + 11 * (range_col_idx-1),range_row_idx], 
      xlab = "derivative using finite diff", ylab = "derivative computed using formula", 
      main = paste("test for the derivative of ", range_row_idx, 
                   "-th Vecchia approx \n when its range is moved at column", range_col_idx, "\n smoothness = ", sm)
    )
    abline (a=0, b=1)
    
    
    # testing impact when impacted range is parent in the DAG
    
    idx_impacted_child = row(NNarray)[which(NNarray == range_row_idx)[6]]
    any(NNarray[idx_impacted_child,] == range_row_idx) # testing that the chosen index is present among parents of the impacted child
    idx_among_parents_of_impacted_child = col(NNarray)[which(NNarray == range_row_idx)[6]]-1
    (NNarray[idx_impacted_child,idx_among_parents_of_impacted_child + 1] == range_row_idx)# testing that the chosen index is at the right place among parents of the impacted child
    
    
    plot(
      (tatato_[,1,idx_impacted_child]-tatato[,1,idx_impacted_child])*1000,
      tatato[,2 + 11 * (range_col_idx-1) + idx_among_parents_of_impacted_child,idx_impacted_child],
      xlab = "derivative using finite diff", ylab = "derivative computed using formula", 
      main = paste("test for the derivative of ", idx_impacted_child, 
                   "-th Vecchia approx \n when the range of its", idx_among_parents_of_impacted_child, 
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

# vecchia approx 
NNarray = GpGp::find_ordered_nn(locs, 10)
NNarray[is.na(NNarray)] = 0
smoothness = c(0.5, 1.5)
for(sm in smoothness){
  # computing Vecchia 
  tatato = vecchia(
    log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = T, smoothness = sm)
  # computing new Vecchia approx with a new range parameter
  range_row_idx = 30
  range_col_idx = 1
  
  # computing new Vecchia approx and derivatrives
  log_range_ = log_range 
  log_range_[range_row_idx, range_col_idx] = log_range_[range_row_idx, range_col_idx] + .001
  tatato_ = vecchia(
    log_range = t(log_range_), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = T, smoothness = sm)
  
  # testing impact when impacted range is child in the DAG
  plot(
    (tatato_[,1,range_row_idx] - tatato[,1,range_row_idx])*1000, 
    tatato[,2 + 11 * (range_col_idx-1),range_row_idx], 
    xlab = "derivative using finite diff", ylab = "derivative computed using formula", 
    main = paste("test for the derivative of ", range_row_idx, 
                 "-th Vecchia approx \n when its range is moved at column", range_col_idx, "\n smoothness = ", sm)
  )
  abline (a=0, b=1)
  
  
  # testing impact when impacted range is parent in the DAG
  
  idx_impacted_child = row(NNarray)[which(NNarray == range_row_idx)[6]]
  any(NNarray[idx_impacted_child,] == range_row_idx) # testing that the chosen index is present among parents of the impacted child
  idx_among_parents_of_impacted_child = col(NNarray)[which(NNarray == range_row_idx)[6]]-1
  (NNarray[idx_impacted_child,idx_among_parents_of_impacted_child + 1] == range_row_idx)# testing that the chosen index is at the right place among parents of the impacted child
  
  
  plot(
    (tatato_[,1,idx_impacted_child]-tatato[,1,idx_impacted_child])*1000,
    tatato[,2 + 11 * (range_col_idx-1) + idx_among_parents_of_impacted_child,idx_impacted_child],
    xlab = "derivative using finite diff", ylab = "derivative computed using formula", 
    main = paste("test for the derivative of ", idx_impacted_child, 
                 "-th Vecchia approx \n when the range of its", idx_among_parents_of_impacted_child, 
                 "-th parent is moved \n smoothness = ", sm )
  )
  abline(a = 0, b=1)
  
}



# testing derivative sandwich

set.seed(1)
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
tatato = vecchia(
  log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 6, compute_derivative = T, smoothness = 1.5)
# computing new Vecchia approx with a new range parameter
range_row_idx = 30
range_col_idx = 1
for(range_col_idx in seq(3)){
  print(paste("range_col_idx = ", range_col_idx))
  # computing new Vecchia approx and derivatrives
  log_range_ = log_range 
  log_range_[range_row_idx, range_col_idx] = log_range_[range_row_idx, range_col_idx] + .00001
  tatato_ = vecchia(
    log_range = t(log_range_), locs = t(locs), NNarray = t(NNarray), num_threads = 6, compute_derivative = T, smoothness = 1.5)
  
  
  v_left = rnorm(n)
  v_right = rnorm(n)
  
  t1 = Sys.time()
  tatata = derivative_sandwiches(vecchia = tatato, left_vector = v_left, right_vector = v_right, NNarray = t(NNarray), num_threads = 4, sauce_determinant_chef = F)
  print(paste("computing derivative sandwich took", Sys.time()-t1, "s"))

  print("analytical sandwich")
  print(((
    v_left %*% Matrix::sparseMatrix(i = row(NNarray)[!is.na(NNarray)], j= NNarray[!is.na(NNarray)], x= t(tatato_[,1,])[!is.na(NNarray)]) %*% v_right -
      v_left %*% Matrix::sparseMatrix(i = row(NNarray)[!is.na(NNarray)], j= NNarray[!is.na(NNarray)], x= t(tatato [,1,])[!is.na(NNarray)]) %*% v_right
  )*100000)[1,1])
  print("finite diff sandwich")
  print(tatata[range_col_idx, range_row_idx])
  
  
  t1 = Sys.time()
  tatata = derivative_sandwiches(vecchia = tatato, left_vector = v_left, right_vector = v_right, NNarray = t(NNarray), num_threads = 4, sauce_determinant_chef = T)
  print(paste("computing derivative sandwich with determinant sauce  took", Sys.time()-t1, "s"))
  print("analytical sandwich")
  print(((
    v_left %*% Matrix::sparseMatrix(i = row(NNarray)[!is.na(NNarray)], j= NNarray[!is.na(NNarray)], x= t(tatato_[,1,])[!is.na(NNarray)]) %*% v_right -
      v_left %*% Matrix::sparseMatrix(i = row(NNarray)[!is.na(NNarray)], j= NNarray[!is.na(NNarray)], x= t(tatato [,1,])[!is.na(NNarray)]) %*% v_right 
    + sum(log(tatato_[1,1,])) - sum(log(tatato[1,1,]))
  )*100000)[1,1])
  print("finite diff sandwich with determinant sauce")
  print(tatata[range_col_idx, range_row_idx])
  print("=========================")
}



Rcpp::sourceCpp("src/vecchia.cpp")



# testing efficiency wrt old packagfe Bidart
set.seed(1)
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
# anisotropic
mean(sapply(seq(30), function(i){
  t1 = Sys.time()
  tatato = Bidart::nonstat_vecchia_Linv(
    log_range = (log_range), locs = (locs), NNarray = (NNarray), num_threads = 6, compute_derivative = T, covfun_name = "nonstationary_matern_anisotropic", sphere = F)
  Sys.time()-t1
})) / 
mean(sapply(seq(30), function(i){
  t1 = Sys.time()
  t1 = Sys.time()
  tatato = vecchia(
  log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 6, compute_derivative = T, smoothness = 1.5)
  Sys.time()-t1
}))



