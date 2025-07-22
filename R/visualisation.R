#' Create visualisations for a PP object, according to the vecchia_approx that produced it.
#' 
#' @param PP an object of class PP, create with `createPP`
#' @param vecchia the vecchia approx object that the PP object was created from. 
#' @param separate logical(default to FALSE). Should the plots be printed separatly ?
#'
#' @returns NULL
#' @export
#'
#' @examples
#' vecchia_approx = createVecchia(observed_locs  = cbind(runif(10000), runif(10000)), 10)
#' pepito = createPP(vecchia_approx, plot=FALSE)
#' plot.PP(pepito, vecchia_approx)
plot.PP <- function(PP, vecchia = NULL, mar_var_loss = TRUE, separate=FALSE) {
  def.par <- par(no.readonly = TRUE)
  par(mar=c(3, 3, 3, 1) + 0.5)
  par(mgp=c(2, 1, 0))
  if(mar_var_loss & !separate) {
    layout(matrix(c(1,1,2), 1, 3, byrow = TRUE))
  }
  mar_var <- NULL
  if(mar_var_loss) {
    mar_var = var_loss_percentage.PP(PP)
  }
  plot_knots.PP(x = PP, locs = vecchia_approx$locs, mar_var_loss = mar_var)
  if(mar_var_loss) {
    hist(mar_var,
         xlab = "percentage of lost variance", 
         main = "Histogram of lost marginal\nvariance between the PP and\nthe full GP",
         cex.main=1)
  }
  par(def.par)
}


#' @title Plot the knots and the spatial locations of a PP
#' @param x an object of class \code{PP}
#' @param locs locations
#' @param mar_var_loss optional, the var loss computed by var_loss_percentage.PP
#' @param show_knots logical, default to TRUE. Should the knots be plotted ?
#' @examples
#' observed_locs = cbind(runif(1000), runif(1000))
#' observed_locs = observed_locs[ceiling(nrow(observed_locs)*runif(3000)),]
#' vecchia_approx = createVecchia(observed_locs)
#' pepito = createPP(vecchia_approx)
#' plot_knots.PP(pepito, vecchia_approx$locs)
#' vPP <- var_loss_percentage.PP(pepito)
#' plot_knots.PP(pepito, vecchia_approx$locs, mar_var_loss = vPP)
plot_knots.PP = function(x, locs, mar_var_loss = NULL, show_knots = TRUE, cex = c(.5, .5)){
  # Graphical parameters
  def.par <- par(no.readonly = TRUE)
  omar <- def.par[["mar"]]
  par(mgp=c(2, 1, 0))
  if(omar[2]>=omar[4] & (show_knots | !is.null(mar_var_loss))) {
    par(mar = omar + c(1, 1, 1, 4*1.5), xpd = TRUE)
  }
  
  nlocs <- nrow(locs)
  nknots <- nrow(x$knots)
  legendtitle <- NULL

  
  ### Colors of locations
  locs_col = "#CDCDCDCC"
  if(!is.null(mar_var_loss)) {
    # Remove knots from mar_var_loss
    legendtitle <- "\n\n\nLoss of\nmarginal\nvariance\n(%)"
    cols <- mar_var_loss[-seq(nknots)]
    # Reorder to plot worst last
    ord <- order(cols)
    cols <- cols[ord]
    locs <- locs[ord,]
    cut_var <- cut(cols, breaks = c(0,1,2,5,10,50,100), include.lowest = T)
    base_colors <- c('#F3D97CCC', '#F5C46FCC', '#EDAB65CC', '#DD8A5BCC', '#C2604FCC', '#A02F42CC')
    locs_col = base_colors[as.numeric(cut_var)]
  }
  
  #### Plot locations and knots
  maxlim <- pmax(apply(locs,2, max), apply(x$knots,2, max))
  minlim <- pmin(apply(locs,2, min), apply(x$knots,2, min))
  # 
  # Plot locations
  pch = 16
  if(nlocs>1e6) pch="."
  title = "Locations of PP"
  if(show_knots) title = "Knot placement of PP"
  plot(locs, 
       cex = 1,
       col = locs_col, 
       pch = pch, 
       xlab = "1st spatial coordinate",
       ylab = "2nd spatial coordinate",
       main = title,
       xlim=c(minlim[1], maxlim[2]),
       ylim=c(minlim[2], maxlim[2]),
       cex.main=1
  )
  # Plot knots
  if(show_knots) points(x$knots, pch = 10, cex=1, col=1)
  
  # ### compute size of legend
  mylegend <- legend(x="right",
                     legend="knots",
                     pch =10,
                     col=1,
                     title="Knots",
                     plot = FALSE)
  if(!is.null(mar_var_loss)) {
    ### Marginal variance loss scale
    legend(x="topright",
           legend=levels(cut_var), 
           fill=base_colors, 
           title=legendtitle,
           inset=c(-1.2*mylegend$rect$w, 0),
           bty="n")
  }
  if(show_knots) {
    legend(x="topright",
           legend="knots", 
           pch =10,
           col=1,
           inset=c(-mylegend$rect$w, 0),
           bty="n")
  }
  
  #### Restaure margins
  par(mar=omar, mgp=def.par[["mgp"]])
}

#' Plots range ellipses for nonstationary covariance functions.
#' @param locs ellipses centers. A matrix with 2 columns, 1 row  for each ellipse
#' @param log_range log_range at the ellipse centers, can have 1 or 3 columns. A matrix with 1 row (? TODO) and 1 or 3 columns.
#' @param main character, main title of the plot. Not used if add is TRUE.
#' @param add logical, add on existing plot ? cancels main. Default to FALSE
#' @param shrink numeric, shrinks or inflates the ellipses. shrink = 1 gives the
#' Mahalanobis distance = 1. Shrink = sqrt(8*nu) gives the ellipses corresponding to
#' correlation  = .1 (rho following INLA's terminology)
#'
#' @returns a plot object
#' @export
#'
#' @examples
#' locs <- matrix(rnorm(20), ncol=2)
#' log_range <- matrix(rnorm(30), ncol=3)
#' plot_ellipses(locs, log_range, shrink=0.1)
#'
#' locs <- matrix(rnorm(20), ncol=2)
#' log_range <- matrix(rnorm(10), ncol=1)
#' plot_ellipses(locs, log_range, shrink=0.1)
plot_ellipses = function(locs,
                         log_range,
                         shrink = .1,
                         main = "ellipses",
                         add  = F)
{
  if (ncol(log_range) == 3) {
    #to match parametrization in compute sparse chol
    log_range = log_range %*% matrix(c(1 / sqrt(2), 1 / sqrt(2), 0, 1 /
                                         sqrt(2), -1 / sqrt(2), 0, 0, 0, 1), 3) * sqrt(2)
    matrices = lapply(split(log_range, row(log_range)), expmat)
  }
  if (ncol(log_range) == 1) {
    matrices = lapply(log_range, function(x)
      diag(exp(x#/sqrt(2))
               , 2))
    )
  }
  if (!add) {
    plot(
      locs,
      type = "n",
      xlab = "",
      ylab = "",
      main = main
    )
  }
  for (i in seq(nrow(locs)))
  {
    # 2.447747 must be some bivariate confidence interval
    # shrink = 1 gives the package's Mahalanobis distance
    matrices[[i]] = eigen(matrices[[i]])
    matrices[[i]] = (matrices[[i]])$vec %*% diag(matrices[[i]]$val^2) %*% t(matrices[[i]]$vec) /
      (2.447747)^2
    ell = ellipse::ellipse(matrices[[i]]) * shrink
    ell[, 1] = ell[, 1] + locs[i, 1]
    ell[, 2] = ell[, 2] + locs[i, 2]
    lines(ell)
  }
}


# # test
# nu = .5
# locs = as.matrix(expand.grid(seq(0, 1, .02), seq(0, 1, .02)))
# locs = locs[GpGp::order_maxmin(locs),]
# locs = rbind(c(.501, .501), locs)
# range_beta = .3*matrix(rnorm(3), 1)
# range_beta[1,1] = -3
# NNarray = GpGp::find_ordered_nn(locs, 10)
#
# sparse_chol_aniso =
#   Matrix::sparseMatrix(
#     i = row(NNarray)[!is.na(NNarray)],
#     j = NNarray[!is.na(NNarray)],
#     x =
#       compute_sparse_chol(
#         range_beta = range_beta,
#         NNarray = NNarray, locs = locs,
#         anisotropic = T,
#         sphere = F,
#         PP = NULL, use_PP = F,
#         range_X = matrix(1, nrow(locs)),
#         compute_derivative = F,
#         nu = nu,
#         locs_idx = NULL,
#         num_threads = 10
#       )[[1]][!is.na(NNarray)],
#     triangular = T
#   )
# sparse_chol_iso =
#   Matrix::sparseMatrix(
#     i = row(NNarray)[!is.na(NNarray)],
#     j = NNarray[!is.na(NNarray)],
#     x =
#       compute_sparse_chol(
#         range_beta = range_beta[,1,drop=F],
#         NNarray = NNarray, locs = locs,
#         anisotropic = F,
#         sphere = F,
#         PP = NULL, use_PP = F,
#         range_X = matrix(1, nrow(locs)),
#         compute_derivative = F,
#         nu = nu,
#         locs_idx = NULL,
#         num_threads = 10
#       )[[1]][!is.na(NNarray)],
#     triangular = T
#   )
#
# log_range =
#   (
#     X_PP_mult_right(
#       X = matrix(1, nrow(locs)),
#       PP = 0, use_PP = F,
#       locs_idx = NULL,
#       Y = range_beta
#     )
#   )
#
#
# # tatato = GpGp::fast_Gp_sim(covparms = c(.5, .1, 1.5, 0), locs = locs, m = 10, covfun_name = "matern_isotropic")-4
# # log_range = cbind(tatato, tatato, 0)
#
# z = rnorm(nrow(locs))
# w_aniso = as.vector(Matrix::solve(sparse_chol_aniso, z))
# cor_aniso = Matrix::tcrossprod(Matrix::solve(sparse_chol_aniso))[,1]
# w_iso = as.vector(Matrix::solve(sparse_chol_iso, z))
# cor_iso = Matrix::tcrossprod(Matrix::solve(sparse_chol_iso))[,1]
# plot_pointillist_painting(locs, log_range[,1])
# plot_pointillist_painting(locs, w_aniso)
# plot_pointillist_painting(locs, w_iso)
# plot_pointillist_painting(locs, cor_aniso)
# plot_pointillist_painting(locs, cor_iso)
#
# plot(locs, pch = 16, col = 1+(cor_iso < .1), cex=  .5, main  = "cor = .1 ellipse for anisotropic")
# plot(locs, pch = 16, col = 1+(cor_aniso < .1), cex=  .5)
# plot_ellipses(
#   locs = locs[1,,drop=F], log_range = log_range[1,,drop=F],
#   shrink = sqrt(8*nu), add=  T)
# legend("topleft", legend = c("cor > .1", "cor < .1"), fill = c(1,2))
# plot(locs, pch = 16, col = 1+(cor_iso < .1), cex=  .5, main  = "cor = .1 ellipse for isotropic")
# plot_ellipses(
#   locs = locs[1,,drop=F], log_range = log_range[1,1,drop=F],
#   shrink = sqrt(8*nu), add=  T)
# legend("topleft", legend = c("cor > .1", "cor < .1"), fill = c(1,2))



## test with GpGp, the empirical rho is always greater than theoretical rho
# ra = .1
# locs = cbind(seq(0, 8*ra, .01), 0)
# nu = .5
# cors = GpGp::matern_isotropic(c(1, ra, nu, 0), locs)[,1]
# plot(locs[,1], cors)
# abline(h = .1)
# abline(v = locs[match(T, cors<.1), 1])
# print(paste("empirical rho = ", locs[match(T, cors<.1), 1]))
# print(paste("theorhetical rho = ", round(ra * sqrt(8*nu), 3)))


#' Title get_colors TODO
#'
#' @param x a vector
#'
#' @returns a vector of colors
#' @export
#'
#' @examples
#' get_colors(c(1,2,3))
get_colors = function(x) {
  colors = rep(1, length(x))
  colors[!is.na(x)] = heat.colors(100)[round((x[!is.na(x)] - min(x[!is.na(x)])) /
                                               (max(x[!is.na(x)]) - min(x[!is.na(x)])) * 90) + 1]
  colors
}


#' Plots a spatial variable like a pointillist painting using R base's points. Stupid, but handy.
#'
#' @param locs spatial locations
#' @param field interest variable
#' @param main main title
#' @param cex shrinks or inflates the points
#' @param add logical (default to FALSE)
#'
#' @returns a plot component (a `plot` if `add == FALSE`, `points` if `add ==TRUE`)
#' @export
#'
#' @examples
#' locs <- matrix(rnorm(2000), ncol=2)
#' plot_pointillist_painting(locs=locs, field=locs[,1])
#' plot_pointillist_painting(locs=locs, field=rnorm(1000))
plot_pointillist_painting = function(locs,
                                     field,
                                     cex = 1,
                                     main = NULL,
                                     add = FALSE)
{
  if (!add)
    plot(
      locs,
      col = get_colors(field),
      main = main,
      pch = 15,
      cex = cex,
      xlab  = "",
      ylab = ""
    )
  if (add)
    points(
      locs,
      col = get_colors(field),
      pch = 15,
      cex = cex,
      xlab  = "",
      ylab = ""
    )
}

#' Plot the scale of colors for a given variable
#'
#' @param field interest variable, a vector
#' @param scalename vector, default to ""
#'
#' @returns a plot
#'
#' @examples
#' pointillist_colorscale(field=c(1,2,3))
pointillist_colorscale = function(field, scalename = "") {
  origmar <- par("mar")
  par(mar = c(0, 0, 0, 0))
  plot(
    0,
    0,
    type = "n",
    xaxt = 'n',
    yaxt = 'n',
    xlim = c(-.007, .005),
    ylim = c(0, 100),
    frame = FALSE,
    xlab = "",
    ylab = ""
  )
  barplot(
    rep(.005, 100),
    col = (heat.colors(101)),
    width = rep(1, 100),
    space = 0,
    xlab = scalename,
    ylab = "",
    main = "",
    border = T,
    horiz = T,
    add = T,
    axes = F
  )
  text(y = 99, -.004, format(signif(max(field)), trim = 4, width = 4))
  text(y = 75, -.004, format(signif(min(field) + (
    max(field) - min(field)
  ) * 3 / 4), trim = 4, width = 4))
  text(y = 50, -.004, format(signif(min(field) + (
    max(field) - min(field)
  ) / 2), trim = 4, width = 4))
  text(y = 25, -.004, format(signif(min(field) + (
    max(field) - min(field)
  ) * 1 / 4), trim = 4, width = 4))
  text(y = 1, -.004, format(signif(min(field)), trim = 4, width = 4))
  par(mar = origmar)
}


#' Cumulative sum of an array
#'
#' @param x an array of more than 2 dimensions
#'
#' @returns an array, the same dimension as x
#' @export
#'
#' @examples
#' myarray <- array(abs(rnorm(200)), dim = c(10, 10, 2))
#' array_cumsum(myarray)
array_cumsum = function(x)
{
  res = array(0, dim = dim(x))
  for (i in seq(dim(x)[1]))
  {
    for (j in seq(dim(x)[2]))
    {
      res[i, j, ] = cumsum(x[i, j, ])
    }
  }
  res
}


#' Title TODO
#'
#' @param x  a 3 dimensional array
#' @param M a Matrix (TODO)
#'
#' @returns an 3 dimensional array, of dimension `(dim(x)[1], dim(x)[2], M@Dim[2])`
#' @export
#'
#' @examples
#' x <- array(1:12, dim = c(2, 2, 3))
#' M <- Matrix::Matrix(matrix(1:9, nrow = 3, ncol = 3), sparse = TRUE)
#' result <- array_multiply_3(x, M)
array_multiply_3 = function(x, M)
{
  res = array(0, dim = c(dim(x)[c(1, 2)], M@Dim[2]))
  for (i in seq(dim(res)[1]))
  {
    for (j in seq(dim(res)[2]))
    {
      res[i, j, ] = as.vector(x[i, j, ] %*% M)
    }
  }
  res
}

#' Title TODO
#'
#' @param x an array of 3 dimensions
#' @param M a sparse matrix ? (TODO)
#'
#' @returns a 3 dimansional array of dimensions `(dim(x)[1], dim(M)[2], dim(x)[3])`
#' @export
#'
#' @examples
#' x <- array(1:24, dim = c(2, 3, 4))
#' M <- matrix(1:6, nrow = 2, ncol = 3)
#' #result <- array_multiply_2(x, M)
array_multiply_2 = function(x, M)
{
  # res = array(0, dim = c(dim(x)[1], dim(M)[1], dim(x)[3])) # TODO ? ça fail sinon
  res = array(0, dim = c(dim(x)[1], dim(M)[2], dim(x)[3]))
  for (i in seq(dim(res)[1]))
  {
    for (j in seq(dim(res)[3]))
    {
      res[i, , j] = as.vector(M %*% x[i, , j])
    }
  }
  res
}

#' Multiply a matrix M with each column of a 3 dimensional array.
#'
#' @param x an 3 dimensional array
#' @param M a matrix
#'
#' @returns a 3 dimensional array of dimensions `(dim(M)[1], dim(x)[2], dim(x)[3]).`
#' @export
#'
#' @examples
#' M <- matrix(1:6, nrow = 3, ncol = 2)
#' x <- array(1:24, dim = c(2, 3, 4))
#' result <- array_multiply_1(x, M)
array_multiply_1 = function(x, M)
{
  res = array(0, dim = c(dim(M)[1], dim(x)[2], dim(x)[3]))
  for (i in seq(dim(res)[2]))
  {
    for (j in seq(dim(res)[3]))
    {
      res[, i , j] = as.vector(M %*% x[, i, j])
    }
  }
  res
}


#' Title TODO ? Retourne les valeurs nécessaires pour les diagnostics ?
#'
#' @param record_arrays a list of arrays, each of dimension 3
#' @param iterations TODO
#' @param burn_in a numeric, default to 0.5 TODO
#' @param starting_proportion a numeric, default to 0.5 TODO
#'
#' @returns a list
#' @export
#'
#' @examples
#' myarrays <- list(
#'   array(rnorm(100000)+10, dim = c(10, 10, 1000)),
#'   array(rnorm(100000)+10, dim = c(10, 10, 1000)),
#'   array(rnorm(100000)+10, dim = c(10, 10, 1000))
#' )
#' test <- grb_diags_field(myarrays, iterations = seq(1000)*10)
grb_diags_field = function(record_arrays,
                           iterations,
                           burn_in = .5,
                           starting_proportion = .5)
{
  cumsums = lapply(record_arrays, array_cumsum)
  sqcumsums = lapply(record_arrays, function(x)
    array_cumsum(x^2))
  iter_start_idx = which(iterations > iterations[length(iterations)] * starting_proportion)[1]
  diff_array = cbind(0, seq(iter_start_idx, length(iterations)))
  n = rep(0, nrow(diff_array))
  lower_bound_idx = 1
  for (i in seq(nrow(diff_array)))
  {
    while (iterations[lower_bound_idx] < iterations[diff_array[i, 2]] * burn_in) {
      lower_bound_idx = lower_bound_idx + 1
    }
    diff_array[i, 1] = lower_bound_idx
    n[i] = i + iter_start_idx - lower_bound_idx - 1
  }
  diff_matrix = Matrix::sparseMatrix(
    j = c(seq(nrow(diff_array)), seq(nrow(diff_array))),
    i = c(diff_array),
    x = rep(c(-1, 1), each = nrow(diff_array)),
  )
  mean_estimators = lapply(cumsums, function(x)
    array_multiply_3(x, M = diff_matrix %*% Matrix::Diagonal(x = 1 / n)))
  sqmean_estimators = lapply(sqcumsums, function(x)
    array_multiply_3(x, M = diff_matrix %*% Matrix::Diagonal(x = 1 / n)))
  var_estimators = mapply(
    function(x, y)
      array_multiply_3(x - y^2, M = Matrix::Diagonal(x = n / (n - 1))),
    x = sqmean_estimators,
    y = mean_estimators,
    SIMPLIFY = F
  )
  within_mean = Reduce("+", mean_estimators) / length(mean_estimators)
  within_var = Reduce("+", var_estimators) / length(var_estimators)
  between_var = (Reduce("+", lapply(mean_estimators, function(x)
    x^2)) / length(mean_estimators) - within_mean^2) * length(mean_estimators) / (length(mean_estimators) -
                                                                                    1)
  PSRF =   (length(var_estimators) + 1) / (length(var_estimators)) * array_multiply_3(x = between_var /
                                                                                        within_var, M = Matrix::Diagonal(x = (n - 1) / n))
  for (i in seq(dim(PSRF)[3]))
  {
    PSRF[, , i] = PSRF[, , i] + (n[i] + 1) / n[i]
  }
  PSRF_quantiles  = apply(PSRF, 3, function(x)
    quantile(x, probs = c(1, .99, .9, .5), na.rm = T))
  list(
    "iterations" = iterations[diff_array[, 2]],
    "mean" = mean_estimators,
    "var" = var_estimators,
    "PSRF" = PSRF,
    "PSRF_quantiles" = PSRF_quantiles
  )
}


#' Title TODO
#'
#' @param PSRF TODO
#' @param individual_varnames TODO !! NOT USED IN CODE
#' @param varname name of the variable
#'
#' @returns a plot
#' @export
#'
#' @examples
#' myarrays = list(
#'   array(rnorm(600), dim = c(4, 2, 100)),
#'   array(rnorm(600), dim = c(4, 2, 100)),
#'   array(rnorm(600), dim = c(4, 2, 100))
#' )
#' test = grb_diags_field(record_arrays = myarrays, iterations = seq(100), burn_in = .5, starting_proportion = .5)
#' plot_PSRF(test, varname = "Example")
plot_PSRF = function(PSRF,
                     individual_varnames = NULL,
                     varname = "")
{
  if (any(is.infinite(PSRF$PSRF)) | any(is.na(PSRF$PSRF)))
  {
    plot(
      0,
      0,
      type = "n",
      main = paste(
        "PSRF of",
        varname,
        "not represented because of infinite GRB diags"
      ),
      xlab = "iterations",
      ylab = "PSRF quantiles"
    )
    return()
  }
  
  plot(
    PSRF$iterations,
    PSRF$PSRF_quantiles[1, ],
    ylim = c(1, max(PSRF$PSRF_quantiles[1, ])),
    main = paste("Quantiles of PSRF of", varname),
    type = "l",
    xlab = "iterations",
    ylab = "PSRF quantiles",
    log = "y"
  )
  for (i in seq(1, dim(PSRF$PSRF)[1]))
  {
    for (j in seq(1, dim(PSRF$PSRF)[2]))
    {
      lines(PSRF$iterations,
            PSRF$PSRF[i, j, ],
            col = scales::alpha("lightgray", .4),
      )
    }
  }
  lines(PSRF$iterations, PSRF$PSRF_quantiles[2, ], col = 2)
  lines(PSRF$iterations, PSRF$PSRF_quantiles[3, ], col = 4)
  lines(PSRF$iterations, PSRF$PSRF_quantiles[4, ], col = 6)
  legend(
    bg = "white",
    "topleft",
    legend = c("max", .99, .9, .5),
    fill = c(1, 2, 4, 6)
  )
  
  abline(h = 1)
}

#' Title TODO
#'
#' @param log_scale_arrays list of arrays TODO
#' @param iterations a vector TODO
#' @param starting_proportion numeric value, default to 0.5 TODO
#' @param varname name of the variable TODO
#'
#' @returns a plot
#' @export
#'
#' @examples
#' arrays = list(
#'    array(rnorm(6000)+ rep(c(10,10,0,10,0,0), 1000), dim = c(6, 1, 1000)),
#'    array(rnorm(6000)+ rep(c(10,10,0,10,0,0), 1000), dim = c(6, 1, 1000)),
#'    array(rnorm(6000)+ rep(c(10,10,0,10,0,0), 1000), dim = c(6, 1, 1000))
#' )
#' plot_log_scale(arrays, iterations = seq(1000), starting_proportion = .5, varname="Example")
plot_log_scale = function(log_scale_arrays,
                          iterations,
                          starting_proportion = .5,
                          varname)
{
  kept_iterations = which(iterations > starting_proportion * iterations[length(iterations)])
  if (dim(log_scale_arrays[[1]])[1] == 1 &
      dim(log_scale_arrays[[1]])[2] == 1)
  {
    log_scale = sapply(log_scale_arrays, function(x)
      x[, , kept_iterations])
    plot(
      iterations[kept_iterations],
      iterations[kept_iterations],
      type = "n",
      xlab = "iteration",
      ylab = "log scale",
      main = varname,
      ylim = c(min(unlist(log_scale)), max(unlist(log_scale)))
    )
    # loop over chains
    for (i in seq(ncol(log_scale)))
    {
      lines(iterations[kept_iterations], log_scale[, i])
    }
  }
  if (dim(log_scale_arrays[[1]])[1] == 6)
  {
    log_scale = sapply(log_scale_arrays, function(x)
      x[, , kept_iterations])
    #marginal_logvars = lapply(log_scale_arrays, function(x)
    #  apply(x[,,kept_iterations], 2, function(x) diag(
    #    rbind(c(1/sqrt(2), 1/sqrt(2), 0), c(1/sqrt(2), -1/sqrt(2), 0), c(0, 0, 1)) %*%
    #      symmat(x) %*%
    #      cbind(c(1/sqrt(2), 1/sqrt(2), 0), c(1/sqrt(2), -1/sqrt(2), 0), c(0, 0, 1))
    #    ))
    #)
    marginal_logvars = lapply(log_scale_arrays, function(x)
      apply(x[, , kept_iterations], 2, function(x)
        diag(symmat(x))))
    plot(
      iterations[kept_iterations],
      iterations[kept_iterations],
      type = "n",
      xlab = "iteration",
      ylab = "log scale",
      main = paste(varname, "split by component"),
      ylim = c(min(unlist(
        marginal_logvars
      )), max(unlist(
        marginal_logvars
      )))
    )
    for (i in seq(length(marginal_logvars)))
    {
      for (j in seq(nrow(marginal_logvars[[i]])))
      {
        lines(iterations[kept_iterations], marginal_logvars[[i]][j, ], col = c(1, 2, 4)[j])
      }
    }
    legend(
      bg = "white",
      "topleft",
      legend = c("Determinant", "Anisotropy", "Anisotropy"),
      fill = c(1, 2, 4)
    )
  }
}


#' Title TODO
#'
#' @param beta_arrays list of arrays TODO
#' @param iterations numerical vector TODO
#' @param starting_proportion numeric value, default to 0.5 TODO
#' @param varname name of variable
#' @param var_names ? TODO
#'
#' @returns several plots
#' @export
#'
#' @examples
#' beta_arrays = list(
#'   array(rnorm(400), dim = c(2, 2, 100)),
#'   array(rnorm(400), dim = c(2, 2, 100)),
#'   array(rnorm(400), dim = c(2, 2, 100))
#' )
#' plot_beta(beta_arrays, seq(100), starting_proportion = .5, varname = "Example", var_names = c(1, 2))
#' beta_arrays = lapply(seq(3), function(i){
#'    res = array(data = 0, dim = c(10, 3, 100))
#'    res[,1,] = rnorm(length(res[,1,]))
#'    res
#'  })
#' plot_beta(beta_arrays, seq(100), starting_proportion = .5, varname = "Example", var_names = c(1, 2))
plot_beta = function(beta_arrays,
                     iterations,
                     starting_proportion = .5,
                     varname,
                     var_names = NULL)
{
  #if(dim(beta_arrays[[1]])[2]==3)beta_arrays = lapply(beta_arrays, function(x)array_multiply_2(x, matrix(c(1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2), -1/sqrt(2), 0, 0, 0, 1), 3)))
  kept_iterations = which(iterations > starting_proportion * iterations[length(iterations)])
  for (i in seq(dim(beta_arrays[[1]])[2]))
  {
    upper_window = max(sapply(beta_arrays, function(x)
      max(x[, i, kept_iterations])))
    lower_window = min(sapply(beta_arrays, function(x)
      min(x[, i, kept_iterations])))
    if (dim(beta_arrays[[1]])[2] == 3)
    {
      if (i == 1)
        main = paste("Determinant component of", varname)
      if (i != 1)
        main = paste("anisotropy component", i - 1, "of", varname)
    } else
      main = varname
    plot(
      iterations[kept_iterations],
      iterations[kept_iterations],
      type = "n",
      main = main,
      ylab = "parameter value",
      xlab = "iteration",
      ylim = c(lower_window, upper_window)
    )
    for (j in seq(dim(beta_arrays[[1]])[1])) {
      for (x in beta_arrays)
        lines(iterations[kept_iterations], x[j, i, kept_iterations], col = c(1, 3, 4, 6, 7, 8)[(j) %%
                                                                                                 6])
    }
    if (!is.null(var_names)) {
      legend(
        bg = "white",
        "topleft",
        legend = var_names,
        fill = c(1, 3, 4, 6, 7, 8)[(seq(dim(beta_arrays[[1]])[1])) %% 6]
      )
    }
  }
}


#' Represents the samples and the Gelman-Rubin-Brooks diagnostics.
#' The proportion of iterations used for the plotting and the proportion of the burn-in are adjustable. They multiply.
#' Note that GRB curves are given in order to avoid spurious green lights.
#' @param mcmc_nngp_list a mcmc_nngp_list created using mcmc_nngp_isitialize and run using mcmc_nngp_run
#' @param plot_PSRF_fields logical(default to FALSE) TODO
#' @param burn_in between 0.01 and .99, the proportion of samples discarded for the burn-in
#' @param starting_proportion between 0.01 and .99, the proportion of iterations that is used
#' @export
#'
#' @returns plots TODO
#' @export
#'
#' @examples
#' \dontrun{TODO}
diagnostic_plots = function(mcmc_nngp_list,
                            plot_PSRF_fields = FALSE,
                            burn_in = .5,
                            starting_proportion = .5)
{
  records_names = names(mcmc_nngp_list$records$chain_1)
  if (!plot_PSRF_fields)
    records_names = records_names[-grep("field", records_names)]
  range_names = records_names[grep("range", records_names)]
  records_names = setdiff(records_names, range_names)
  noise_names = records_names[grep("noise", records_names)]
  records_names = setdiff(records_names, noise_names)
  scale_names = records_names[grep("scale", records_names)]
  records_names = setdiff(records_names, scale_names)
  response_names = records_names
  if (length(mcmc_nngp_list$records) > 1) {
    for (i in list(range_names, noise_names, scale_names, response_names))
    {
      par(mfrow = c(1, 1))
      for (j in i)
      {
        PSRF = grb_diags_field(
          record_arrays = lapply(mcmc_nngp_list$records, function(x)
            x[[j]]),
          iterations = mcmc_nngp_list$iterations$thinning,
          burn_in = burn_in,
          starting_proportion = starting_proportion
        )
        plot_PSRF(
          PSRF = PSRF,
          varname = j,
          individual_varnames = row.names(mcmc_nngp_list$states$chain_1$params[[j]])
        )
      }
    }
  }
  for (i in list(range_names, noise_names, scale_names))
  {
    #print(i)
    #par(mfrow = c(length(grep("log_scale", i)) + dim(mcmc_nngp_list$records$chain_1 [[i[grep("beta", i)]]])[2], 1))
    par(mfrow = c(1, 1))
    if (length(grep("log_scale", i)) > 0)
      plot_log_scale(
        log_scale_arrays = lapply(mcmc_nngp_list$records, function(x)
          x[[i[grep("log_scale", i)]]]),
        iterations = mcmc_nngp_list$iterations$thinning,
        starting_proportion = starting_proportion,
        varname = i[grep("log_scale", i)]
      )
    plot_beta(
      beta_arrays = lapply(mcmc_nngp_list$records, function(x)
        x[[i[grep("beta", i)]]]),
      iterations = mcmc_nngp_list$iterations$thinning,
      starting_proportion = starting_proportion,
      varname = i[grep("beta", i)],
      var_names = row.names(mcmc_nngp_list$states$chain_1$params[[i[grep("beta", i)]]])
    )
  }
}


#' Prints the Effective Sample Size of a MCMC run.
#' The proportion of the burn-in is adjustable.
#' @param mcmc_nngp_list a mcmc_nngp_list created and run by the package
#' @param burn_in between 0.01 and .99, the proportion of samples discarded for the burn-in
#'
#' @returns a list
#' @export
#'
#' @examples
#' \dontrun{TODO}
ESS = function(mcmc_nngp_list, burn_in = .5) {
  iterations = mcmc_nngp_list$iterations
  iter_start_idx = match(TRUE,
                         iterations$thinning > (iterations$checkpoints[nrow(iterations$checkpoints), 1] *
                                                  burn_in))
  varnames = c("range_beta", "scale_beta", "noise_beta", "beta")
  if ("range_log_scale" %in% names(mcmc_nngp_list$records$chain_1))
    varnames = c(varnames, "range_log_scale")
  if ("noise_log_scale" %in% names(mcmc_nngp_list$records$chain_1))
    varnames = c(varnames, "noise_log_scale")
  if ("scale_log_scale" %in% names(mcmc_nngp_list$records$chain_1))
    varnames = c(varnames, "scale_log_scale")
  ESSs = lapply(varnames, function(name) {
    res = as.matrix(Reduce("+", lapply(mcmc_nngp_list$records, function(record)
      apply(record[[name]][, , -seq(iter_start_idx), drop = F], c(1, 2), function(x)
        coda::effectiveSize(c(x))))))
    row.names(res) = row.names(mcmc_nngp_list$states$chain_1$params[[name]])
    res
  })
  names(ESSs) = varnames
  ESSs
}
