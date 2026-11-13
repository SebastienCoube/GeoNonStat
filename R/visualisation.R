#' Create visualisations for a PP object, according to the vecchia_approx that produced it.
#' @param x an object of class PP, create with `createPP`
#' @param ... additional arguments (unused)
#' @param mar_var_loss boolean, default to TRUE. Should loss of margine variance be computed and plotted ?
#' @param samples logical, default to TRUE. Also show samples that could be produced by the PP ? 
#' @param separate logical(default to FALSE). Should the plots be printed separatly ?
#'
#' @returns NULL
#' @rdname PP
#' @aliases PP
#' @method plot PP
#' @export
#'
#' @examples
#' vecchia_approx <- createVecchia(observed_locs = cbind(runif(10000), runif(10000)), 10)
#' pepito <- createPP(vecchia_approx, plot = FALSE)
#' plot(pepito)
plot.PP <- function(x,
                    ...,
                    mar_var_loss = TRUE,
                    samples = TRUE,
                    separate = FALSE) {
  oldpar <- par(c("mar", "mgp", "mfrow"))
  on.exit(par(mar=oldpar$mar, mgp=oldpar$mgp, mfrow=oldpar$mfrow))
  par(mar = c(3, 3, 3, 1) + 0.5)
  par(mgp = c(2, 1, 0))
  if (mar_var_loss & !separate) {
    layout(matrix(c(1, 1, 2), 1, 3, byrow = TRUE))
  }
  mar_var <- NULL
  if (mar_var_loss) {
    mar_var <- varLossPP(x)
  }
  plotKnotsPP(x = x, mar_var_loss = mar_var)
  if (mar_var_loss) {
    hist(mar_var,
      xlab = "percentage of lost variance",
      main = "Histogram of lost marginal\nvariance between the PP and\nthe full GP",
      cex.main = 1
    )
  }
  if(samples){
    if(!separate) par(mfrow = c(2,2))
    plotPointillistPainting(
      x$vecchia_locs,  legend = F, main = "Example of PP sample (1/4)",
      Matrix::solve(x$sparse_chol, c(rnorm(x$n_knots), rep(0, nrow(x$vecchia_locs))))[-seq(x$n_knots)]
    )
    plotPointillistPainting(
      x$vecchia_locs,  legend = F, main = "Example of PP sample (2/4)",
      Matrix::solve(x$sparse_chol, c(rnorm(x$n_knots), rep(0, nrow(x$vecchia_locs))))[-seq(x$n_knots)]
    )
    plotPointillistPainting(
      x$vecchia_locs,  legend = F, main = "Example of PP sample (3/4)",
      Matrix::solve(x$sparse_chol, c(rnorm(x$n_knots), rep(0, nrow(x$vecchia_locs))))[-seq(x$n_knots)]
    )
    plotPointillistPainting(
      x$vecchia_locs,  legend = F, main = "Example of PP sample (4/4)",
      Matrix::solve(x$sparse_chol, c(rnorm(x$n_knots), rep(0, nrow(x$vecchia_locs))))[-seq(x$n_knots)]
    )
  }
}


#' @title Plot the knots and the spatial locations of a PP
#' @param x an object of class PP, create with `createPP`
#' @param mar_var_loss optional, the var loss computed by `varLossPP`
#' @param show_knots logical, default to TRUE. Should the knots be plotted ?
#' @export
#' @keywords internal
#'
#' @examples
#' observed_locs <- cbind(runif(1000), runif(1000))
#' observed_locs <- observed_locs[ceiling(nrow(observed_locs) * runif(3000)), ]
#' vecchia_approx <- createVecchia(observed_locs)
#' pepito <- createPP(vecchia_approx, plot = FALSE)
#' plotKnotsPP(pepito)
plotKnotsPP <- function(x,
                        mar_var_loss = NULL,
                        show_knots = TRUE) {
  # Graphical parameters
  oldpar <- par(c("mar", "mgp", "xpd"))
  on.exit(par(mar = oldpar$mar, mgp = oldpar[["mgp"]], xpd=oldpar$xpd))
  
  omar <- oldpar$mar
  par(mgp = c(2, 1, 0))
  if (omar[2] >= omar[4] & (show_knots | !is.null(mar_var_loss))) {
    par(mar = omar + c(1, 1, 1, 4 * 1.5), xpd = TRUE)
  }

  locs <- x$vecchia_locs
  nlocs <- nrow(locs)
  nknots <- nrow(x$knots)
  legendtitle <- NULL


  ### Colors of locations
  locs_col <- "#CDCDCDCC"
  if (!is.null(mar_var_loss)) {
    # Remove knots from mar_var_loss
    legendtitle <- "\n\n\nLoss of\nmarginal\nvariance\n(%)"
    cols <- mar_var_loss[-seq_len(nknots)]
    # Reorder to plot worst last
    ord <- order(cols)
    cols <- cols[ord]
    locs <- locs[ord, ]
    # Get
    cut_var <- cut(cols,
      breaks = c(0, 1, 2, 5, 10, 50, 100),
      include.lowest = TRUE
    )
    base_colors <- getColors(1:6, alpha = TRUE)
    locs_col <- base_colors[as.numeric(cut_var)]
  }

  #### Plot locations and knots
  maxlim <- pmax(apply(locs, 2, max), apply(x$knots, 2, max))
  minlim <- pmin(apply(locs, 2, min), apply(x$knots, 2, min))
  #
  # Plot locations
  pch <- 16
  if (nlocs > 1e6) {
    pch <- "."
  }
  title <- "Locations of PP"
  if (show_knots) {
    title <- "Knot placement of PP"
  }
  plot(
    locs,
    cex = 1,
    col = locs_col,
    pch = pch,
    xlab = "1st spatial coordinate",
    ylab = "2nd spatial coordinate",
    main = title,
    xlim = c(minlim[1], maxlim[2]),
    ylim = c(minlim[2], maxlim[2]),
    cex.main = 1
  )

  # ### compute size of legend
  mylegend <- legend(
    x = "right",
    legend = "knots",
    title = "Knots",
    plot = FALSE
  )
  if (!is.null(mar_var_loss)) {
    ### Marginal variance loss scale
    legend(
      x = "topright",
      legend = levels(cut_var),
      fill = base_colors,
      title = legendtitle,
      inset = c(-1.2 * mylegend$rect$w, 0),
      bty = "n"
    )
  }
  # Plot knots
  if (show_knots) {
    points(x$knots,, pch = 10, cex = 1, col = 1)
    legend(
      x = "topright",
      legend = "knots",
      pch = 10,
      col = 1,
      inset = c(-mylegend$rect$w, 0),
      bty = "n"
    )
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
#         anisotropic = TRUE,
#         sphere = F,
#         PP = NULL, use_PP = F,
#         range_X = matrix(1, nrow(locs)),
#         compute_derivative = F,
#         nu = nu,
#         locs_idx = NULL,
#         num_threads = 10
#       )[[1]][!is.na(NNarray)],
#     triangular = TRUE
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
#     triangular = TRUE
#   )
#
# log_range =
#   (
#     xPPMultRight(
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
# plotPointillistPainting(locs, log_range[,1])
# plotPointillistPainting(locs, w_aniso)
# plotPointillistPainting(locs, w_iso)
# plotPointillistPainting(locs, cor_aniso)
# plotPointillistPainting(locs, cor_iso)
#
# plot(locs, pch = 16, col = 1+(cor_iso < .1), cex=  .5, main  = "cor = .1 ellipse for anisotropic")
# plot(locs, pch = 16, col = 1+(cor_aniso < .1), cex=  .5)
# plotEllipses(
#   locs = locs[1,,drop=FALSE], log_range = log_range[1,,drop=FALSE],
#   shrink = sqrt(8*nu), add=  TRUE)
# legend("topleft", legend = c("cor > .1", "cor < .1"), fill = c(1,2))
# plot(locs, pch = 16, col = 1+(cor_iso < .1), cex=  .5, main  = "cor = .1 ellipse for isotropic")
# plotEllipses(
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


#' Function to get a color vector to map values of a numerical vector. 
#'
#' @param x a numerical vector
#' @param alpha logical, default to TRUE. Add transparency to colors ?
#'
#' @returns a vector of colors of length `length(x)`
#' @export
#' @keywords internal
#'
#' @examples
#' getColors(c(1, 2, 3))
getColors <- function(x, alpha = FALSE) {
  colors <- rep(1, length(x))
  xnoNA <- x[!is.na(x)]
  col1 <- "#B10026"
  col2 <- "#FFFFB2"
  if (alpha) {
    col1 <- paste0(col1, "CC")
    col2 <- paste0(col2, "CC")
  }
  cols <- colorRampPalette(c(col1, col2))(100)
  colors[!is.na(x)] <- cols[round((xnoNA - min(xnoNA)) /
    (max(xnoNA) - min(xnoNA)) * 99) + 1]
  return(colors)
}

#' Colors for a categorical variable
#'
#' @param n integer, number of colors to get
#' @param alpha logical, default to TRUE. Add transparency to colors ?
#'
#' @returns a vector of colors of length `length(x)`
#' @export
#' @keywords internal
#'
#' @examples
#' getColorsCat(3)
getColorsCat <- function(n, alpha = FALSE) {
  if (n > 9) {
    cols <- c("#000000", grDevices::rainbow(n - 1))
  } else {
    base_colors <- c("#000000", "#E73F74", "#F1CE63", "#77AADD", "#009988", "#9467BD", "#FF9D9A", "#99DDFF", "#AAAA00")
    cols <- base_colors[seq_len(n)]
  }
  if (alpha) {
    cols <- paste0(cols, "CC")
  }
  return(cols)
}

#' Plots a spatial variable like a \href{https://en.wikipedia.org/wiki/Pointillism}{pointillist painting}, the color intensity being the variable value.
#'
#' @param locs numeric matrix of spatial locations
#' @param field numerical vector, interest variable to define color of points
#' @param add logical, default to FALSE. Add on existing plot ?
#' @param alpha logical, default to TRUE. Add transparency to colors ?
#' @param legend logical, default to TRUE. Add legend ?
#' @param ... additional graphical parameters, sent to `plot` or,
#' if `add=TRUE`, to `points`
#' @family {Visualisations}
#' @returns a plot component (a `plot` if `add == FALSE`, `points` if `add ==TRUE`)
#' @export
#'
#' @examples
#' locs <- matrix(runif(2000), ncol = 2)
#' plotPointillistPainting(locs = locs, field = rnorm(1000))
plotPointillistPainting <- function(locs,
                                    field,
                                    add = FALSE,
                                    alpha = TRUE,
                                    legend = TRUE,
                                    ...) {
  args <- list(...)
  if (is.null(args$pch))  args[["pch"]]  <- "."
  if (is.null(args$cex))  args[["cex"]]  <- 4
  if (is.null(args$xlab)) args[["xlab"]] <- NA
  if (is.null(args$ylab)) args[["ylab"]] <- NA
  args$col <- getColors(field, alpha = alpha)
  
  if (add) {
    do.call("points", c(list(x = locs), args))
    return(invisible(NULL))
  }
  
  # margin for legend, if needed
  par(mar = c(2, 2, 4, if (legend) 8 else 0))
  do.call("plot", c(list(x = locs), args))
  
  if (legend) pointillistColorscale(field)
  invisible(NULL)
}

#' Plot the scale of colors for a given variable
#'
#' @param field numerical vector, interest variable to define color of points
#' @param scalename vector, default to ""
#'
#' @returns a plot
#' @export
#' @keywords internal
#'
#' @examples
#' pointillistColorscale(field = c(1, 2, 3))
pointillistColorscale <- function(field, scalename = "",
                                  n = 50,
                                  rel_x = c(0.01, 0.07)) {
  field <- field[!is.na(field)]
  if (length(field) == 0) stop("only NA values in field")
  
  usr <- par("usr")
  dx  <- usr[2] - usr[1]
  x0  <- usr[2] + rel_x[1] * dx
  x1  <- usr[2] + rel_x[2] * dx
  ys  <- seq(usr[3], usr[4], length.out = n + 1)
  
  old_xpd <- par(xpd = NA)          # to draw in margin
  on.exit(par(old_xpd))

  # Color bar
  rect(x0, ys[-(n + 1)], x1, ys[-1], col = getColors(1:n), border = 0.3)
  
  # text
  pos    <- seq(0, 1, by = 0.25)
  values <- stats::quantile(field, probs = pos)
  text(x1 + 0.005 * dx,
       usr[3] + pos * (usr[4] - usr[3]),
       labels = round(values, 3),
       adj = c(0, 0.5), cex = 0.8)
  
  # scale name
  if (nzchar(scalename))
    text((x0 + x1) / 2, usr[3] - 0.03 * (usr[4] - usr[3]),
         labels = scalename, adj = c(0.5, 1), cex = 0.8)
}



#' Plots range ellipses for nonstationary covariance functions. Default is 0.1 correlation.
#' @param object a GeoNonstat object
#' @param prediction_block block of prediction on wich to draw ellipses, usually an element of 
#' the output of `predict.GeoNonStat`
#' @param correlation_level numeric value, between 0 and 1, resize the ellipses according 
#' to correlation of the spatial kernel. Default to .1
#' @family {Visualisations}
#' @returns a plot object
#' @export
#'
#' @examples
#' # Recreating new data with less observations but same structure as example data
#' locs <-processedGnsDemo$vecchia_approx$observed_locs
#' X <- processedGnsDemo$covariates$X$arg
#' suppressWarnings(
#'   testdata <- splitData(locs,  list("X" = X, "noise_X"=X), plot=FALSE)
#' )
#' testpredictions <- predict(processedGnsDemo, testdata)
#' plotPointillistPainting(locs = locs, field = X[,1], cex=20)
#' drawEllipses(processedGnsDemo, testpredictions$test_close, correlation_level=0.7)
drawEllipses = function(object, prediction_block, correlation_level = .1)
{
  log_range <- cbind(
    prediction_block$summaries$range[,"median"],
    prediction_block$summaries$aniso1[,"median"],
    prediction_block$summaries$aniso2[,"median"]
  )
  if(ncol(log_range) ==3){
    #to match parametrization in compute sparse chol
    log_range = log_range %*% matrix(c(1 / sqrt(2), 1 / sqrt(2), 0, 1 / sqrt(2), -1 / sqrt(2), 0, 0, 0, 1), 3) * sqrt(2) * 2
  } else if (ncol(log_range) == 1) {
    log_range <- log_range * 2
  }
  
  range <- t(expRanges(t(log_range)))
  matrices = lapply(split(range, row(range)), function(x){
    res <- matrix(0,2,2)
    res[1] <- x[1];
    res[4] <- x[2]; 
    res[2] <- res[3] <- x[3]
    res
  })
  locs <- prediction_block$locs
  shrink <-
    stats::uniroot(
      f = function(x) maternThingy_(x, 
                                    smoothness15 = object$hierarchical_model$matern_smoothness == 1.5) - 
        correlation_level,
      lower = 0, 
      upper = 1e12)$root
  for(i in seq(nrow(locs))){
    ell = ellipseMahalanobis(matrices[[i]], shrink^2)
    ell[,1] = ell[,1]+locs[i, 1]
    ell[,2] = ell[,2]+locs[i, 2]
    lines(ell)
  }
}


ellipseMahalanobis <- function(Sigma, x) {
  eig   <- eigen(Sigma, symmetric = TRUE)
  A     <- eig$vectors %*% diag(sqrt(eig$values * x), nrow = 2)
  
  theta  <- seq(0, 2 * pi, length.out = 200)
  circle <- rbind(cos(theta), sin(theta))      # 2 x n
  
  pts <- t(A %*% circle)
  colnames(pts) <- c("y1", "y2")
  pts
}


#' Plot the prediction of a 'GeoNonStat' object on new locations
#'
#' @param prediction_block an element of the output of the 
#' `predict.GeoNonStat()` function. 
#' @param which_parameter parameter to plot, default to NULL.
#' Should be one of `names(prediction_block)`
#' @param which_summary summary to plot, default to "mean". On of the available 
#' statistics of the summaries of the prediction object : 
#' "mean", "sd", "median", "q 2.5%" or "q 97.5%".
#' @param main character, title of the plot. 
#' @family {Visualisations}
#' @returns a plot
#' @export
#' 
#' @examples
#' # Recreating new data with less observations but same structure as example data
#' locs <- cbind(rnorm(200), rnorm(200))
#' X <- data.frame(matrix(rnorm(200*4), ncol=4))
#' colnames(X) <- paste0("V", 1:4)
#' testdata <- list("locs" = locs, "X" = X, noise_X=X)  
#' testpredictions <- predict(processedGnsDemo, list("testdata" = testdata))
#' plotGeoNonStatPrediction(testpredictions$testdata,  "field", "mean")
plotGeoNonStatPrediction <- function(
    prediction_block, which_parameter, which_summary, main = "automatic"){
  sums <- prediction_block$summaries
  if(is.null(which_parameter) || !(which_parameter %in% names(sums))) {
    stop("`which_parameter` should be on of ", paste(names(sums), collapse=", "))
  }
  if(is.null(which_summary) || !(which_summary %in% colnames(sums[[which_parameter]]))) {
    stop("`which_summary` should be on of ", 
         paste(colnames(sums[[which_parameter]]), collapse=", "))
  }
  if(main =="automatic") main <- paste("Predicted", which_summary, "of", which_parameter)
  plotPointillistPainting(
    locs = prediction_block$locs, 
    field = prediction_block$summaries[[which_parameter]][,which_summary],
    main = main
  )
}

#' Plot the estimation of a GeoNonStat object
#'
#' @param object an object of class \code{GeoNonStat}
#' @param estimation an output of the estimation function on object. 
#' If time consuming, it is adviced to pre-compute estimation and pass it as 
#' an argument. Default to NULL otherwise (recomputed before plot)
#' @param which_parameter parameter to plot, default to NULL.
#' Should be one of `names(estimation)`
#' @param which_summary summary to plot, default to "mean". On of the available 
#' statistics of estimation object : "mean", "sd", "median", "q 2.5%" or "q 97.5%".
#' @param main character, title of the plot. 
#' @family {Visualisations}
#' @returns a plot
#' @export
#'
#' @examples
#' plotGeoNonStatEstimation(object = processedGnsDemo, 
#'                          which_parameter = "field", which_summary = "mean")
plotGeoNonStatEstimation <- function(
    object, 
    estimation = NULL, 
    which_parameter = NULL, 
    which_summary = "mean", 
    main = "automatic"){
  if(is.null(estimation)) estimation <- estimate(object)
  est <- estimation$summaries$at_observed_locs
  if(is.null(which_parameter) || !(which_parameter %in% names(est))) {
    stop("`which_parameter` should be on of ", paste(names(est), collapse=", "))
  }
  if(is.null(which_summary) || !(which_summary %in% colnames(est[[which_parameter]]))) {
    stop("`which_summary` should be on of ", 
         paste(colnames(est[[which_parameter]]), collapse=", "))
  }
  if(main =="automatic") main <- paste("Estimated", which_summary, "of", which_parameter)
  plotPointillistPainting(
    locs = object$vecchia_approx$observed_locs, 
    field = est[[which_parameter]][,which_summary], 
    main = main
  )
}