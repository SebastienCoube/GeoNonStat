#' Create visualisations for a PP object, according to the vecchia_approx that produced it.
#' @param x an object of class PP, create with `createPP`
#' @param ... additional arguments (unused)
#' @param mar_var_loss boolean, default to TRUE. Should loss of margine variance be computed and plotted ?
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
  def.par <- par(no.readonly = TRUE)
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
    par(mfrow = c(2,2))
    plotPointillistPainting(
      x$vecchia_locs,  legend = F, main = "A PP sample",
      Matrix::solve(x$sparse_chol, c(rnorm(x$n_knots), rep(0, nrow(x$vecchia_locs))))[-seq(x$n_knots)]
    )
    plotPointillistPainting(
      x$vecchia_locs,  legend = F, main = "Another PP sample",
      Matrix::solve(x$sparse_chol, c(rnorm(x$n_knots), rep(0, nrow(x$vecchia_locs))))[-seq(x$n_knots)]
    )
    plotPointillistPainting(
      x$vecchia_locs,  legend = F, main = "Yet another PP sample",
      Matrix::solve(x$sparse_chol, c(rnorm(x$n_knots), rep(0, nrow(x$vecchia_locs))))[-seq(x$n_knots)]
    )
    plotPointillistPainting(
      x$vecchia_locs,  legend = F, main = "A last PP sample",
      Matrix::solve(x$sparse_chol, c(rnorm(x$n_knots), rep(0, nrow(x$vecchia_locs))))[-seq(x$n_knots)]
    )
  }
  par(def.par)
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
  def.par <- par(no.readonly = TRUE)
  omar <- def.par[["mar"]]
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
    points(x$knots,
      pch = 10,
      cex = 1,
      col = 1
    )
    legend(
      x = "topright",
      legend = "knots",
      pch = 10,
      col = 1,
      inset = c(-mylegend$rect$w, 0),
      bty = "n"
    )
  }

  #### Restaure margins
  par(mar = omar, mgp = def.par[["mgp"]])
}

#' Plots range ellipses for nonstationary covariance functions. Default is 0.1 correlation.
#' @param locs ellipses centers coordinates. A matrix with 2 columns, 1 row for each ellipse
#' @param log_range log_range at the ellipse centers, can have 1 or 3 columns. A matrix with 1 row (? TODO) and 1 or 3 columns.
#' @param add logical, default to FALSE. Add on existing plot ?
#' @param shrink numeric, shrinks or inflates the ellipses. shrink = 1 gives the
#' Mahalanobis distance = 1. Shrink = sqrt(8*nu) gives the ellipses corresponding to
#' correlation  = .1 (rho following INLA's terminology)
#' @param ... additional graphical parameters, sent to `plot` if `add=FALSE`
#'
#' @returns a plot object
#' @export
#'
#' @examples
#' locs <- matrix(rnorm(20), ncol = 2)
#' log_range <- matrix(rnorm(30), ncol = 3)
#' plotEllipses(locs, log_range, shrink = 0.1)
#'
#' locs <- matrix(rnorm(20), ncol = 2)
#' log_range <- matrix(rnorm(10), ncol = 1)
#' plotEllipses(locs, log_range, shrink = 0.1)
plotEllipses <- function(locs,
                         log_range,
                         shrink = .1,
                         add = FALSE,
                         ...) {
  if (ncol(log_range) == 3) {
    # to match parametrization in compute sparse chol
    log_range <- log_range %*% matrix(c(1 / sqrt(2), 1 / sqrt(2), 0, 1 /
      sqrt(2), -1 / sqrt(2), 0, 0, 0, 1), 3) * sqrt(2)
    matrices <- lapply(split(log_range, row(log_range)), expmat)
  }
  if (ncol(log_range) == 1) {
    matrices <- lapply(log_range, function(x) {
      diag(exp(x), 2)
    })
  }
  if (!add) {
    plot(locs, type = "n", ...)
  }
  for (i in seq_len(nrow(locs)))
  {
    # 2.447747 must be some bivariate confidence interval
    # shrink = 1 gives the package's Mahalanobis distance
    matrices[[i]] <- eigen(matrices[[i]])
    matrices[[i]] <- (matrices[[i]])$vec %*% diag(matrices[[i]]$val^2) %*% t(matrices[[i]]$vec) /
      (2.447747)^2
    ell <- ellipse::ellipse(matrices[[i]]) * shrink
    ell[, 1] <- ell[, 1] + locs[i, 1]
    ell[, 2] <- ell[, 2] + locs[i, 2]
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


#' Title getColors TODO
#'
#' @param x a vector
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
  col1 <- "#FFFFB2"
  col2 <- "#B10026"
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
#'
#' @returns a plot component (a `plot` if `add == FALSE`, `points` if `add ==TRUE`)
#' @export
#'
#' @examples
#' locs <- matrix(rnorm(2000), ncol = 2)
#' plotPointillistPainting(locs = locs, field = rnorm(1000))
plotPointillistPainting <- function(locs,
                                    field,
                                    add = FALSE,
                                    alpha = TRUE,
                                    legend = TRUE,
                                    ...) {
  args <- list(...)
  if (is.null(args$pch)) args[["pch"]] <- "."
  if (is.null(args$cex)) args[["cex"]] <- 4
  if (is.null(args$xlab)) args[["xlab"]] <- NA
  if (is.null(args$ylab)) args[["ylab"]] <- NA
  args$col <- getColors(field, alpha = TRUE)

# gestion de la mise en page
  if (legend && !add) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))

    layout(
      matrix(c(1, 2), nrow = 1),
      widths = c(6, 1)
    )
  }

  if (add) {
    do.call("points", c(list(x = locs), args))
  } else {
    par(mar = c(2, 2, 4, 0))
    do.call("plot", c(list(x = locs), args))
  }

# légende
  if (legend && !add) {
    pointillistColorscale(field)
  }
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
pointillistColorscale <- function(field, scalename = "") {
  origmar <- par("mar")
  origlwd <- par("lwd")
  par(mar = c(2, 0, 4, 0))
  par(lwd = 0.5)
  field <- field[!is.na(field)]
  if (length(field) == 0) {
    stop("only NA values in field")
  }
  barplot(
    rep(0.3, 50),
    col = getColors(1:50),
    width = rep(1, 50),
    space = 0,
    # lwd = 0.1,
    xlab = scalename,
    ylab = "",
    main = "",
    border = FALSE,
    ylim = c(0, 50),
    horiz = TRUE,
    xlim = c(-0.2, 0.1),
    axes = FALSE,
    mgp = c(0, 0, 0),
  )
  pos <- seq(0, 1, by = 0.25)
  values <- stats::quantile(field, probs = pos)
  text(
    y = 1 + pos * (50 - 1),
    x = -0.01,
    adj = c(0.95, 1),
    # format(signif(values), digits=4),
    round(values, 3),
    cex = 0.8
  )
  par(mar = origmar, lwd = origlwd)
}
