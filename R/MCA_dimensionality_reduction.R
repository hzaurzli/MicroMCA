#' @title  MicroMCA
#
#' @description  Dimensionality reduction of high-dimensional microbiome data using MCA
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
#
#' @param X Dataframe, Rows represent different OTUs or microorganisms classification entry and columns represent samples.
#' @param nmcs In order to estimate number of right singular vectors, In this case nmcs + 1 indicates number of right singular vectors.
#' @param reduction.name Method of the reduction default set to 'MCA' for feature data and mca.
#'
#' @details The selected parameters are adapted from Griffiths & Steyvers (2004).
#' @import Rcpp
#' @import stringr
#' @import tictoc
#' @import pbapply
#' @import data.table
#' @importFrom matrixStats rowVars
#' @importFrom irlba irlba
#' @export run_mca_dim
#'
#' @examples
#'


run_mca_dim <- function(X, nmcs = 10, reduction.name = "MCA") {
  # preprocessing matrix
  # ----------------------------------------------------

  library(matrixStats)
  library(tictoc)
  library(stringr)
  library(vegan)


  if (is.null(X)) {
    stop("Input X is null")
  }

  X_Flattening <- as.data.frame(t(rrarefy(t(X), min(colSums(X)))))
  X_Flattening <- as.matrix(t(X_Flattening))
  X_Flattening <- X_Flattening[rowVars(X_Flattening) != 0, ]
  X_Flattening <- X_Flattening[str_length(rownames(X_Flattening)) > 0, ]
  X_Flattening <- X_Flattening[!duplicated(rownames(X_Flattening)), ]
  otusN <- colnames(X_Flattening)
  samplesN <- rownames(X_Flattening)
  tic()
  message("Computing Fuzzy Matrix")
  MCAPrepRes <- MCAStep1(X_Flattening)
  toc()
  message("Computing SVD")
  tic()
  SVD <- irlba::irlba(A = MCAPrepRes$Z, nv = nmcs + 1, nu = 1)[seq(3)]
  toc()
  message("Computing Coordinates")
  tic()
  MCA <- MCAStep2(Z = MCAPrepRes$Z, V = SVD$v[, -1], Dc = MCAPrepRes$Dc)
  component <- paste0(reduction.name, "_", seq(ncol(MCA$otusCoordinates)))
  colnames(MCA$otusCoordinates) <- component
  colnames(MCA$samplesCoordinates) <- component
  rownames(MCA$otusCoordinates) <- otusN
  rownames(MCA$samplesCoordinates) <- samplesN
  MCA$stdev <- SVD$d[-1]
  class(MCA) <- "MCA"
  toc()
  return(MCA)
}


#' @rdname calculate_importance
#' @export calculate_importance
#'
#' @param dat list, MCA redution result.
#' @param info dataframe, Sample metadata for group.
#' @param reduction Method of the reduction default set to 'MCA' for feature data and mca.
#' @param dims A vector of integers indicating which dimensions to use with reduction embedding and loading for distance calculation.
#'
calculate_importance <- function(dat, info, reduction = "mca", dims = c(1:50)) {

  if (is.null(dat)) {
    stop("Input dat is null")
  }

  if (is.null(info)) {
    stop("Input info is null")
  }

  otusEmb = dat$otusCoordinates
  # samplesEmb = dat$samplesCoordinates

  dat_tmp = cbind(info,dat$samplesCoordinates)
  type = unique(info[,2])
  for (i in type) {
    assign(i, colMeans(dat_tmp[which(dat_tmp[,2] == i),3:ncol(dat_tmp)]))
  }

  if (length(type) == 2) {
    samplesEmb = t(as.matrix(cbind(get(type[1]), get(type[2]))))
    row.names(samplesEmb) = type
  } else if (length(type) > 2) {
    tmp = cbind(get(type[1]), get(type[2]))
    for (i in 3:length(type)) {
      tmp = cbind(tmp, get(type[i]))
    }
    samplesEmb = t(as.matrix(tmp))
    row.names(samplesEmb) = type
  } else {
    stop('info format error!')
  }

  message("\ncalculating distance\n")
  # calculate distances
  Distance <- pairDist(otusEmb, samplesEmb)
  FeatureRanking <- DistSort(Distance)
  return(FeatureRanking)
}



#' Importance Calculation
#'
#' @description Quickly calculate the importance in different OTUs
#'
#' @param x a matrix
#' @param y a matrix
#'
#' @return A Distance Matrix
pairDist <- function(x, y) {
  z <- fastPDist(y, x)
  rownames(z) <- rownames(y)
  colnames(z) <- rownames(x)
  return(z)
}


#' Sort OTUs/samples distance Matrix
#'
#' @param distance distance matrix with samples at rows and OTUs at columns
#'
#' @return list of ranking of features by OTUs
#'
DistSort <- function(distance) {
  library(data.table)
  message("\ncreating ranking\n")
  GroupOTUDistance_DF <- as.data.table(distance)
  GroupOTURanking <-
    pbapply::pblapply(
      GroupOTUDistance_DF,
      FUN = function(x) {
        names(x) <- rownames(distance)
        x <- sort(x, method = "quick")
        return(x)
      }
    )
  return(GroupOTURanking)
}
