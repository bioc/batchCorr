#' ThreeBatchData
#'
#' The untargeted LC-MS metabolomics data was collected as per the associated
#' publication. In brief, the data was preprocessed using the \pkg{xcms}
#' package and biological samples were removed.
#'
#' @format A list with three objects encompassing three batches, 48 QC samples,
#' 42 long-term reference samples and 5000 features
#' \describe{
#'   \item{PTnofill}{matrix, peak table with missing values}
#'   \item{PTfill}{matrix, peak table without missing values}
#'   \item{meta}{data.frame, sample metadata including batch, group and
#'   injection order}
#' }
#'
#' @name ThreeBatchData
#' @docType data
#' @author Carl Brunius
#' @references
#' Carl Brunius, Lin Shi, Rikard Landberg
#' Large-scale untargeted LC-MS metabolomics data correction using
#' between-batch feature alignment and cluster-based within-batch signal
#' intensity drift correction.
#' Metabolomics, 12:173. \url{https://doi.org/10.3390/metabo10040135}
NULL

#' @name PTnofill
#' @rdname ThreeBatchData
NULL

#' @name PTfill
#' @rdname ThreeBatchData
NULL

#' @name meta
#' @rdname ThreeBatchData
NULL
