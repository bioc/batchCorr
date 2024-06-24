#' Extract m/z and rt from peak table
#'
#' Extract features from peak table and report their m/z and rt values
#'
#' @param PT a peak table with variables as columns
#' @param sep character separating mz from rt, e.g. "_"
#' @param start character from which to start the read of peakInfo (from PT colnames)
#' @param timepos Which position carries info about rt (1 for before separator; 2 for after separator)
#' @return a matrix with m/z and rt of features as columns
#'
#' @examples
#' data("ThreeBatchData")
#' # Extract peakinfo (i.e. m/z and rt of features). These column names have 2
#' # leading characters describing LC-MS mode -> start at 3
#' peakIn <- peakInfo(PT = PTnofill, sep = "@", start = 3)
#'
#' @export
peakInfo <- function(PT, sep = "_",
                     timepos = 2,
                     start = 1) {
    # Making sure separator is present in colnames of peak table (PT)
    if (!any(grepl(sep, colnames(PT)))) {
        stop("Separator not present in column names of peak table.")
    }
    # Making sure type of start is numeric in nature
    if (!is.numeric(start) || !is.numeric(timepos)) {
        whichNotNumeric <- c("'start'", "'timepos'")
        stop(
            "Arguments ",
            paste(whichNotNumeric[c(
                !is.numeric(start),
                !is.numeric(timepos)
            )], collapse = " & "),
            " not numeric."
        )
    }

    # Checking that separator is not longer than 1
    if (length(sep) > 1) {
        stop("Only use of one separator possible")
    }

    # Making matrix of mz and rt after splitting by separator
    peakInfo <- matrix(unlist(strsplit(colnames(PT), sep)), ncol = 2, byrow = TRUE)
    peakInfo[, 1] <- substr(peakInfo[, 1], start, max(nchar(peakInfo[, 1])))
    peakInfo <- matrix(as.numeric(peakInfo), ncol = 2)
    # Reverse column order if time is not after separator
    if (timepos != 2) {
        peakInfo <- peakInfo[, 2:1]
    }
    colnames(peakInfo) <- c("mz", "rt")
    rownames(peakInfo) <- paste("feature", seq_len(nrow(peakInfo)), sep = "_")
    return(peakInfo)
}
