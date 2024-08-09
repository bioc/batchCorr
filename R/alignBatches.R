#' Multi-batch alignment merging features artificially split between batches
#'
#' @param peakInfo matrix with mz and rt in columns 1:2 (see e.g. ?peakInfo)
#' @param PeakTabNoFill Multi-batch peak table including missing values
#' @param PeakTabFilled Multi-batch peak table without missing values
#' @param batches Vector (length=nrow(PeakTab)) of batch identifiers
#' @param sampleGroups Vector (length=nrow(PeakTab)) of sample type
#' (e.g. "sample", "QC", "Ref)
#' @param selectGroup Which sample type to base alignment on 
#' (e.g. "sample", "QC" or "Ref")
#' @param NAhard proportion of NAs within batch for feature to be 
#' considered missing
#' @param mzdiff Tolerance for difference in m/z
#' @param rtdiff Tolerance for difference in retention time
#' @param report Whether to export diagnostic plots into your work 
#' directory (defaults to TRUE)
#' @param reportPath directory path for report
#'
#' @return batchAlign Object
#' @return Aligned peak table available under $PTalign
#' @return Returns NULL if no alignment candidates can be found
#' @export
#'
#' @examples
#' \dontshow{
#' .old_wd <- setwd(tempdir())
#' }
#' data("ThreeBatchData")
#' # Extract peakinfo (i.e. m/z and rt of features). These column names have 2
#' # leading characters describing LC-MS mode -> start at 3
#' peakIn <- peakInfo(PT = PTnofill, sep = "@", start = 3)
#' # Perform multi-batch alignment
#' alignBat <- alignBatches(
#'     peakInfo = peakIn, PeakTabNoFill = PTnofill,
#'     PeakTabFilled = PTfill, batches = meta$batch,
#'     sampleGroups = meta$grp, selectGroup = "QC",
#'     reportPath = "drift_report/"
#' )
#' # Extract new peak table
#' PT <- alignBat$PTalign
#' \dontshow{
#' setwd(.old_wd)
#' }
alignBatches <- function(peakInfo,
                        PeakTabNoFill,
                        PeakTabFilled,
                        batches,
                        sampleGroups,
                        selectGroup = "QC",
                        NAhard = 0.8,
                        mzdiff = 0.002,
                        rtdiff = 15,
                        report = TRUE,
                        reportPath = NULL) {
    if (report & is.null(reportPath)) {
        stop("Argument 'reportPath' is missing")
    }
    if (anyNA(PeakTabFilled)) {
        stop("Missing values in PeakTabFilled")
    }
    if (!selectGroup %in% sampleGroups) {
        stop("Group ", selectGroup, " not in metadata.")
    }
    if (report & is.null(reportPath)) {
        stop("Argument 'reportPath' is missing")
    }
    if (report & !is.null(reportPath)) {
        if (!endsWith(reportPath, "/")) {
            message("Adding slash to path to allow proper folder structure")
            reportPath <- paste0(reportPath, "/")
        }
        if (!file.exists(reportPath)) {
            message("Creating folder ", reportPath)
            dir.create(reportPath, recursive = TRUE)
        }
    }
    bFlag <- batchFlag(
        PTnofill = PeakTabNoFill,
        batch = batches,
        sampleGroup = sampleGroups,
        peakInfo = peakInfo,
        NAhard = NAhard
    )
    aIQ <- alignIndex(
        batchflag = bFlag,
        grpType = selectGroup,
        mzdiff = mzdiff,
        rtdiff = rtdiff,
        report = report,
        reportPath = reportPath
    )
    if (is.null(aIQ)) {
        return(NULL)
    }
    if (report) {
        plotAlign(batchflag = bFlag, alignindex = aIQ, reportPath = reportPath)
    }

    bAlign <- batchAlign(
        batchflag = bFlag,
        alignindex = aIQ,
        peaktable_filled = PeakTabFilled,
        batch = batches
    )
}
