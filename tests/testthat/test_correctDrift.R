test_that("Results from correctDrift has right dimensions",{
  # Check that all samples are preseved with and without reference samples
  expect_equal(nrow(batchB$meta), length(BCorr$TestInjs))
  expect_equal(nrow(batchF$meta), length(FCorr$TestInjs))
  # Check that supplying RefID results in object of longer length
  expect_lt(length(BCorr), length(FCorr))
})
                         
test_that("driftWrap extracts samples correctly", {
  expect_true(length(QCObject$inj) == nrow(QCObject$Feats))
})                         
                         
# test_that("clust throws warning when string is supplied to G", {
# expect_warning(clust(QCObject$inj, QCObject$Feats, 
#                       modelNames = c("VVV", "VEE"),
#                       G = c(1, "SOS"), report = FALSE))
# }) 
                      
test_that("clust throws error when invalid string is supplied to modelNames", {
expect_error(clust(QCObject$inj, QCObject$Feats, 
                   modelNames = c("SOS", "VEE"),
                   G = c(1, 5, 10), report = FALSE))
})


test_that("QC sample number and indexing is unchanged after clustering", {
    expect_identical(clustered$QCInjs, QCObject$inj)
})
                                 
test_that("Drift calculation results are different depending on smoothFunc", {
expect_false(isTRUE(all.equal(calculated_loess, calculated_spline)))
})

test_that("Bad string results error in drift calculation", {
expect_error(driftCalc(clustered, smoothFunc = "lel", report= FALSE))

})

test_that("QC sample number and indexing is same after drift calculation", {
expect_identical(calculated_spline$QCInjs, QCObject$inj)
})
                      
test_that("driftCorr and cleanVar work with RefID", {
  batchRef <- .getGroup(peakTable=batchB$peakTable, meta=meta_test, 
                        sampleGroup=batchB$meta$grp, select="Ref")
  RefObject <- makeBatchObject(peakTable = batchRef$peakTable, inj = batchRef$meta$inj, QCObject = QCObject)
  
  corrected_ref <- driftCorr(QCDriftCalc = calculated_spline, refList=RefObject,
                         refType = "one", CorrObj=BatchObject, report = FALSE)
  expect_lt(length(corrected_none), length(corrected_ref))
  
  cleaned_ref <- cleanVar(corrected_ref, report = FALSE)
  
  expect_lt(length(cleaned_none), length(cleaned_ref))
})
