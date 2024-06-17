test_that("Results from correctDrift has right dimensions",{
  # Check that all samples are preseved with and without reference samples
  expect_equal(nrow(batchB$meta), length(BCorr$TestInjs))
  expect_equal(nrow(batchF$meta), length(FCorr$TestInjs))
  # Check that supplying RefID results in object of longer length
  expect_lt(length(BCorr), length(FCorr))
})

test_that("correctDrift throws error when supplied various values", {
  # When invalid string is supplied to modelNames"
  expect_error(clust(QCObject$inj, QCObject$Feats, 
                   modelNames = c("SOS", "VEE"),
                   G = c(1, 5, 10), report = FALSE))       
  # When bad string is supplied to smoothFunc
  expect_error(driftCalc(clustered, smoothFunc = "less", report= FALSE))
})
                         
test_that("driftWrap extracts samples correctly", {
  expect_true(length(QCObject$inj) == nrow(QCObject$Feats))
})
                                               
# Random note: Clustering is not deterministic despite setting seed
test_that("Output of clust is correct", {
  # Check that QC sample number and indexing is unchanged after clustering
  expect_identical(clustered$QCInjs, QCObject$inj)
})

test_that("Output of driftCalc is correct", {                                
  # Check that outcome is different depending on smoothFunc value
  expect_false(isTRUE(all.equal(calculated_loess, calculated_spline)))
  # Check that QC sample number and indexing is same after drift calculation
  expect_identical(calculated_spline$QCInjs, QCObject$inj)
})
                     
test_that("Output of driftCorr and cleanVar is correct", {
  batchRef <- .getGroup(peakTable=batchB$peakTable, meta=meta_test, 
                        sampleGroup=batchB$meta$grp, select="Ref")
  RefObject <- makeBatchObject(peakTable = batchRef$peakTable, inj = batchRef$meta$inj, QCObject = QCObject)
  corrected_ref <- driftCorr(QCDriftCalc = calculated_spline, refList=RefObject,
                             refType = "one", CorrObj=BatchObject, 
                             report = FALSE)
  # Check that driftCorr result list is longer when using reference samples
  expect_lt(length(corrected_none), length(corrected_ref))
  # Check that QC sample number and indexing is same after correction
  expect_identical(corrected_ref$QCInjs, QCObject$inj)
  # Check that reference sample number and indexing is same after correction
  expect_identical(corrected_ref$RefInjs, RefObject$inj)

  cleaned_ref <- cleanVar(corrected_ref, CVlimit = 0.3, report = FALSE)   
  # Check that cleanVar result list is longer when using reference samples
  expect_lt(length(cleaned_none), length(cleaned_ref))
  # Check that QC sample number and indexing is same after cleaning
  expect_identical(cleaned_ref$QCInjs, QCObject$inj)
  # Check that reference sample number and indexing is same after cleaning
  expect_identical(corrected_ref$RefInjs, RefObject$inj)
})
