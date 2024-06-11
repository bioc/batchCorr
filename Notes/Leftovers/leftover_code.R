# Make data with 1000 features
data("ThreeBatchData") # old data

PTnofill <- PTnofill[, 1:5000]
PTfill <- PTfill[, 1:5000]
meta <- meta

save(PTnofill, PTfill, meta, file = "ThreeBatchData.rda")

# Make TwoBatchData
data("ThreeBatchData") # old data

PTfill <- PTfill[meta$batch != "F", ]
PTnofill <- PTnofill[meta$batch != "F", ]
meta <- meta[meta$batch != "F", ]

save(PTnofill, PTfill, meta, file = "TwoBatchData.rda")

# Test data
test_that("Data reading works", {
  data("ThreeBatchData")
  expect_true(is(meta, "data.frame"))
  expect_true(is(PTnofill, "matrix"))
  expect_true(is(PTfill, "matrix"))
})

