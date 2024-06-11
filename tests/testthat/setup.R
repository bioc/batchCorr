batchB <- getBatch(peakTable = PTfill, meta = meta, 
                   batch = meta$batch, select = 'B')
batchF <- getBatch(peakTable = PTfill, meta = meta, 
                   batch = meta$batch, select = 'F')
# # Drift correction using QCs
BCorr <- correctDrift(peakTable = batchB$peakTable, 
                      injections = batchB$meta$inj, 
                      sampleGroups = batchB$meta$grp, QCID = 'QC', 
                      G = seq(5,35,by=3), modelNames = c('VVE', 'VEE'),
                      report = FALSE)
# More unbiased drift correction using QCs & external reference samples
FCorr <- correctDrift(peakTable = batchF$peakTable, 
                      injections = batchF$meta$inj,
                      sampleGroups = batchF$meta$grp, QCID = 'QC',
                      RefID='Ref', G = seq(5,35,by=3), 
                      modelNames = c('VVE', 'VEE'),
                      report = FALSE)

# correctDrift low-level walkthrough
meta_test <- data.frame(injections = batchB$meta$inj,sampleGroups = batchB$meta$grp)
batchQC <- .getGroup(peakTable=batchB$peakTable, meta=meta_test, 
                     sampleGroup=batchB$meta$grp, select = "QC")
QCObject <- makeQCObject(peakTable = batchQC$peakTable, 
                         inj = batchQC$meta$inj)
clustered <- clust(QCObject$inj, QCObject$Feats, 
                   modelNames = c("VVE", "VEE"),
                   G = seq(5,35,by=10), report = FALSE)               
calculated_loess <- driftCalc(clustered, smoothFunc = "loess", 
                              report= FALSE)
calculated_spline <- driftCalc(clustered, smoothFunc = "spline", 
                               report= FALSE)
BatchObject = makeBatchObject(peakTable = batchB$peakTable, 
                              inj = batchB$meta$inj, QCObject = QCObject)  
corrected_none <- driftCorr(QCDriftCalc = calculated_spline, 
                            CorrObj = BatchObject, report = FALSE)
cleaned_none <- cleaned_none <- cleanVar(corrected_none, report = FALSE)

