# batchCorr
Within and between batch correction of LC-MS metabolomics data

## Description
This is package contains functions within three areas of batch correction. These algorithms were originally developed 
to increase quality and information content in data from LC-MS metabolomics. However, the algorithms should be applicable to 
other data structures/origins, where within and between batch irregularities occur.

The three areas indicated are:
- alignBatches(): Align features systematically misaligned between batches
- correctDrift(): Perform within-batch intensity drift correction
- normalizeBatches(): Perform between-batch normalization

### Batch alignment 
Batch alignment is achieved based on three concepts:
- Aggregation of feature presence/missingness on batch level.
- Identifying features with missingness within "the box", i.e. sufficiently similar in retention time and m/z.
- Ensuring orthogonal batch presence among feature alignment candidates.

### Drift correction
Drift correction is achieved based on:
- Clustering is performed on features in observation space (as opposed to the normally used observations in feature space)
- Clustering provides a tradeoff between 
  - modelling detail (multiple drift patterns within data set)
  - power per drift pattern
- Unbiased clustering is achieved using the Bayesian `mclust` R package

### Batch normalisation
Batch normalisation is achieved based on:
- QC/Reference (standard normalisation) or population (median normalisation)
- The choice between the two is based on a quality heuristic determining whether the QC/Ref samples are suitable for normalization. Otherwise population normalization is performed instead.

## Reference
The development and inner workings of these algorithms are reported in:
*Brunius C, Shi L, Landberg R, 2016. Large-scale untargeted LC-MS metabolomics data correction using between-batch feature alignment and cluster-based within-batch signal intensity drift correction. Metabolomics 12:173, doi: 10.1007/s11306-016-1124-4*

