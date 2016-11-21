Sir_UnMixALot uses Microarray or RNAseq data derived from heterogeneous brain cortical samples to estimate the relative balance of each cell type across samples. To do this, we use a database(CellTypeSpecificGenes_Master3) of genes that have been previously-indicated to have cell type specific expression in either the forebrain or cortex of humans or mice.

Input: A dataframe containing the gene expression data for the samples (RNAseq or microarray data that has already been variance stabilized and received appropriate quality control to remove outlier samples and large-scale technical artifacts), including one column of gene symbols. Defaults to Error.

Output: A list containing two data frames: PublicationSpecific_CellTypeIndex and AveragePrimary_CellTypeIndex. These data frames provide estimates for the relative balance of each cell type across samples. The first data frame provides estimates based on the cell type specific gene lists provided by particular publications ("cell type indices"), the second data frame averages each of the publication-specific cell type indices to create an average cell type index for each primary cell type.

