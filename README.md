Installation instructions: To install this package you must have devtools installed. Devtools allows the user to install packages directly from github. Once devtools is installed, run the code: install_github("hagenaue/BrainInABlender"). This will give you access to the Sir_UnMixALot function.

Dependencies: This function makes use of the plyr package for its join function. You must also have single character string indicating the species from which the data were derived, currently allows the values "mouse", "Mouse", "human", or "Human". Defaults to "mouse".

Sir_UnMixALot uses Microarray or RNAseq data derived from heterogeneous brain cortical samples to estimate the relative balance of each cell type across samples. To do this, we use a database(CellTypeSpecificGenes_Master3) of genes that have been previously-indicated to have cell type specific expression in either the forebrain or cortex of humans or mice.

Input: A dataframe containing the gene expression data for the samples (RNAseq or microarray data that has already been variance stabilized and received appropriate quality control to remove outlier samples and large-scale technical artifacts), including one column of gene symbols. Defaults to Error.

Output: A list containing two data frames: PublicationSpecific_CellTypeIndex and AveragePrimary_CellTypeIndex. These data frames provide estimates for the relative balance of each cell type across samples. The first data frame provides estimates based on the cell type specific gene lists provided by particular publications ("cell type indices"), the second data frame averages each of the publication-specific cell type indices to create an average cell type index for each primary cell type.

Example function call: Sir_UnMixALot(userInput=ZhangRNAseqData, dataColumns=c(2:18), geneColumn=1, species="mouse")

