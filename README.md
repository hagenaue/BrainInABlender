This package was created using R version 3.3.1

Installation instructions: To install this package you must have devtools installed. Devtools allows the user to install packages directly from github. Once devtools is installed, install the BrainInABlender package from Github using the code Below. This will give you access to the Sir_UnMixALot function.

Installation Code:

Install.packages("devtools")

library("devtools")

install_github("hagenaue/BrainInABlender")

library("BrainInABlender")

Dependencies: This function makes use of the plyr package for its join function. 

Sir_UnMixALot uses Microarray or RNAseq data derived from heterogeneous brain cortical samples to estimate the relative balance of each cell type across samples. To do this, we use a database(CellTypeSpecificGenes_Master3) of genes that have been previously-indicated to have cell type specific expression in either the forebrain or cortex of humans or mice.

Input: A dataframe containing the gene expression data for the samples (RNAseq or microarray data that has already been variance stabilized and received appropriate quality control to remove outlier samples and large-scale technical artifacts), including one column of gene symbols. Defaults to Error. You must also have single character string indicating the species from which the data were derived, currently allows the values "mouse", "Mouse", "human", or "Human". Defaults to "mouse".

Output: A list containing two data frames: PublicationSpecific_CellTypeIndex and AveragePrimary_CellTypeIndex. These data frames provide estimates for the relative balance of each cell type across samples. The first data frame provides estimates based on the cell type specific gene lists provided by particular publications ("cell type indices"), the second data frame averages each of the publication-specific cell type indices to create an average cell type index for each primary cell type.

Example function call: Sir_UnMixALot(userInput=ZhangRNAseqData, dataColumns=c(2:18), geneColumn=1, species="mouse")


Full List of Citations included in the Database:
Cahoy JD, Emery B, Kaushal A, Foo LC, Zamanian JL, Christopherson KS, et al. A transcriptome database for astrocytes, neurons, and oligodendrocytes: a new resource for understanding brain development and function. J Neurosci Off J Soc Neurosci. 2008 Jan 2;28(1):264–78. 
Daneman R, Zhou L, Agalliu D, Cahoy JD, Kaushal A, Barres BA. The mouse blood-brain barrier transcriptome: a new resource for understanding the development and function of brain endothelial cells. PloS One. 2010;5(10):e13741. 

Darmanis S, Sloan SA, Zhang Y, Enge M, Caneda C, Shuer LM, et al. A survey of human brain transcriptome diversity at the single cell level. Proc Natl Acad Sci U S A. 2015 Jun 9;112(23):7285–90.

Doyle JP, Dougherty JD, Heiman M, Schmidt EF, Stevens TR, Ma G, et al. Application of a translational profiling approach for the comparative analysis of CNS cell types. Cell. 2008 Nov 14;135(4):749–62.

GeneCards®: The Human Gene Database:  http://www.genecards.org

Sugino K, Hempel CM, Miller MN, Hattox AM, Shapiro P, Wu C, et al. Molecular taxonomy of major neuronal classes in the adult mouse forebrain. Nat Neurosci. 2006 Jan;9(1):99–107. 

Zeisel A, Muñoz-Manchado AB, Codeluppi S, Lönnerberg P, La Manno G, Juréus A, et al. Brain structure. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. Science. 2015 Mar 6;347(6226):1138–42.

Zhang Y, Chen K, Sloan SA, Bennett ML, Scholze AR, O’Keeffe S, et al. An RNA-sequencing transcriptome and splicing database of glia, neurons, and vascular cells of the cerebral cortex. J Neurosci Off J Soc Neurosci. 2014 Sep 3;34(36):11929–47. 
