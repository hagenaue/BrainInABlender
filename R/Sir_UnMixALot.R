#' Sir_UnMixALot
#'
#' This function uses Microarray or RNAseq data derived from heterogeneous brain cortical samples to estimate the relative balance of each cell type across samples. To do this, we use a database(CellTypeSpecificGenes_Master3) of genes that have been previously-indicated to have cell type specific expression in either the forebrain or cortex of humans or mice. 

#' @param userInput A dataframe containing the gene expression data for the samples (RNAseq or microarray data that has already been variance stabilized and received appropriate quality control to remove outlier samples and large-scale technical artifacts), including one column of gene symbols.  Defaults to Error.
#' @param dataColumns A vector containing the column numbers for the gene expression data in the dataframe. 
#' @param geneColumn A single integer indicating the column number for the column that includes the gene symbols. 
#' @param species A single character string indicating the species from which the data were derived, currently allows the values  "mouse", "Mouse", "human", or "Human". 
#' @return A list containing two data frames: PublicationSpecific_CellTypeIndex and AveragePrimary_CellTypeIndex. These data frames provide estimates for the relative balance of each cell type across samples. The first data frame provides estimates based on the cell type specific gene lists provided by particular publications ("cell type indices"), the second data frame averages each of the publication-specific cell type indices to create an average cell type index for each primary cell type.
#' @importFrom plyr join
#' @keywords brain cortex cell type matrix deconvolution microarray RNAseq
#' @export
#' @examples 
#' Sir_UnMixALot(userInput=ZhangRNAseqData, dataColumns=c(2:18), geneColumn=1, species="mouse")


Sir_UnMixALot<- function(userInput, dataColumns, geneColumn, species){
  
  #required packages
  #install.packages("plyr")
  #library(plyr)
  
  
  #DATASET READ IN
  #userInput <- read.table(file = fileName, header = T, sep = ",", stringsAsFactors = F)
  #userInput <- read.table(file = "JoinedZhang_AsNumMatrix_Log2.csv", header = T, sep = ",", stringsAsFactors = F)
  
  
  #MH: I'm adding some code to double check that the parameters fit requirements:
  #Alek - this would be a good place to stop the code if there is an error.
  
  print("I like big brains and I cannot lie...")
  
  if(is.data.frame(userInput)==F){print("Error: Your userInput is not a data.frame")}else{print("Good: The userInput is a data frame")}
  
  if(is.integer(dataColumns)==F){print("Warning: Your dataColumns are not integer values")}else{print("Good: The dataColumns are integer values")}
  
  if(is.integer(geneColumn)==F){print("Warning: Your geneColumn is not an integer value")}else{print("Good: The geneColumn an integer value")}
  #MH: perhaps I should just check to make sure it is numeric?
  
  if(length(geneColumn)>1){print("Error: You have specified more than one geneColumn")}else{print("Good: You have specified only one geneColumn")}
  
  if(species%in%c("Mouse", "mouse", "Human", "human")){print("Good: Your species selection is one of our included species.")}else{print("Error: You have specified a different species than the two included in our database (mouse and human).")}
  
  
  #MH: I generalized this to specify the gene column number:
  GeneNamesForJoinedInput <- userInput[,geneColumn]
  GeneNamesForJoinedInput <- as.matrix(GeneNamesForJoinedInput)
  
  #MH: I generalized this to specify the data column numbers:
  TempJoinedInput_AsNum <- userInput[,dataColumns]
  
  #MH: I added a little feedback to the user here:
  print("The dimensions of the matrix for the user's gene expression data:")
  print(dim(TempJoinedInput_AsNum))
  
  ##################################################
  
  #correlation matrices
  #This code is all commented out - should we just delete it? It is used to detect large trends or outliers in the dataset.
  
  #temp<-cor(TempJoinedInput_AsNum)
  #row.names(temp)<-colnames(userInput)
  
  #png("09 Sample Sample Correlations Heatmap.png")
  #  heatmap(temp, main="Visualizing correlations between entire samples", xlab="Red=Less correlated, Light yellow=Highly correlated")
  #dev.off()
  #Note that the heatmap can be tailored to focus on a certain level of correlation by using the command zlim=c(lower limit, upper limit)
  
  #Visualize the sample-sample correlations using a boxplot:
  #png("09 Boxplot Sample Sample Correlations.png", width=2000, height=300)
  #  boxplot(data.frame(cor(TempJoinedInput_AsNum)), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of sample-sample correlations", xlab="Subject", ylab="Sample-Sample Correlations")
  #  Median10thQuantile<-median(apply((cor(TempJoinedInput_AsNum)), 1, quantile, 0.1))
  #  MedianQuantile<-median(apply((cor(TempJoinedInput_AsNum)), 1, quantile, 0.5))
  #  abline(a=Median10thQuantile, b=0, col=2)
  #  abline(a=MedianQuantile, b=0, col=3)
  #  mtext(paste("Median Sample-Sample Correlation=", round(MedianQuantile, digits=3), sep=" ")) 
  #dev.off()
  
  
  
  ##################################################
  
  #MH: We should probably double check whether there are NA columns/rows too - they may crash the function, since the z-score function depends on mean and sd calculations. 
  #Actually, it may actually be o.k. - we should probably test it though.
  
  
  
  ##################################################
  
  #MH: This code determines how many of the rows (probes/genes) have zero variability and therefore have to be removed from the   calculations - I've made it so that the user gets some feedback about what is happening:
  JoinedInput_StDev<-apply(TempJoinedInput_AsNum, 1, function(y) sd(y, na.rm=T)) 
  
  print("The number of rows in the dataset that show no variability (sd=0):")
  #MH: this originally said "return" but I'm not sure why. I changed it to print
  print(sum(JoinedInput_StDev==0))
  #5314
  print("The percentage of rows in the dataset that show no variability (sd=0):")
  print((sum(JoinedInput_StDev==0)/length(JoinedInput_StDev))*100)
  print("These rows will be removed.")
  
  JoinedInput_AsNumMatrix_Log2_NoSD0<-TempJoinedInput_AsNum[JoinedInput_StDev>0,]
  temp<-GeneNamesForJoinedInput
  GeneNamesForJoinedInput_NoSD0<-temp[JoinedInput_StDev>0]
  
  GeneNamesForJoinedInput_NoSD0 <-as.matrix(GeneNamesForJoinedInput_NoSD0)
  
  ZscoreInput<-t(scale(t(JoinedInput_AsNumMatrix_Log2_NoSD0), center=T, scale=T))#Zscores the data 
  write.csv(ZscoreInput, "ZscoreInput.csv")
  print("Output added to working directory: ZscoreInput.csv")
  print("...ZscoreInput.csv is a data frame that includes a a z-scored version of your input (centered and scaled by row), with all rows removed that had zero variability.")
  
  sum(is.na(ZscoreInput))
  # ZERO WOOOOOOOOOO
  
  #############################################################
  #Code for determining if the gene symbols are unique, and if not averaging by gene symbol.
  
  #MH: I added a little feedback to the user here:
  print("The number of unique gene symbols (rows) included in the user's input after filtering out genes with expression that completely lacked variability (sd=0):")
  print(length(unique(GeneNamesForJoinedInput_NoSD0)))
  
  if(length(unique(GeneNamesForJoinedInput_NoSD0))==length((GeneNamesForJoinedInput_NoSD0[,1]))){print("All gene symbols are now unique.")}else{
    print("The gene symbols are not unique, the z-scored data will be averaged by gene symbol and then re-z-scored")
    
    temp<-matrix(0, length(names(table(GeneNamesForJoinedInput_NoSD0[,1]))), ncol(ZscoreInput))
    
    for(i in c(1:ncol(ZscoreInput))){
      temp[,i]<-tapply(ZscoreInput[,i], GeneNamesForJoinedInput_NoSD0[,1], mean)
    }
    row.names(temp)<-names(table(GeneNamesForJoinedInput_NoSD0))
    colnames(temp)<-colnames(ZscoreInput)
    ZscoreInput<-t(scale(t(temp), center=T, scale=T))
    GeneNamesForJoinedInput_NoSD0<-as.matrix(as.character(row.names(temp)))
    
    write.csv(ZscoreInput, "ZscoreInput_AveragedByGeneSymbol.csv")
    print("Output added to working directory: ZscoreInput_AveragedByGeneSymbol.csv")
    print("...ZscoreInput_AveragedByGeneSymbol.csv is a data frame that includes a a z-scored version of your input, with all rows removed that had zero variability. Data was then averaged by gene symbol, and z-scored again so that the data for each gene symbol is centered and scaled.")
  }
  

  
  #############################################################
  
  #MH: We will probably make the cell type gene database part of the R package, so eventually reading in this data.frame won't be necessary.
  
  #I think loading in CellTypeSpecificGenesMaster3 isn't necessary now that I placed it in the R package directory.
  #CellTypeSpecificGenes_Master3<-read.csv("CellTypeSpecificGenes_Master3.csv", header=T)
  colnames(CellTypeSpecificGenes_Master3)[4]<-"GeneSymbol_Human"
  colnames(CellTypeSpecificGenes_Master3)[5]<-"GeneSymbol_Mouse"
  
  #I'm testing out whether removing the Doyle gene lists improves the signatures.
  CellTypeSpecificGenes_Master3<-CellTypeSpecificGenes_Master3[CellTypeSpecificGenes_Master3$Citation!="Doyle_Cell_2008", ]
  CellTypeSpecificGenes_Master3<-droplevels(CellTypeSpecificGenes_Master3)
  #MH - the removing NA values code was moved later.
  
  #############################################################
  
  dir.create("Detailed_CellTypeAnalysisOutput")
  setwd("Detailed_CellTypeAnalysisOutput")
  
  #joining celltype to zscore data
  tempForJoin <- data.frame(GeneNamesForJoinedInput_NoSD0, ZscoreInput, stringsAsFactors=F)
  
  #########
  #choosing which gene symbol to reference in the cell type specific gene database based on user Input "species":
  #MH: this code wasn't working originally, I think because it wanted the else statement to immediately follow the {} (no new line)
  if (species == "Mouse" || species == "mouse"){
    CellTypeSpecificGenes_Master3NoNA <- CellTypeSpecificGenes_Master3[is.na(CellTypeSpecificGenes_Master3$GeneSymbol_Mouse)==F,]
    colnames(CellTypeSpecificGenes_Master3NoNA)[5] <- "GeneNamesForJoinedInput_NoSD0"
    ZscoreInput_Expression_CellType<-join(CellTypeSpecificGenes_Master3NoNA, tempForJoin, by="GeneNamesForJoinedInput_NoSD0", type="inner")
    write.csv(ZscoreInput_Expression_CellType, "ZscoreInput_Expression_CellType.csv")
  } else if(species == "Human" || species == "human"){
    CellTypeSpecificGenes_Master3NoNA <- CellTypeSpecificGenes_Master3[is.na(CellTypeSpecificGenes_Master3$GeneSymbol_Human)==F,]
    colnames(CellTypeSpecificGenes_Master3NoNA)[4] <- "GeneNamesForJoinedInput_NoSD0"
    ZscoreInput_Expression_CellType<-join(CellTypeSpecificGenes_Master3NoNA, tempForJoin, by="GeneNamesForJoinedInput_NoSD0", type="inner")
    write.csv(ZscoreInput_Expression_CellType, "ZscoreInput_Expression_CellType.csv")
  }else{print("Error: The designated species is not in our database.")} 
  
  print("Output added to working directory: ZscoreInput_Expression_CellType.csv")
  print("...ZscoreInput_Expression_CellType.csv is a data frame listing all of the cell type specific genes in your dataset and their respective expression for each sample (in z-scores).") 
  print("...Please note that the cell type specific genes included in this output have not yet been filtered based on whether they were identified as specifically expressed in different cell types in different publications (e.g. a gene that has been identified as specifically expressed in astrocytes in one publication but identified as specifically expressed in neurons in another publication).")
  
  ####################################
  #///////////////////////////////////
  ####################################
  
  #do these need to get output or are they testing?
  #MH: I commented it out because it was an earlier version of the analysis that doesn't work as well: 
  
  #CELLTYPE PRIMARY BY SAMPLE
  # AVE_Expression_CellType_Primary_bySample<-matrix(NA, nrow=length(names(table(ZscoreInput_Expression_CellType$CellType_Primary))), ncol=(ncol(ZscoreInput_Expression_CellType)-14))
  # row.names(AVE_Expression_CellType_Primary_bySample)<-names(table(ZscoreInput_Expression_CellType$CellType_Primary))
  # colnames(AVE_Expression_CellType_Primary_bySample)<-colnames(tempForJoin)[-1]
  # 
  # for(i in c(15:ncol(ZscoreInput_Expression_CellType))){
  #   AVE_Expression_CellType_Primary_bySample[,(i-14)]<-tapply(ZscoreInput_Expression_CellType[,i], ZscoreInput_Expression_CellType$CellType_Primary, function(y) mean(y, na.rm=T))
  # }
  # 
  # png("CorrMatrixCellTypeVsCellType_HeatMap.png")
  # heatmap(cor(t(AVE_Expression_CellType_Primary_bySample[,is.na(AVE_Expression_CellType_Primary_bySample[1,])==F])), margins = c(15, 15), cex.lab=0.5)
  # dev.off()
  # 
  # CorrelationMatrixCellTypeVsCellType<-cor(t(AVE_Expression_CellType_Primary_bySample[,is.na(AVE_Expression_CellType_Primary_bySample[1,])==F]))
  # 
  # write.csv(CorrelationMatrixCellTypeVsCellType, "CorrelationMatrixCellTypeVsCellType.csv")
  # #Huh - all correlations are positive. Perhaps because some samples simply have less reads or more artifacts? 
  
  #####################################
  # CELL TYPE TAG BY SAMPLE
  #Note: this version of the analysis does not remove genes identified as "specific" to more than one cell type yet.
  
  write.csv(table(ZscoreInput_Expression_CellType$Tag), "NumberOfGenesInYourDatasetFoundInEachPublicationsCellTypeSpecificGeneList.csv")
  
  print("A file has been outputted to your working directory NumberOfGenesInYourDatasetFoundInEachPublicationsCellTypeSpecificGeneList.csv.csv that tells you the number of gene symbols included in your dataset that were found in each publication's cell type specific gene list.")
  
  AVE_Expression_CellType_Tag_bySample<-matrix(NA, nrow=length(names(table(ZscoreInput_Expression_CellType$Tag))), ncol=ncol(ZscoreInput_Expression_CellType)-14)
  row.names(AVE_Expression_CellType_Tag_bySample)<-names(table(ZscoreInput_Expression_CellType$Tag))
  colnames(AVE_Expression_CellType_Tag_bySample)<-colnames(tempForJoin)[-1]
  
  for(i in c(15:ncol(ZscoreInput_Expression_CellType))){
    AVE_Expression_CellType_Tag_bySample[,(i-14)]<-tapply(ZscoreInput_Expression_CellType[,i], ZscoreInput_Expression_CellType$Tag, mean)
  }
  
  png("CorrMatrixCellIndexVsCellIndex_HeatMap.png", width=1000, height=1000)
  heatmap(cor(t(AVE_Expression_CellType_Tag_bySample[,is.na(AVE_Expression_CellType_Tag_bySample[1,])==F])), cex.lab=0.3, margins = c(20, 20))
  dev.off()
  
  CorrelationMatrixCellIndexVsCellIndex<-cor(t(AVE_Expression_CellType_Tag_bySample[,is.na(AVE_Expression_CellType_Tag_bySample[1,])==F]))
  
  write.csv(CorrelationMatrixCellIndexVsCellIndex, "CorrelationMatrixCellIndexVsCellIndex.csv")
  
  write.csv(AVE_Expression_CellType_Tag_bySample, "AVE_Expression_CellType_Tag_bySample.csv")
  
 print("Output added to working directory: AVE_Expression_CellType_Tag_bySample.csv, CorrelationMatrixCellIndexVsCellIndex.csv")
 print("...The file AVE_Expression_CellType_Tag_bySample.csv includes the average z-score for each sample for all genes identified as cell type specific for each cell type from each publication - a cell type index. This can be treated as one type of estimate of the relative balance of each cell type across all samples.")
print("...The file CorrelationMatrixCellIndexVsCellIndex.csv shows the correlation between each of these cell type indices.")
print("...Please note that the cell type specific genes included in this output have not yet been filtered based on whether they were identified as specifically expressed in different cell types in different publications (e.g. a gene that has been identified as specifically expressed in astrocytes in one publication but identified as specifically expressed in neurons in another publication).")
  
  
  #############################################
  #////////////////////////////////////////////
  #############################################
  
  #Alright, so part of the trouble here is that there hasn't been any removal of overlapping probes yet, and we aren't averaging by tag. Let's go ahead and do that.
  

  #Making a storage matrix to store information about overlap between primary indices:
  CellTypeSpecificGenes_Master3_Overlap<-matrix(0, length(table(ZscoreInput_Expression_CellType$CellType_Primary)), length(table(ZscoreInput_Expression_CellType$CellType_Primary)) )
  
  colnames(CellTypeSpecificGenes_Master3_Overlap)<-names(table(ZscoreInput_Expression_CellType$CellType_Primary))
  row.names(CellTypeSpecificGenes_Master3_Overlap)<-names(table(ZscoreInput_Expression_CellType$CellType_Primary))
  
  #Quantifying overlap between primary cell type indices:
  for(i in 1: length(table(ZscoreInput_Expression_CellType$CellType_Primary))){
    for(j in 1: length(table(ZscoreInput_Expression_CellType$CellType_Primary))){
      
      CellTypeSpecificGenes_Master3_Overlap[i,j]<-sum(ZscoreInput_Expression_CellType[ZscoreInput_Expression_CellType$CellType_Primary==names(table(ZscoreInput_Expression_CellType$CellType_Primary)[i]), 4]%in%ZscoreInput_Expression_CellType[ZscoreInput_Expression_CellType$CellType_Primary==names(table(ZscoreInput_Expression_CellType$CellType_Primary)[j]), 4])/length(ZscoreInput_Expression_CellType[ZscoreInput_Expression_CellType$CellType_Primary==names(table(ZscoreInput_Expression_CellType$CellType_Primary)[i]), 4])
      
    }
  }
  
  write.csv(CellTypeSpecificGenes_Master3_Overlap, "CellTypeSpecificGenes_Master3_Overlap.csv")
  
  print("Output added to working directory: CellTypeSpecificGenes_Master3_Overlap.csv")
  print("...CellTypeSpecificGenes_Master3_Overlap.csv indicates how many of the gene symbols included in your dataset are indicated to be cell type specific in cell type specific gene lists from different publications or for different cell types")
  
  
  #############################################
  #What happens if we eliminate overlap between primary categories and then make master indices:
  
  dim(ZscoreInput_Expression_CellType)
  # [1] 2914   31
  
  #Making an empty first row for the storage matrix:
  #
  ZscoreInput_Expression_CellType_NoPrimaryOverlap<-matrix(0, 1, (length(ZscoreInput_Expression_CellType[1,])))
  colnames(ZscoreInput_Expression_CellType_NoPrimaryOverlap)<-colnames(ZscoreInput_Expression_CellType)
  
  for(i in 1: length(table(ZscoreInput_Expression_CellType$CellType_Primary))){
    
    #Choosing all data for a particular primary cell type:
    TempCurrentIndexAllInfo<-ZscoreInput_Expression_CellType[ZscoreInput_Expression_CellType$CellType_Primary==names(table(ZscoreInput_Expression_CellType$CellType_Primary)[i]), ] 
    
    #All of the gene symbols within the current primary cell type:
    TempCurrentIndex<-ZscoreInput_Expression_CellType[ZscoreInput_Expression_CellType$CellType_Primary==names(table(ZscoreInput_Expression_CellType$CellType_Primary)[i]), 4] 
    
    #All of the gene symbols within all other primary cell types:
    TempAllOtherIndices<-ZscoreInput_Expression_CellType[ZscoreInput_Expression_CellType$CellType_Primary%in%names(table(ZscoreInput_Expression_CellType$CellType_Primary)[-i]), 4]
    
    #Grabs only rows of data with gene symbols not found in other primary cell type indices:
    ZscoreInput_Expression_CellType_NoPrimaryOverlap<-rbind(ZscoreInput_Expression_CellType_NoPrimaryOverlap, TempCurrentIndexAllInfo[(TempCurrentIndex%in%TempAllOtherIndices)==F,])
    
  }
  
  dim(ZscoreInput_Expression_CellType_NoPrimaryOverlap)
  # [1] 2435   31
  
  #removing that one dummy row:
  ZscoreInput_Expression_CellType_NoPrimaryOverlap<-ZscoreInput_Expression_CellType_NoPrimaryOverlap[-1,]
  
  dim(ZscoreInput_Expression_CellType_NoPrimaryOverlap)
  # [1] 2434   31
  
  write.csv(ZscoreInput_Expression_CellType_NoPrimaryOverlap, "ZscoreInput_Expression_CellType_NoPrimaryOverlap.csv")

  print("Output added to working directory: ZscoreInput_Expression_CellType_NoPrimaryOverlap.csv")
  print("...A data frame listing all of the cell type specific genes in your dataset and their respective expression for each sample (in z-scores).") 
  print("...The cell type specific genes included in this output have been filtered so that the data is now removed that was associated with gene symbols that were identified as specifically expressed in different cell types in different publications (e.g. a gene that has been identified as specifically expressed in astrocytes in one publication but identified as specifically expressed in neurons in another publication).")
  
  
  write.csv(table(ZscoreInput_Expression_CellType_NoPrimaryOverlap$Tag), "NumberOfGenesInYourDatasetFoundInEachPublicationsCellTypeSpecificGeneList_NoNonSpecific.csv")
  
  CellTypeSpecificGenes_Master3_PrimaryTable_NoPrimaryOverlap<-table(ZscoreInput_Expression_CellType_NoPrimaryOverlap$CellType_Primary)
  write.csv(CellTypeSpecificGenes_Master3_PrimaryTable_NoPrimaryOverlap, "NumberOfGenesInYourDatasetSpecificToEachPrimaryCellType.csv")
  
  print("Output added to working directory: NumberOfGenesInYourDatasetFoundInEachPublicationsCellTypeSpecificGeneList_NoNonSpecific.csv and  NumberOfGenesInYourDatasetSpecificToEachPrimaryCellType.csv")
print("...NumberOfGenesInYourDatasetFoundInEachPublicationsCellTypeSpecificGeneList_NoNonSpecific.csv that tells you the number of gene symbols included in your dataset that were found in each publication's cell type specific gene list. The cell type specific genes included in this output have been filtered so that the data is now removed that was associated with gene symbols that were identified as specifically expressed in different cell types in different publications (e.g. a gene that has been identified as specifically expressed in astrocytes in one publication but identified as specifically expressed in neurons in another publication).") 
print("...The file  NumberOfGenesInYourDatasetSpecificToEachPrimaryCellType.csv tells you how many genes included in your dataset were found to be specific to each primary cell type.")
  
    
  ######################################################################

  #This code just creates an object aligning each of the publication's cell type specific gene lists with their primary cell type. 
  temp<-data.frame(ZscoreInput_Expression_CellType_NoPrimaryOverlap$Tag, ZscoreInput_Expression_CellType_NoPrimaryOverlap$CellType_Primary) 
  dim(temp)
  # [1] 2434    2
  
  CellTypePrimaryVsTag<-unique(temp)
  dim(CellTypePrimaryVsTag)
  # [1] 38  2
  
  colnames(CellTypePrimaryVsTag)<-c("Tag","CellType_Primary")
  head(CellTypePrimaryVsTag)
  
  
  #####################################################
  
  #This code outputs the average z-score for each publication's cell type specific gene lists, after having filtered out the genes that are found to be "specific" to different kinds of cells in different publications (i.e., actually non-specific!):
  
  ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag<-matrix(0, nrow=length(table(ZscoreInput_Expression_CellType_NoPrimaryOverlap$Tag)), ncol=(length(ZscoreInput_Expression_CellType_NoPrimaryOverlap[1,])-14))
  
  row.names(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag)<-names(table(ZscoreInput_Expression_CellType_NoPrimaryOverlap$Tag))
  
  colnames(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag)<-colnames(ZscoreInput_Expression_CellType_NoPrimaryOverlap)[-c(1:14)]
  
  for(i in c(15:length(ZscoreInput_Expression_CellType_NoPrimaryOverlap[1,]))){
    ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag[,i-14]<-tapply(ZscoreInput_Expression_CellType_NoPrimaryOverlap[,i], ZscoreInput_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  head(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag)
  
  setwd("..")
  getwd()
  
  write.csv(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag, "Zscore_Expression_CellType_NoPrimaryOverlap_MeanTag.csv")
  
  print("Output added to working directory: ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag")
  print("... ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag is a data frame that includes the publication-specific cell type indices for each sample: the average z-score for each publication's cell type specific gene lists, after having filtered out the genes that are found to be specific to different kinds of cells in different publications (i.e., genes that have expression that is actually non-cell type specific!).")
  
  setwd("Detailed_CellTypeAnalysisOutput")
  
  #Making histograms for each cell type tag:
  dir.create("Cell Type Histograms")
  setwd("Cell Type Histograms")
  for(i in 1:nrow(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag)){
    png(paste("Histogram_", row.names(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag)[i], ".png", sep=""))
    hist(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag[i,], main=row.names(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag)[i], breaks=16, col=i)
    dev.off()
  }
  setwd("..")
  getwd()
  
  print("A subdirectory named Cell Type Histograms has been created within your working directory.")
print("...This directory contains histograms illustrating the distribution of each of your publication-specific cell type indices across samples. You probably want to examine these for extreme outliers.")

setwd("..")
getwd()
  
  ############################
  
  png("Heatmap_CellType_NoPrimaryOverlap_MeanTag.png", height=1000, width=1000)
  heatmap(cor(t(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag[,is.na(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag[1,])==F])), cex.lab=0.3, margins = c(20, 20))
  dev.off()
  
  CellType_NoPrimaryOverlap_MeanTag_CorrMatrix<-cor(t(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag[,is.na(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag[1,])==F]))
  
  head(CellType_NoPrimaryOverlap_MeanTag_CorrMatrix)
  
  write.csv(CellType_NoPrimaryOverlap_MeanTag_CorrMatrix, "CellType_NoPrimaryOverlap_MeanTag_CorrMatrix.csv")
  
  print("Output added to working directory: CellType_NoPrimaryOverlap_MeanTag_CorrMatrix.csv and Heatmap_CellType_NoPrimaryOverlap_MeanTag.png.")
print("...Both of these files illustrate the correlations between your publication-specific cell type indices.")
  
  ###############################################
  
  #This code outputs the average z-score across each of the publication-specific cell type indices:
  
  temp2<-tapply(ZscoreInput_Expression_CellType_NoPrimaryOverlap[,15], ZscoreInput_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  Tag<-names(temp2)
  
  CellTypePrimaryVsTag2<-join(as.data.frame(Tag), as.data.frame(CellTypePrimaryVsTag), by="Tag")
  
  ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean<-matrix(0, nrow=length(table(ZscoreInput_Expression_CellType_NoPrimaryOverlap$CellType_Primary)), ncol=(length(ZscoreInput_Expression_CellType_NoPrimaryOverlap[1,])-14))
  
  row.names(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean)<-names(table(ZscoreInput_Expression_CellType_NoPrimaryOverlap$CellType_Primary))
  
  colnames(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean)<-colnames(ZscoreInput_Expression_CellType_NoPrimaryOverlap)[-c(1:14)]
  
  for(i in c(1:length(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag[1,]))){
    ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean[,i]<-tapply(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag[,i], CellTypePrimaryVsTag2[,2], mean)
  }
  
  head(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean)
  
  write.csv(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean, "Zscore_Expression_CellType_NoPrimaryOverlap_Mean.csv")
  
  print("Output added to working directory: Zscore_Expression_CellType_NoPrimaryOverlap_Mean.csv.")
print("...This file includes the average cell type indices for each primary cell type for each sample: the average of the publication- specific cell type indices.")
  
  ##############################################

setwd("Detailed_CellTypeAnalysisOutput")

  #Making histograms for each primary cell type:
  dir.create("Primary cell type Histograms")
  setwd("Primary cell type Histograms")
  for(i in 1:nrow(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean)){
    png(paste("Histogram_", row.names(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean)[i], ".png", sep=""))
    hist(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean[i,], main=row.names(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean)[i], breaks=16, col=i)
    dev.off()
  }
  
print("A subdirectory named Primary cell type Histograms has been created within your working directory.")
print("...This directory contains histograms illustrating the distribution of each of the averaged cell type indices for each primary cell type across samples. You probably want to examine these for extreme outliers.")

  #MH: do we need to reset the working directory? I just added this code because I didn't see it anywhere else.
  setwd("..")
  getwd()
  
  setwd("..")
  getwd()
  
  #MH: Is this necessary?
  is.numeric(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean)
  
  png("Heatmap_CorMatrixPrimaryCellsNoOverlap.png", height=1000, width=1000)
  heatmap(cor(t(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean[,is.na(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean[1,])==F])), cex.lab=0.3, margins = c(20, 20))
  dev.off()
  
  CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix<-cor(t(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean[,is.na(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean[1,])==F]))
  
  head(CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix)
  
  write.csv(CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix, "CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix.csv")
  
  print("Output added to working directory: CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix.csv and Heatmap_CorMatrixPrimaryCellsNoOverlap.png.")
print("...Both of these files illustrate the correlations between the averaged cell type indices for each primary cell type across samples.")
  
  ##################################################
#Returning the final cell type indices to the workspace:
print("The output returned to your workspace is a list containing two data frames: PublicationSpecific_CellTypeIndex and AveragePrimary_CellTypeIndex. These dataframes provide estimates for the relative balance of each cell type across samples.")

  CellTypeIndexOutput<-list(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag, ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean)
  names(CellTypeIndexOutput)<-c("PublicationSpecific_CellTypeIndex", "AveragePrimary_CellTypeIndex")
  return(CellTypeIndexOutput)

}

