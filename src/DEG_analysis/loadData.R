#' @name loadData
#' @param dataPaht
#'
#'
#' @import dplyr
#'
#' @return dataCounts
#'
#' @author Alberto Rodriguez-Izquierdo, 2021

loadData <- function(configFile){

  dirsCounts <- list.dirs(path=configFile$dataPath, full.names=FALSE )

  for (nameSample in dirsCounts){

    if (!nameSample==""){

      mappingFile <- paste0(configFile$dataPath, nameSample, '/rescounts.csv')

      csv <- read.csv(mappingFile, sep=";", header=FALSE)
      df_csv <- data.frame(csv)
      clean_df_csv <- na.omit(df_csv)

      names(clean_df_csv)[names(clean_df_csv)=='V1'] <- 'gene_id'

      if (exists("dataCounts")==FALSE){

        dataCounts <- clean_df_csv
        names(dataCounts)[names(dataCounts)=='V2'] <- nameSample

      }else{

        dataCounts  <- cbind(dataCounts, V2 = clean_df_csv$V2)
        names(dataCounts)[names(dataCounts)=='V2'] <- nameSample
      }
    }
  }

  return(dataCounts)
}

#' @name creatingCountsByFactorsNoiseq
#' @param dataCounts
#' @param factor
#' @import dplyr
#'
#' @return selDataCounts
#'
#' @author Alberto Rodriguez-Izquierdo, 2021
#'

creatingCountsByFactorsNoiseq <- function(dataCounts,factor){
  
  SelDataCounts <- dplyr::select(dataCounts,factor$sample)
  
  row.names(SelDataCounts) <- dataCounts$gene_id

#  SelDataCounts <- cbind(gene_id=dataCounts$gene_id, SelDataCounts)

  return(SelDataCounts)
}


#' @name creatingCountsByFactorsDeseq
#' @param dataCounts
#' @param factor
#' @import dplyr
#'
#' @return selDataCounts
#'
#' @author Alberto Rodriguez-Izquierdo, 2021
#'

creatingCountsByFactorsDeseq <- function(dataCounts,factor){

  SelDataCounts <- dplyr::select(dataCounts,factor$sample)

  #SelDataCounts <- cbind(gene_id=dataCounts$gene_id, SelDataCounts)

  countData     <- as.matrix(SelDataCounts)

  rownames(countData) <- dataCounts$gene_id

  colnames(countData) <- factor$sample

  return(countData)
}




#' @name ReadNOISeqFactors
#' @param myCounts
#' @param myFactor
#' @param lengthGene
#' @import NOISeq
#' @import dplyr
#'
#' @return myData
#'
#' @author Alberto Rodriguez-Izquierdo, 2021
#'



ReadNOISeqFactors <- function(myCounts, lengthGene, myFactor){
  
  myData <- NOISeq::readData(data=myCounts, length = lengthGene, factors=myFactor)

  return(myData)
}

#' @name ReadDeseqFactors
#' @param myCounts
#' @param myFactor
#'
#' @import DESeq2
#' @import dplyr
#'
#' @return myData
#'
#' @author Alberto Rodriguez-Izquierdo, 2021
#'



ReadDeseqFactors <- function(countData,myFactor, configFile, factors){
  
  myFactor <- data.frame(myFactor)
  
  usedFactor <- eval(parse(text= paste0('configFile$', factors,'$factorDESeq')))
  
  #rownames(colData) <- 1
  
  if (length(unique(usedFactor)) > 1){
    
    for(columnFactor in usedFactor){
    
      designParameters <- columnFactor
      
      designParameters <- paste(designParameters, columnFactor, sep=' + ')
      
      }      
  
    colData <- myFactor
      
    myData    <- eval(parse(text=paste0('DESeqDataSetFromMatrix(countData=countData, colData=colData, design= ~ ',designParameters,')')))
      
    myData <- myData[rowSums(counts(myData)) > 5, ]
      
  
  }else{
    
    colData <- myFactor
    
    myData    <- eval(parse(text=paste0('DESeqDataSetFromMatrix(countData=countData, colData=colData, design= ~ ',usedFactor,')')))
    
    myData <- myData[rowSums(counts(myData)) > 5, ]
  
  }
  
    
#    if (!columnFactor == 'batch'){
#      if (!columnFactor == "sample"){

#        if (eval(parse(text=paste0('length(unique(myFactor$',columnFactor,'))')))==1){

#          eval(parse(text=paste0('myFactor$',columnFactor,' <- NULL')))

#        }

#      if (!exists('designParameters')){

#        designParameters <- columnFactor

#      }else{

#        designParameters <- paste(designParameters, columnFactor, sep=' + ')

#        }
#      }
#    }
  
  return(myData)
}


#' @name ReadNOISeqFactorsPCA
#' @param myCounts
#' @param myFactor
#' @param lengthGene
#' @import NOISeq
#'
#' @return myFactors
#' @author Alberto Rodriguez-Izquierdo, 2021


ReadNOISeqFactorsPCA <- function(myCounts, lengthGene, myFactor){
 
  myCounts2 <- myCounts

  myCounts2 <- myCounts2 + runif(nrow(myCounts2) * 4,3,5)

  mydata2   <- NOISeq::readData(data=myCounts2, length = lengthGene, factors=myFactor)


  return(mydata2)
}
