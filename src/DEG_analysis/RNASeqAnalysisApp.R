#' @name noiseqApp
#' @param root
#' @param datapath
#' @import dplyr
#' @author Alberto Rodriguez-Izquierdo, 2021

RNASeqAnalysisApp <- function(root){


#  library(logging)

#  logging::loginfo('####--------Starting app---------------####',logger = logs)

  configFile <- getConfigFile(root)

  dataCounts <- loadData(configFile)
########---------------Processing data------------------########
####---------------------Building factors-------------------####

  factorFile <- readFactor(configFile)

  bFactors <- configFile[grepl("bFactor",names(configFile))]

  for (factors in names(bFactors)){

    eval(parse(text=paste0('myFactors',factors,' <- loadFactors(factorFile,configFile$',factors,')')))

  }
####----------------Read data with NOISeq App--------------####

##---------------Calculating gene length---------------------##

  if (!file.exists(paste0(configFile$geneLength$geneLengthDir,configFile$geneLength$geneLengthFilename))){

    geneLength <- calculateGeneLength(configFile$geneLength$geneLengthAnnotation,
                                      configFile$geneLength$geneLengthExtAnnot,
                                      configFile$geneLength$geneLengthDir,
                                      configFile$geneLength$geneLengthOutputFilename)

  }else{
    geneLength <- read.csv(paste0(configFile$geneLength$geneLengthDir,configFile$geneLength$geneLengthFilename), sep= ';')
  }

####-------------Building data in NOISeq object------------####

  for (factors in names(bFactors)){

    if (isTRUE(configFile$analysis$noiseq)){
      
      ##---------------------Loading data from Factors------------------------##

      eval(parse(text=paste0('myCounts', factors,'noiseq <- creatingCountsByFactorsNoiseq(dataCounts,myFactors',factors,')')))

      eval(parse(text=paste0('myData',factors,' <- ReadNOISeqFactors(myCounts',factors,'noiseq,geneLength,myFactors',factors,')')))

      eval(parse(text=paste0('myData',factors,'_2 <- ReadNOISeqFactorsPCA(myCounts',factors,'noiseq,geneLength,myFactors',factors,')')))

      ##--------------Define variables to put into the analysis---------------##

#      eval(parse(text=paste0('mySpecimen1 <- configFile$',factors,'$specimen1')))

#      eval(parse(text=paste0('mySpecimen2 <- configFile$',factors,'$specimen2')))

#      eval(parse(text=paste0('myComparison <- configFile$',factors,'$comparison')))

      ##--------------Processing data and calculating DEGenes-----------------##

#      eval(parse(text=paste0(mySpecimen1,'_vs_', mySpecimen2,'_', myComparison,' <- ProcessingNOISeqbio(myData',factors,',factors,myFactors',factors,',configFile)')))

#      eval(parse(text=paste0(mySpecimen1,'_vs_', mySpecimen2,'_', myComparison,'.deg <- calculateDEGenes(',mySpecimen1,'_vs_', mySpecimen2,'_', myComparison,',configFile)')))

      ##----------------------------Output Results----------------------------##
      
      outputPathName <- eval(parse(text=paste0('outputPathNameBuild(myFactors',factors,',factors)')))
      
      eval(parse(text=paste0('resProcessingNOISeq <- ProcessingNOISeqbio(myData',factors,',factors,myFactors',factors,',configFile)')))
      
      resProcessingNOISeqDEG <- calculateDEGenes(resProcessingNOISeq,configFile)

      eval(parse(text=paste0('outputResultsNOISeq(resProcessingNOISeqDEG,resProcessingNOISeq,myData',factors,'_2, outputPathName, configFile, configFile$',factors,'$treatment)')))


    }

    if (isTRUE(configFile$analysis$deseq)){

      eval(parse(text=paste0('myCounts', factors,'deseq <- creatingCountsByFactorsDeseq(dataCounts,myFactors',factors,')')))

      eval(parse(text=paste0('myData', factors,' <- ReadDeseqFactors(myCounts',factors,'deseq,myFactors',factors,', configFile, factors)')))
      
      for(usedFactor in eval(parse(text=paste0('configFile$',factors,'$factorDESeq')))){
        
        #usedFactor <- eval(parse(text= paste0('configFile$', factors,'$factorDESeq')))
          
        factorsToDeseq <- eval(parse(text=paste0('unique(myFactors',factors,'$',usedFactor,')')))
          
        factorsToDeseq <- append(usedFactor,factorsToDeseq)
          
        eval(parse(text=paste0('deseqResults',factors,' <- processingDeseq(myData',factors,', configFile, factorsToDeseq)')))
          
        outputPathName <- eval(parse(text=paste0('outputPathNameBuild(myFactors',factors,',usedFactor)')))
          
        eval(parse(text=paste0('outputResultsDESeq(deseqResults',factors,',outputPathName,configFile, factors, myData',factors,')')))
          
      }
    }
    eval(parse(text=paste0('union_data_analysis <- coincidences(resProcessingNOISeqDEG,deseqResults',factors,', myCounts',factors,'noiseq, configFile)')))
    
    eval(parse(text=paste0('outputCoincidences(union_data_analysis, outputPathName, configFile)')))
  }

  print('App finished successfully!')
}



