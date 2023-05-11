#' @name outputResultsNOISeq
#' @param data
#' @param configFile
#'
#' @import NOISeq
#' @import dplyr
#'
#' @return myData
#'
#' @author Alberto Rodriguez-Izquierdo, 2021
#'


outputResultsNOISeq <- function(data, data_noiseq, dataPCA, outputPathName,configFile, treatment){
  
  dirOutput         <- paste0(configFile$output$outputDirNOISeq,outputPathName)
  
  graph             <- configFile$output$graph

  graphTypeExpr     <- configFile$output$graphType$expr

  graphTypeMD       <- configFile$output$graphType$MD

  graphTypeBiotypes <- configFile$output$graphType$biotypes

  graphTypePCA      <- configFile$output$graphType$PCA

  graphTypeq        <- configFile$output$graphType$q

  if (!dir.exists(dirOutput)){

    dir.create(dirOutput)

  }else{

    unlink(dirOutput, recursive = TRUE)

    dir.create(dirOutput)

  }
  
  write.table(data, file=paste0(dirOutput,'/', configFile$output$outputNameNOISeq), sep=';')

  if (isTRUE(graph)){

    dir.create(paste0(dirOutput,'/figures'))

    if (isTRUE(graphTypeExpr)){

      tiff(filename=paste0(dirOutput,'/figures/',outputPathName,'_expressionPlot.tiff'),units="in", width=5, height=5, res=300)

      DE.plot(data_noiseq, q = graphTypeq, graphic = "expr", log.scale = TRUE)

      dev.off()

    }

    if (isTRUE(graphTypeMD)){

      tiff(filename=paste0(dirOutput,'/figures/',outputPathName,'_MDPlot.tiff'),units="in", width=5, height=5, res=300)

      DE.plot(data_noiseq, q = graphTypeq, graphic = "MD")

      dev.off()

    }

    if (isTRUE(graphTypeBiotypes)){

      tiff(filename=paste0(dirOutput,'/figures/',outputPathName,'_BiotypesPlot.tiff'),units="in", width=5, height=5, res=300)

      DE.plot(data_noiseq, chromosomes = c(1, 2), q = graphTypeq, graphic = "distr")

      dev.off()
    }
    if(isTRUE(graphTypePCA)){

      myPCA <- dat(dataPCA, type='PCA')

      tiff(filename=paste0(dirOutput,'/figures/',outputPathName,'_PCA.tiff'),units="in", width=5, height=5, res=300)

      explo.plot(myPCA, factor=treatment)

      dev.off()

    }
  }


}


#' @name outputResultsNOISeq
#' @param data
#' @param configFile
#'
#' @import NOISeq
#' @import dplyr
#'
#' @return myData
#'
#' @author Alberto Rodriguez-Izquierdo, 2021
#'



outputResultsDESeq <- function(myResults, outputPathName,configFile, factors,myData){
  
  
  dirOutput         <- paste0(configFile$output$outputDirDESeq,outputPathName)
  
  if (!dir.exists(dirOutput)){
    
    dir.create(dirOutput)
    
  }else{
    
    unlink(dirOutput, recursive = TRUE)
    
    dir.create(dirOutput)
    
  }
  write.table(myResults, file=paste0(dirOutput,'/', configFile$output$outputNameDESeq), sep=';')
  
  dir.create(paste0(dirOutput,'/figures'))
    
  
  ####-----------------------DESeq heatMap--------------------------------#####
  
  ProcessingDESeqHeatMap(myData,configFile, dirOutput, outputPathName)
  
  #eval(parse(text=paste0('ProcessingDESeqPCA(myData,configFile$',factors,'$factor, dirOutput, outputPathName)')))
  
  
}


#' @name outputPathNameBuild
#' @param myfactor
#'
#' @import NOISeq
#' @import dplyr
#'
#' @return outputPathName
#'
#' @author Alberto Rodriguez-Izquierdo, 2021
#'


outputPathNameBuild <- function(myFactors,factors){
  
  if (length(unique(myFactors$strain))==1){
    
    outputSpecimen <- unique(myFactors$strain)
    
  }else{
    
    specvalues <-unique(myFactors$strain)
    
    outputSpecimen <- paste(specvalues[[1]],specvalues[[2]], sep = "_vs_")
    
  }
  if (length(unique(myFactors$factor))==1){
    
    outputfactor <- unique(myFactors$factor)
    
  }else{
    
    factorvalues <- unique(myFactors$factor)
    
    outputfactor <- paste(factorvalues[[1]],factorvalues[[2]], sep = "_vs_")
    
  }
  if (length(unique(myFactors$days))==1){
    
    outputdays <- unique(myFactors$days)
    
  }else{
    
    daysvalues <-unique(myFactors$days)
    
    outputdays <- paste(daysvalues[[1]],daysvalues[[2]], sep = "_vs_")
    
  }
  
  if (length(unique(myFactors$treatment))==1){
    
    outputtreatment <- unique(myFactors$treatment)
    
  }else{
    
    treatvalues <-unique(myFactors$treatment)
    
    outputtreatment <- paste(treatvalues[[1]],treatvalues[[2]], sep = "_vs_")
  }
  
  if (length(unique(myFactors$organ))==1){
    
    outputorgan <- unique(myFactors$organ)
    
  }else{
    
    organvalues <-unique(myFactors$organ)
    
    outputorgan <- paste(organvalues[[1]],organvalues[[2]], sep = "_vs_")
    
  }
  
  if (length(unique(myFactors$typeTreatment))==1){
    
    outputtypeTreatment <- unique(myFactors$typeTreatment)
    
  }else{
    
    typeTreatmentvalues <-unique(myFactors$typeTreatment)
    
    outputtypeTreatment <- paste(typeTreatmentvalues[[1]],typeTreatmentvalues[[2]], sep = "_vs_")
    
  }
  
  outputPathName <- paste0(outputSpecimen,'_',outputfactor,'_',outputdays,'_',outputtreatment,'_',outputorgan,'_',outputtypeTreatment,'_', factors)
  
  return(outputPathName)
  
}


#' @name outputCoincidences
#' @param data
#' @param configFile
#'
#' @import NOISeq
#' @import dplyr
#'
#' @return outputPathName
#'
#' @author Alberto Rodriguez-Izquierdo, 2021
#'


outputCoincidences <- function(myResults, outputPathName,configFile){
  
  dirOutput         <- paste0(configFile$output$outputDirCoincidences,outputPathName)
  
  if (!dir.exists(dirOutput)){
    
    dir.create(dirOutput)
    
  }else{
    
    unlink(dirOutput, recursive = TRUE)
    
    dir.create(dirOutput)
    
  }
  
  write.table(myResults, file=paste0(dirOutput,'/', configFile$output$outputNameCoincidences), sep=';', row.names = FALSE)
  
}