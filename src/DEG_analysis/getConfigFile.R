#' App structure
#' @name getConfigFile
#' @param root
#' @description App to start geoapp
#'
#' @import logging
#' @return configFile
#'
#' @author Alberto Rodriguez-Izquierdo, 2021

getConfigFile <- function(root){

  #loginfo('Reading configFile...')
  configFile <- readConfigFile(root)


  #####-------------------Validations----------------------######

  #loginfo('Validating configFile...')


  #------------------Validation nodes----------------------#

  configFile <- nodesValidation(configFile)


  ######--------------Validation nodes content --------------######

#  error <- FALSE

  #Validate generalParameters


  validateDataPath    <- validateCharacter(configFile$dataPath)

  configFile$dataPath <- validateDataPath

  validateSamplesFactorPath    <- validateCharacter(configFile$samplesFactorFilename)

  configFile$samplesFactorFilename <- validateSamplesFactorPath

  validateNOISeq <- validateCharacter(configFile$analysis$noiseq)

  configFile$analysis$noiseq <- validateNOISeq

  validateDESeq <- validateCharacter(configFile$analysis$deseq)

  configFile$analysis$deseq <- validateDESeq


  bFactors <- configFile[grepl("bFactor",names(configFile))]

  for (factors in names(bFactors)){

    validateNameSamples1  <- eval(parse(text=paste0('validateCharacter(configFile$', factors,'$name_samples_1)')))

    eval(parse(text=paste0('configFile$', factors,'$name_samples_1 <- validateNameSamples1')))

    validateNameSamples2 <- eval(parse(text=paste0('validateCharacter(configFile$', factors,'$name_samples_2)')))

    eval(parse(text=paste0('configFile$', factors,'$name_samples_2 <- validateNameSamples2')))

    validateFactorNOISeq <- eval(parse(text=paste0('validateCharacter(configFile$', factors,'$factorNOISeq)')))
    
    eval(parse(text=paste0('configFile$', factors,'$factorNOISeq <- validateFactorNOISeq')))
    
    validateFactorDESeq <- eval(parse(text=paste0('validateCharacter(configFile$', factors,'$factorDESeq)')))

    eval(parse(text=paste0('configFile$', factors,'$factorDESeq <- validateFactorDESeq')))

    }


####-------------Validate condAnalysisNOISeqbio-------------------####

  validatenclust <- validateNumber(configFile$condAnalysisNOISeqbio$nclust)

  configFile$condAnalysisNOISeqbio$nclust <- validatenclust

  validate_k <- validateNumber(configFile$condAnalysisNOISeqbio$k)

  configFile$condAnalysisNOISeqbio$k <- validate_k

  validatenorm <- validateCharacter(configFile$condAnalysisNOISeqbio$norm)

  configFile$condAnalysisNOISeqbio$norm <- validatenorm

  validatefilter <- validateNumber(configFile$condAnalysisNOISeqbio$filter)

  configFile$condAnalysisNOISeqbio$filter <- validatefilter

  validate_r <- validateNumber(configFile$condAnalysisNOISeqbio$r)

  configFile$condAnalysisNOISeqbio$r <- validate_r

  validaterandom.seed <- validateNumber(configFile$condAnalysisNOISeqbio$random.seed)

  configFile$condAnalysisNOISeqbio$random.seed <- validaterandom.seed

####-------------------Validate geneLength------------------------####

  validategeneLengthDir <- validateCharacter(configFile$geneLength$geneLengthDir)

  configFile$geneLength$geneLengthDir <- validategeneLengthDir

  validategeneLengthFilename <- validateCharacter(configFile$geneLength$geneLengthFilename)

  configFile$geneLength$geneLengthFilename <- validategeneLengthFilename

  validategeneLengthAnnotation <- validateCharacter(configFile$geneLength$geneLengthAnnotation)

  configFile$geneLength$geneLengtAnnotation <- validategeneLengthAnnotation

  validategeneLengthExtAnnot <- validateCharacter(configFile$geneLength$geneLengtExtAnnot)

  configFile$geneLength$geneLengtExtAnnot <- validategeneLengthExtAnnot

  validategeneLengthOutputFilename <- validateCharacter(configFile$geneLength$geneLengtOutputFilename)

  configFile$geneLength$geneLengtOutputFilename <- validategeneLengthOutputFilename

####-------------------Validate degenesParam------------------------####

  validatedegenesParamq <- validateNumber(configFile$degenesParam$q)

  configFile$degenesParam$q <- validatedegenesParamq

  validatedegenesParamm <- validateCharacter(configFile$degenesParam$m)

  configFile$degenesParam$m <- validatedegenesParamm
  
  
####------------------Validate deseqProperties Nodes---------------####
  
  validatedeseqPoisson <- validateCharacter(configFile$deseqProperties$poisson)
  
  configFile$deseqProperties$poisson <- validatedeseqPoisson

  
####--------------------ValidateCoincidences-----------------------####
  
  validateCoincidencesFileName <- validateCharacter(configFile$coincidences$coincidencesFileName)
  
  configFile$coincidences$coincidencesFileName <- validateCoincidencesFileName
  
  
####-------------------Validate Output------------------------####

  validateGraph <- validateCharacter(configFile$output$graph)

  configFile$output$graph <- validateGraph

  validateGraphExpr <- validateCharacter(configFile$output$graphType$expr)

  configFile$output$graphType$expr <- validateGraphExpr

  validateGraphMD <- validateCharacter(configFile$output$graphType$MD)

  configFile$output$graphType$MD <- validateGraphMD

  validateGraphBiotypes <- validateCharacter(configFile$output$graphType$biotypes)

  configFile$output$graphType$biotypes <- validateGraphBiotypes

  validateGraphPCA <- validateCharacter(configFile$output$graphType$PCA)

  configFile$output$graphType$PCA <- validateGraphPCA

  validateGraphq <- validateNumber(configFile$output$graphType$q)

  configFile$output$graphType$q <- validateGraphq

  validateOutputNameNOISeq  <- validateCharacter(configFile$output$outputNameNOISeq)

  configFile$output$outputNameNOISeq <- validateOutputNameNOISeq

  validateOutputDirNOISeq  <- validateCharacter(configFile$output$outputDirNOISeq)

  configFile$output$outputDirNOISeq <- validateOutputDirNOISeq
  
  validateOutputNameDESeq  <- validateCharacter(configFile$output$outputNameDESeq)
  
  configFile$output$outputNameDESeq <- validateOutputNameDESeq
  
  validateOutputDirDESeq  <- validateCharacter(configFile$output$outputDirDESeq)
  
  configFile$output$outputDirDESeq <- validateOutputDirDESeq
  
  validateOutputDirCoincidences  <- validateCharacter(configFile$output$outputDirCoincidences)
  
  configFile$output$outputDirCoincidences <- validateOutputDirCoincidences
  
  validateOutputNameCoincidences  <- validateCharacter(configFile$output$outputNameCoincidences)
  
  configFile$output$outputNameCoincidences <- validateOutputNameCoincidences
  
  return (configFile)
}


#' App structure
#' @name readConfigFile
#' @param root
#' @description App to start geoapp
#'
#' @import logging, XML
#' @return configFile
#'
#' @author Alberto Rodriguez-Izquierdo, 2021

readConfigFile <- function(root){

  require(XML)

  myDirectoryConfigFile <- paste0(root, 'config/')

  if (!file.exists(paste0(myDirectoryConfigFile,"configFile.xml"))){

    #log_error <- 'File path does not exist! Please review the path'
    #logerror(log_error)

    #stop(log_error)
  }

  configFile <- xmlParse(file = paste0(myDirectoryConfigFile,"configFile.xml"))
  configFile <- xmlToList(configFile)
#  configFile <- data.frame(configFile)

  return (configFile)
}



#' @name nodesValidation
#' @param configFile
#' @return configFile
#'
#' @author Alberto Rodriguez-Izquierdo, 2021


nodesValidation <- function(configFile){

  #Building list with principal and secondary nodes for validation

  principalNodes              <- c('dataPath','samplesFactorFilename','output')

  analysisNodes               <- c('noiseq','deseq')

  bFactors                    <- configFile[grepl("bFactor",names(configFile))]

  geneLengthNodes             <- c('geneLengthDir',
                                   'geneLengthFilename',
                                   'geneLengthAnnotation',
                                   'geneLengthExtAnnot',
                                   'geneLengthOutputFilename')

  condAnalysisNOISeqbioNodes  <- c('nclust',
                                   'k',
                                   'norm',
                                   'filter',
                                   'r',
                                   'random.seed')

  degenesParamNodes           <- c('q','m')

  outputNodes                 <- c('graph',
                                   'graphType',
                                   'outputNameDESeq',
                                   'outputNameNOISeq',
                                   'outputDirNOISeq',
                                   'outputDirDESeq',
                                   'outputDirCoincidences',
                                   'outputNameCoincidences')

  graphNodes                  <- c('expr','MD', 'biotypes', 'PCA','q')
  
  deseqPropertiesNodes        <- c('poisson')
  
  coincidencesNodes           <- c('coincidencesFileName')

  #Validation principal nodes

  generalParametersNodes        <- validateConfigNodes(principalNodes, configFile)

  ValAnalysisNodes              <- validateConfigNodes(analysisNodes, configFile$analysis)

  bFactorsNodes                 <- validateConfigNodes(bFactors,configFile)

  ValgeneLengthNodes            <- validateConfigNodes(geneLengthNodes, configFile$geneLength)

  ValcondAnalysisNOISeqbioNodes <- validateConfigNodes(condAnalysisNOISeqbioNodes, configFile$condAnalysisNOISeqbio)

  ValdegenesParamNodes          <- validateConfigNodes(degenesParamNodes, configFile$degenesParam)
  
  valdeseqPropertiesNodes       <- validateConfigNodes(deseqPropertiesNodes, configFile$deseqProperties)

  ValOutputNodes                <- validateConfigNodes(outputNodes, configFile$output)

  ValGraphNodes                 <- validateConfigNodes(graphNodes, configFile$output$graphType)

  ValCoincidencesNodes                <- validateConfigNodes(coincidencesNodes, configFile$coincidences)
  
  #Validation secondary nodes

  for (factors in names(bFactors)){

    factor <- eval(parse(text=paste0('configFile$',factors)))
    ##-----------------------Define name of samples------------------------##
    nameFactor1 <- strsplit(eval(parse(text=paste0('configFile$',factors, '$name_samples_1'))),',')
    factor$name_samples_1 <- nameFactor1

    nameFactor2 <- strsplit(eval(parse(text=paste0('configFile$',factors, '$name_samples_2'))),',')
    factor$name_samples_2 <- nameFactor2
    
    if (grepl(',', factor$factorDESeq)){
      
      nameFactorFactor <- strsplit(eval(parse(text=paste0('configFile$',factors, '$factorDESeq'))),',')
      factor$factorDESeq <- nameFactorFactor  
      
      eval(parse(text=paste0('configFile$',factors,'$factorDESeq <- ', factor$factorDESeq)))
      
    }

    ##----------------------------Splitting by ','-------------------------##
    eval(parse(text=paste0('configFile$',factors,'$name_samples_1 <- ',factor$name_samples_1)))

    eval(parse(text=paste0('configFile$',factors,'$name_samples_2 <- ',factor$name_samples_2)))
    
    ##------------------------------Validate nodes--------------------------##

    bFactorsLevels1     <- eval(parse(text=paste0('configFile$',factors,'$name_samples_1')))

    bFactorsLevels2     <- eval(parse(text=paste0('configFile$',factors,'$name_samples_2')))

    bFactorsLevelsFact <- eval(parse(text=paste0('configFile$',factors,'$factor')))

    eval(parse(text=paste0('validationNode', factors, 'sampleName <- validateConfigNodes(bFactorsLevels1,configFile$', factors,'$name_samples_1)')))

    eval(parse(text=paste0('validationNode', factors, 'sampleName <- validateConfigNodes(bFactorsLevels2,configFile$', factors,'$name_samples_2)')))

    eval(parse(text=paste0('validationNode', factors, 'factor <- validateConfigNodes(bFactorsLevelsFact,configFile$', factors,'$factor)')))
  }

  return(configFile)
  #loginfo('Nodes validation success!')
}




#' App structure
#' @name validateNodes
#' @param configFile
#' @description App to start geoapp
#'
#' @import logging, XML
#' @return configFile
#'
#' @author Alberto Rodríguez Izquierdo


validateConfigNodes <- function (nodes, configFile){

  if (!is.list(nodes)){
    for (x in nodes){
      if (!(x %in% names(configFile))){
        if (!(x %in% configFile)){

#         log_error <- paste0('Node ',x,' does not exist! Please review configFile.xml')
#         logerror(log_error)
          print(paste0('Warning: Node ', nodes,' does not exist'))
#         stop(log_error)
        }
      }
    }
  }else{
    for (x in names(nodes)){
      if (!(x %in% names(configFile))){

        #       log_error <- paste0('Node ',x,' does not exist! Please review configFile.xml')
        #       logerror(log_error)
        print(paste0('Warning: Node ', nodes,' does not exist'))
        #       stop(log_error)
      }
    }
  }
}


#' App structure
#' @name ValidateCharacter
#' @param configFile
#' @description App to start geoapp
#'
#' @import logging, XML
#' @return configFile
#'
#' @author Alberto Rodríguez Izquierdo


validateCharacter <- function(configFile){

  if (!is.null(configFile)){

    if(configFile == "TRUE"){

      configFile <- TRUE

    }else if(configFile == "FALSE"){

      configFile <- FALSE
    }
  }

  return(configFile)

}


#' App structure
#' @name ValidateNumber
#' @param configFile
#' @description App to start geoapp
#'
#' @import logging, XML
#' @return configFile
#'
#' @author Alberto Rodríguez Izquierdo



validateNumber <- function(configFile){



#  Error <- FALSE

  if (is.null(configFile)){

#    log_error <- paste0('Value ', configFile, ' is null. Please check configFile')
#    logerror(log_error)
    print('Warning: Node is null')
#    Error <- TRUE
  }else if (!as.numeric(configFile)){

#    log_error <- paste0('Value ', configFile, ' is not a number. Please check configFile')
#    logerror(log_error)
    print('Warning: Node is not numeric char')
#    Error <- TRUE
  }else{
    configFile <- as.numeric(configFile)
#  if (any(Error == TRUE)){
#    error <- Error
#    stop('Not possible to validate, please check configFile format')
#  }
  }
  return(configFile)
}

