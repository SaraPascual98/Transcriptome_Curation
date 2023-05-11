
#' @name ProcessingNOISeqbio
#' @param data
#' @param factor_nb
#' @param factor
#' @param configFile
#' @import NOISeq
#' @import dplyr
#'
#' @return value_noiseqbio
#'
#' @author Alberto Rodriguez-Izquierdo, 2021
#'


ProcessingNOISeqbio <- function(data,factor_nb,factor, configFile){
  
  #-----------------specifying variables-------------------#
  nclust_val        <- configFile$condAnalysisNOISeqbio$nclust
  k_val             <- configFile$condAnalysisNOISeqbio$k
  norm_val          <- configFile$condAnalysisNOISeqbio$norm
  factor_val        <- eval(parse(text=paste0('configFile$',factor_nb,'$factorNOISeq')))
  
  conditions_val    <- eval(parse(text=paste0('unique(factor$',factor_val,')')))
  
#  conditions_val    <- eval(parse(text=paste0('c(configFile$',factor_nb,'$samples_1$factor_1,configFile$',factor_nb,'$samples_2$factor_2)')))
  filter_val        <- configFile$condAnalysisNOISeqbio$filter
  r_val             <- configFile$condAnalysisNOISeqbio$r
  random.seed_val   <- configFile$condAnalysisNOISeqbio$random.seed

  value_noiseqbio   <- noiseqbio(data,
                                 nclust = nclust_val,
                                 k = k_val,
                                 norm = norm_val,
                                 factor = factor_val,
                                 conditions = conditions_val,
                                 filter = filter_val,
                                 r = r_val,
                                 random.seed = random.seed_val)
  
  return(value_noiseqbio)

}


#' @name calculateDEGenes
#' @param data
#' @param configFile
#'
#' @import NOISeq
#' @import dplyr
#'
#' @return value_degenes
#'
#' @author Alberto Rodriguez-Izquierdo, 2021
#'


calculateDEGenes <- function(data,configFile){
  
  q_val <- configFile$degenesParam$q
  M_val <- configFile$degenesParam$m
  
  data_degenes <- degenes(data, q = q_val, M = M_val)

  return(data_degenes)
}
