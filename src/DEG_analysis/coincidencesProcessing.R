#' @name loadFactors
#' @param dataNoiseq
#' @param dataDeseq
#' @param dataCounts 
#' @return dataCoincidences
#' @author Alberto Rodriguez-Izquierdo, 2021


coincidences <- function(dataNoiseq, dataDeseq, dataCounts, configFile){
  
  #Selection coincidences between NOISeq and DESeq2
  
  dataNoiseq$gene_id <- row.names(dataNoiseq)
  
  dataDeseq$gene_id <- dataDeseq@rownames
  
  dataDeseq <- data.frame(dataDeseq)
  
  dataCounts$gene_id <- row.names(dataCounts)
  
  coincidencesDataAnalysis <- merge(dataNoiseq,dataDeseq, by = "gene_id")
  
  #Selection coincidences in raw data
  
  rawdataCoincidences <- merge(coincidencesDataAnalysis, dataCounts, by = "gene_id")
  
  listGenes_all <- read.csv(configFile$coincidences$coincidencesFileName, sep=";", header=TRUE)
  
  geneList <- listGenes_all$Final_v3_name_tentative_creation_of_new_genes
  
  geneList <- data.frame(geneList)
  
  geneList$gene_id <- geneList$geneList
  
  geneList$geneList <- NULL
  
  geneList$v1 <- listGenes_all$v1_name
  
  geneList$refseq <- listGenes_all$Refseq_name
  
  rawdataCoincidences$gene_id_v3 <- rawdataCoincidences$gene_id
  
  rawdataCoincidences$gene_id <- gsub ('(?<=)[[:punct:]][a-z][0-9][0-9]', '', rawdataCoincidences$gene_id, perl = TRUE)
  
  geneListCoincidence <- merge(geneList, rawdataCoincidences, by= 'gene_id')
  
  return(geneListCoincidence)
}