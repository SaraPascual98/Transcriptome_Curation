#' @name ProcessingDESeqHeatMap
#' @param myData
#' @param configFile
#' @import DESeq2
#' @import PoiClaClu
#'
#' @return value_noiseqbio
#'
#' @author Alberto Rodriguez-Izquierdo, 2021
#'


ProcessingDESeqHeatMap <- function(myData, configFile, dirOutput, outputPathName){
  
  vsd <- vst(myData, blind = TRUE)
  
  sampleDists <- dist(t(assay(vsd)))
  
  sampleDistsMatrix <- as.matrix (sampleDists)
  
  if (isTRUE(configFile$deseqProperties$poisson)){
    
    poisd <- PoissonDistance(t(counts(myData)))
    
    poisdMatrix <- as.matrix(poisd$dd)
    
    colors <- colorRampPalette(rev(brewer.pal(9,'Blues')))(255)
    
    rownames(poisdMatrix) <- myData$factor
    
    colnames(poisdMatrix) <- NULL
    
    
    tiff(filename=paste0(dirOutput,'/figures/',outputPathName,'_heatMap_poisson.tiff'),units="in", width=5, height=5, res=300)
    
    pheatmap(poisdMatrix,
             clustering_distance_rows = poisd$dd,
             clustering_distance_cols = poisd$dd,
             col = colors)
    
    dev.off()
    
  }else{
    
    tiff(filename=paste0(dirOutput,'/figures/',outputPathName,'_heatMap.tiff'),units="in", width=5, height=5, res=300)
    
    pheatmap(sampleDistsMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors)
    
    dev.off()
    
  }

}

#' @name ProcessingDESeqPCA
#' @param myData
#' @param configFile
#' @import DESeq2
#' @import ggplot2
#'
#' @return value_noiseqbio
#'
#' @author Alberto Rodriguez-Izquierdo, 2021
#'

ProcessingDESeqPCA <- function(myData, configFile, dirOutput, outputPathName){
  
  
  vsd <- vst(myData, blind = TRUE)
  
  mypca <- DESeq2::plotPCA(vsd, intgroup=configFile)
  
  percentVar <- round(100 * attr(data.frame(mypca$data$PC1),which = 'percentVar'))
  
  tiff(filename=paste0(dirOutput,'/figures/',outputPathName,'_PCA.tiff'),units="in", width=5, height=5, res=300)
  
  ggplot(mypca$data, aes(x = mypca$data$PC1, y = mypca$data$PC2, color = configFile)) +
    geom_point(size=3)+
    xlab(paste0('PC1: ', percentVar[1], '% variance')) +
    ylab(paste0('PC2: ', percentVar[1], '% variance')) +
    coord_fixed()
  
  dev.off()
}

#' @name processingDeseq
#' @param myData
#' @param configFile
#' @import DESeq2
#' @import ggplot
#'
#' @return value_noiseqbio
#'
#' @author Alberto Rodriguez-Izquierdo, 2021
#'


processingDeseq <- function(myData, configFile, sampleFactors){
  
  if(length(unique(sampleFactors))>1){
    
    for(columnFactor in sampleFactors){
      
      designParameters <- columnFactor
      
      designParameters <- paste(designParameters, columnFactor, sep=' + ')
      
    }  
  
    deseqProcessing <- DESeq(myData)
    
    deseqResults <- DESeq2::results(deseqProcessing)
    
  } 
  deseqProcessing <- DESeq(myData)
  
  deseqResults <- DESeq2::results(deseqProcessing, contrast = sampleFactors)
  
  return(deseqResults)
  
}
