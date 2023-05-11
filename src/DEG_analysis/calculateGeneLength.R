#' @name calculateGeneLength
#' @param annotationFile
#' @param formatAnnot
#' @param fileNameDir
#' @param fileName
#' @import genomicFeatures, BiocManager
#'
#' @return geneLength
#'
#' @author Alberto Rodriguez-Izquierdo, 2021



calculateGeneLength <- function(annotationFile, formatAnnot, fileNameDir,fileName){

  txdb <- makeTxDbFromGFF(annotationFile,format=formatAnnot)

  exons.list.per.gene <- exonsBy(txdb,by="gene")

  exonic.gene.sizes <- bplapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})

  unlist_geneLength<-unlist(exonic.gene.sizes)

  write.table(unlist_geneLength,paste0(fileNameDir, fileName))

  return(unlist_geneLength)

}
