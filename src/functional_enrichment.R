# package import

library(stringr)
library(tidyr)
library(AnnotationForge)


# Read eggnog-mapper results and transcriptome Vbe v1.1 and Vrp V1.1 gff files 

egg_vbe_new <- read.table("/home/dcarrasco/Resultados/eggnog-mapper/Vbe_new/110_R_Vbe.emapper.annotations", sep= "\t", fill=TRUE) # functional annotation for vbe haplotype 
egg_vrp_new <- read.table("/home/dcarrasco/Resultados/eggnog-mapper/Vrp_new/110_R_Vrp.emapper.annotations", sep= "\t", fill=TRUE) # functional annotation for vrp haplotype

GFF_VRP_NEW <- read.table("home/dcarrasco/Resultados/transdecoder/VRP/R110_hap_vrp.gff3", sep= "\t") # transcriptome of Vrp haplotype
GFF_VRP_NEW_gene <- GFF_VRP_NEW[GFF_VRP_NEW$V3 =="gene",] # filter by gene
GFF_VRP_NEW_gene  <- GFF_VRP_NEW_gene [,-c(2:8)] # get only the gene id and chromosome info
GFF_VRP_NEW_gene$V9 <-str_extract(GFF_VRP_NEW_gene$V9, "(?<=\\=)(.*?)(?=\\;N)") # extract the gene id to match the results of the eggnog-mapper


GFF_VBE_NEW <- read.table("/home/dcarrasco/Resultados/transdecoder/VBE/R110_hap_vbe.gff3", sep= "\t") # transcriptome of Vbe haplotype
GFF_VBE_NEW_gene <- GFF_VBE_NEW[GFF_VBE_NEW$V3 =="gene",] # filter by gene
GFF_VBE_NEW_gene  <- GFF_VBE_NEW_gene [,-c(2:8)] # get only the gene id and chromosome info
GFF_VBE_NEW_gene$V9 <-str_extract(GFF_VBE_NEW_gene$V9, "(?<=\\=)(.*?)(?=\\;N)") # extract the gene id to match the results of the eggnog-mapper


# filter and retain only the genes that have GOs associated

egg_vbe_new_1 <- egg_vbe_new[egg_vbe_new$V10 != "-",]
egg_vbe_new_1$V1 <- str_extract(vbe_new_1$V1, "(?<=\\=)(.*?)(?=\\;N)")

egg_vrp_new_1 <- egg_vrp_new[egg_vrp_new$V10 != "-",]
egg_vrp_new_1$V1 <- str_extract(vrp_new_1$V1, "(?<=\\=)(.*?)(?=\\;N)")


# get the fSym ("GID", "SYMBOL", "GENENAME"), fGo ("GID", "GO", "IEA") and fChr ("GID", "CHROMOSOME") files for V.rupestris


fSym_vrp_new <- cbind(egg_vrp_new_1$V1, egg_vrp_new_1$V2, egg_vrp_new_1$V8) 
colnames(fSym_vrp_new) <- c("GID", "SYMBOL", "GENENAME")

fChr_vrp_new <- GFF_VRP_NEW_gene[GFF_VRP_NEW_gene$V9 %in% egg_vrp_new_1$V1,]
fChr_vrp_new <- cbind(fChr_vrp_new$V9, fChr_vrp_new$V1)
colnames(fChr_vrp_new) <- c("GID", "CHROMOSOME")  


fGO_vrp_new <- cbind(egg_vrp_new_1$V1, egg_vrp_new_1$V10)
colnames(fGO_vrp_new) <- c("GID", "GO") 
fGO_vrp_new <- data.frame(fGO_vrp_new)
fGO_vrp_new <- separate_rows(fGO_vrp_new, GO, sep = ",")
fGO_vrp_new$EVIDENCE <- "IEA"


# get the fSym, fGo and fChr files for V.berlandieri


fSym_vbe_new <- cbind(egg_vbe_new_1$V1, egg_vbe_new_1$V2, egg_vbe_new_1$V8)
colnames(fSym_vbe_new) <- c("GID", "SYMBOL", "GENENAME")

fChr_vbe_new <- GFF_VBE_NEW_gene[GFF_VBE_NEW_gene$V9 %in% egg_vbe_new_1$V1,]
fChr_vbe_new <- cbind(fChr_vbe_new$V9, fChr_vbe_new$V1)
colnames(fChr_vbe_new) <- c("GID", "CHROMOSOME")  

fGO_vbe_new <- cbind(vbe_new_1$V1, vbe_new_1$V10)
colnames(fGO_vbe_new) <- c("GID", "GO") 
fGO_vbe_new <- data.frame(fGO_vbe_new)
fGO_vbe_new <- separate_rows(fGO_vbe_new, GO, sep = ",")
fGO_vbe_new$EVIDENCE <- "IEA"

# Read GOs list from different plants for filter GOs that are not in plants

go_arab <- read.table("/home/dcarrasco/Resultados/GO_to_cross/Athaliana_mart_export.txt", sep= "\t", fill=TRUE)
go_ptri <- read.table("/home/dcarrasco/Resultados/GO_to_cross/Ptrichocarpa_mart_export.txt", sep= "\t", fill=TRUE)
go_vitis <- read.table("/home/dcarrasco/Resultados/GO_to_cross/vitis_mart_export.txt", sep= "\t", fill=TRUE)
go_slyco <- read.table("/home/dcarrasco/Resultados/GO_to_cross/Slycopersicum_mart_export.txt", sep= "\t", fill=TRUE)
go_zmays <- read.table("/home/dcarrasco/Resultados/GO_to_cross/Zmays_mart_export.txt", sep= "\t", fill=TRUE)


# Filter GOs from fsym, fChr and fGO files of V.rupestris with the previous lists

fGO_vrp_new_1 <- fGO_vrp_new[fGO_vrp_new$GO %in% go_arab$V3 | fGO_vrp_new$GO %in% go_ptri$V3 | fGO_vrp_new$GO %in% go_vitis$V3 |
                               fGO_vrp_new$GO %in% go_slyco$V3 | fGO_vrp_new$GO %in% go_zmays$V3, ]

fSym_vrp_new_1 <- fSym_vrp_new[fSym_vrp_new$GID %in% fGO_vrp_new_1$GID,]
fChr_vrp_new_1 <- fChr_vrp_new[fChr_vrp_new$GID %in% fGO_vrp_new_1$GID,]

fGO_vrp_new_1 <- data.frame(fGO_vrp_new_1)
fSym_vrp_new_1  <- data.frame(fSym_vrp_new_1)
fChr_vrp_new_1 <- data.frame(fChr_vrp_new_1)


# Filter GOs from fsym, fChr and fGO files of V.berlandieri with the previous lists

fGO_vbe_new_1 <- fGO_vbe_new[fGO_vbe_new$GO %in% go_arab$V3 | fGO_vbe_new$GO %in% go_ptri$V3 | fGO_vbe_new$GO %in% go_vitis$V3 |
                               fGO_vbe_new$GO %in% go_slyco$V3 | fGO_vbe_new$GO %in% go_zmays$V3, ]

fSym_vbe_new_1 <- fSym_vbe_new[fSym_vbe_new$GID %in% fGO_vbe_new_1$GID,]
fChr_vbe_new_1 <- fChr_vbe_new[fChr_vbe_new$GID %in% fGO_vbe_new_1$GID,]

fGO_vbe_new_1 <- data.frame(fGO_vbe_new_1)
fSym_vbe_new_1  <- data.frame(fSym_vbe_new_1)
fChr_vbe_new_1 <- data.frame(fChr_vbe_new_1)



 
# Join the files of the two haplotypes
fGO_R110 <- rbind(fGO_vbe_new_1, fGO_vrp_new_1)
fsym_R110 <- rbind(fSym_vbe_new_1, fSym_vrp_new_1)
fChr_R110 <- rbind (fChr_vbe_new_1, fChr_vrp_new_1)



# create gene annotation data package of R110 v.1.1

makeOrgPackage(gene_info=fsym_R110, 
                            go=fGO_R110, 
                            chromosome=fChr_R110, 
                            version = "1.1", 
                            maintainer = "sara <sara.pascuale@estudiante.uam.es>",
                            author="sara <sara.pascuale@estudiante.uam.es>",
                            tax_id= "103352",
                            genus= "Vitis", 
                            species = "R110New", 
                            goTable = "go")   
                            



# Read eggnog-mapper results and transcriptome Vbe v1.0 and Vrp V.1.0 gff files 


egg_vbe_ref <- read.table("/media/sara/easystore/RESULTADOS_DEF/eggnog_mapper/resultados/vbe_old_gene/110_R_Vbe.emapper.annotations", sep= "\t", fill=TRUE)
egg_vrp_ref <- read.table("/media/sara/easystore/RESULTADOS_DEF/eggnog_mapper/resultados/vrp_old_gene/110_R_Vrp.emapper.annotations", sep= "\t", fill=TRUE)

GFF_VBE_OLD<- read.table("/media/sara/easystore/RESULTADOS_DEF/genoma_original/VITVix110R_v1.0.pseudomolecules.hap_Vbe.gff3", sep= "\t")
GFF_VBE_OLD_gene <- GFF_VBE_OLD[GFF_VBE_OLD$V3 =="gene",]
GFF_VBE_OLD_gene  <- GFF_VBE_OLD_gene[,-c(2:8)]
GFF_VBE_OLD_gene$V9 <-str_extract(GFF_VBE_OLD_gene$V9, "(?<=\\=)(.*?)(?=\\;N)")

GFF_VRP_OLD<- read.table("/media/sara/easystore/RESULTADOS_DEF/genoma_original/VITVix110R_v1.0.pseudomolecules.hap_Vrp.gff3", sep= "\t")
GFF_VRP_OLD_gene <- GFF_VRP_OLD[GFF_VRP_OLD$V3 =="gene",]
GFF_VRP_OLD_gene  <- GFF_VRP_OLD_gene[,-c(2:8)]
GFF_VRP_OLD_gene$V9 <-str_extract(GFF_VRP_OLD_gene$V9, "(?<=\\=)(.*?)(?=\\;N)")GFF_VRP_OLD<- read.table("/media/sara/easystore/RESULTADOS_DEF/genoma_original/VITVix110R_v1.0.pseudomolecules.hap_Vrp.gff3", sep= "\t")


# filter and retain only the genes that have GOs associated

egg_vbe_ref_1 <- egg_vbe_ref[egg_vbe_ref$V10 != "-",]
egg_vrp_ref_1 <- egg_vrp_ref[egg_vrp_ref$V10 != "-",]

# get the fSym ("GID", "SYMBOL", "GENENAME"), fGo ("GID", "GO", "EVIDENCE") and fChr ("GID", "CHROMOSOME") files for V.rupestris

fSym_vrp_old <- cbind(egg_vrp_ref_1$V1, egg_vrp_ref_1$V2, egg_vrp_ref_1$V8)
colnames(fSym_vrp_old) <- c("GID", "SYMBOL", "GENENAME")


fChr_vrp_old <- GFF_VRP_OLD_gene[GFF_VRP_OLD_gene$V9 %in% fSym_vrp_old$GID,]
fChr_vrp_old  <- cbind(fChr_vrp_old$V9, fChr_vrp_old$V1)
colnames(fChr_vrp_old ) <- c("GID", "CHROMOSOME") 

fGO_vrp_old <- cbind(egg_vrp_ref_1$V1, egg_vrp_ref_1$V10)
colnames(fGO_vrp_old) <- c("GID", "GO") 
fGO_vrp_old<- separate_rows(fGO_vrp_old, GO, sep = ",")
fGO_vrp_old$EVIDENCE <- "IEA"



# get the fSym, fGo and fChr files for V.berlandieri


fSym_vbe_old <- cbind(egg_vbe_ref_1$V1, egg_vbe_ref_1$V2, egg_vbe_ref_1$V8)
colnames(fSym_vbe_old) <- c("GID", "SYMBOL", "GENENAME")

fChr_vbe_old <- GFF_VBE_OLD_gene[GFF_VBE_OLD_gene$V9 %in% fSym_vbe_old$GID,]
fChr_vbe_old <- cbind(fChr_vbe_old$V9, fChr_vbe_old$V1)
colnames(fChr_vbe_old) <- c("GID", "CHROMOSOME")  

fGO_vbe_old <- cbind(egg_vbe_ref_1$V1, egg_vbe_ref_1$V10)
colnames(fGO_vbe_old) <- c("GID", "GO")
fGO_vbe_old<- separate_rows(fGO_vbe_old, GO, sep = ",")
fGO_vbe_old$EVIDENCE <- "IEA"


# Filter GOs of V.rupestris with the GOs from different plants

fGO_vrp_old_1 <- fGO_vrp_old[fGO_vrp_old$GO %in% go_arab$V3 | fGO_vrp_old$GO %in% go_ptri$V3 | fGO_vrp_old$GO %in% go_vitis$V3 |
                               fGO_vrp_old$GO %in% go_slyco$V3 | fGO_vrp_old$GO %in% go_zmays$V3, ]


fSym_vrp_old_1 <- fSym_vrp_old[fSym_vrp_old$GID %in% fGO_vrp_old_1$GID,]
fChr_vrp_old_1 <- fChr_vrp_old[fChr_vrp_old$GID %in% fGO_vrp_old_1$GID,]




# Filter GOs of V.berlandieri from different plants

fGO_vbe_old_1 <- fGO_vbe_old[fGO_vbe_old$GO %in% go_arab$V3 | fGO_vbe_old$GO %in% go_ptri$V3 | fGO_vbe_old$GO %in% go_vitis$V3 |
                               fGO_vbe_old$GO %in% go_slyco$V3 | fGO_vbe_old$GO %in% go_zmays$V3, ]

fSym_vbe_old_1 <- fSym_vbe_old[fSym_vbe_old$GID %in% fGO_vbe_old_1$GID,]
fChr_vbe_old_1 <- fChr_vbe_old[fChr_vbe_old$GID %in% fGO_vbe_old_1$GID,]


fsym_R110_ref <- rbind(fSym_vbe_old_1, fSym_vrp_old_1)
fchr_R110_ref <-rbind(fChr_vbe_old_1, fChr_vrp_old_1)
fGO_R110_ref <- rbind(fGO_vbe_old_1, fGO_vrp_old_1)

# create gene annotation data package of R110 v1.0

 makeOrgPackage(gene_info=fsym_R110_ref, 
                            go=fGO_R110_ref, 
                            chromosome=fchr_R110_ref, 
                            version = "0.1", 
                            maintainer = "sara <sara.pascuale@estudiante.uam.es>",
                            author="sara <sara.pascuale@estudiante.uam.es>",
                            tax_id= "103352",
                            genus= "Vitis", 
                            species = "R110Ref", 
                            goTable = "go")


#install and import the previously created R110 v.1.1 and R110 v.1.0 transcriptome functional annotation packages

install.packages("org.VR110New.eg.db", repos=NULL)
install.packages("org.VR110Ref.eg.db", respos=NULL)

library(org.VR110New.eg.db)
library(org.VR110Ref.eg.db)




### FUNCTIONAL ENRICHMENT ANALYSIS OF NEW, SPECIFIC AND TOTAL GENES OF EACH HAPLOTYPE ###

db_organism <- "org.VR110New.eg.db" 
vector_variables <- c('BP') # Only for biological process ontology
univ <- keys(org.VR110New.eg.db) # Use the annotation database from AnnotationForge

list_dirs <- list.dirs('/home/dcarrasco/Resultados/gene_list') # directory with 3 subdirectories with a list of genes inside each one: total genes of each haplotype, specific and new genes of each haplotype

for(x in list_dirs){
  
  if (x != 'home/dcarrasco/Resultados/gene_list'){
    
    if(file.exists(paste0(x,'/list_genes.txt'))){
      a <- read.table(paste0(x,'/list_genes.txt'), sep='\t', header=TRUE)
      for(y in vector_variables){
        
        b <- new('GOHyperGParams', geneIds=a, universeGeneIds=univ, annotation= org.VberlandieriRef.eg.db, ontology=y)
        c <- hyperGTest(b)
        d <- summary(c)
        dir.create(paste0(x,'/GO_analysis'))
        dir.create(paste0(x,'/GO_analysis/stats'))
        dir.create(paste0(x, '/GO_analysis/images'))
        write.table(d, file=paste0(x,'/GO_analysis/stats/',y, 'stats.txt'), sep='\t', row.names=FALSE)
        
        if(nrow(d)>1){
          simMatrix <- eval(parse(text=paste0('calculateSimMatrix(d$GO', y,'ID, orgdb=db_organism, keytype="GO", ont=y,method="Rel")')))
          print('simMatrix done! Starting scores step...')
          eval(parse(text=paste0('scores <- setNames(-log10(d$Pvalue), d$GO', y,'ID)')))
          print('Scores step done! Starting reduced terms step...')
          reducedTerms <- reduceSimMatrix(simMatrix, scores=scores, threshold=0.7,orgdb=db_organism, keytype="GO")
          print('reduceSimMatrix done! Plotting results...')
          tiff(filename=paste0(x,'/GO_analysis/images/scatterplot_', y, '.tiff'), units='in', width=10, height=5, res=300)
          print('Starting heatmap plot...')
          tiff(filename=paste0(x,'/GO_analysis/images/heatmap_', y, '.tiff'), units='in', width=10, height=5, res=300)
          heatmapPlot(simMatrix, reducedTerms, annotateParent=TRUE, annotationLabel="parentTerm", fontsize=4)
          dev.off()
          print('Heatmap plot done! Starting treemap plot...')
          tiff(filename=paste0(x,'/GO_analysis/images/treemap_', y, '.tiff'), units='in', width=10, height=5, res=300)
          treemapPlot(reducedTerms)  
          dev.off()
          print('Treemap plot done! Starting wordcloudPlot...')
          tiff(filename=paste0(x,'/GO_analysis/images/wordcloudPlot_', y, '.tiff'), units='in', width=10, height=5, res=300)
          wordcloudPlot(reducedTerms, min.freq=1, colors="black")
          dev.off()
        }
      }
    }
  }
}
