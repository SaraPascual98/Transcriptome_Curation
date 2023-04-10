library(stringr)
library(tidyr)
library(AnnotationForge)


# Read eggnog-mapper results and transcriptome gff files

egg_vbe_new <- read.table("/media/sara/easystore/RESULTADOS_DEF/eggnog_mapper/110_R_Vbe.emapper.annotations", sep= "\t", fill=TRUE)
egg_vrp_new <- read.table("/media/sara/easystore/RESULTADOS_DEF/eggnog_mapper/110_R_Vrp.emapper.annotations", sep= "\t", fill=TRUE)

GFF_VRP_NEW <- read.table("/media/sara/easystore/RESULTADOS_DEF/trancriptoma_nuevo/transcripts_VBE_def.fasta.transdecoder.genome_blast.gff3", sep= "\t")
GFF_VRP_NEW_gene <- GFF_VRP_NEW[GFF_VRP_NEW$V3 =="gene",]
GFF_VRP_NEW_gene  <- GFF_VRP_NEW_gene [,-c(2:8)]
GFF_VRP_NEW_gene$V9 <-str_extract(GFF_VRP_NEW_gene$V9, "(?<=\\=)(.*?)(?=\\;N)")


GFF_VBE_NEW <- read.table("/media/sara/easystore/RESULTADOS_DEF/trancriptoma_nuevo/transcripts_VRP_def.fasta.transdecoder.genome_blast.gff3", sep= "\t")
GFF_VBE_NEW_gene <- GFF_VBE_NEW[GFF_VBE_NEW$V3 =="gene",]
GFF_VBE_NEW_gene  <- GFF_VBE_NEW_gene [,-c(2:8)]
GFF_VBE_NEW_gene$V9 <-str_extract(GFF_VBE_NEW_gene$V9, "(?<=\\=)(.*?)(?=\\;N)")


# filter and retain only the genes that have GOs associated

egg_vbe_new_1 <- egg_vbe_new[egg_vbe_new$V10 != "-",]
egg_vbe_new_1$V1 <- str_extract(vbe_new_1$V1, "(?<=\\=)(.*?)(?=\\;N)")

egg_vrp_new_1 <- egg_vrp_new[egg_vrp_new$V10 != "-",]
egg_vrp_new_1$V1 <- str_extract(vrp_new_1$V1, "(?<=\\=)(.*?)(?=\\;N)")


# get the fSym, fGo and fChr files for V.rupestris

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

# Read GOs list from different plants

go_arab <- read.table("/media/sara/easystore/RESULTADOS_DEF/eggnog_mapper/GO_to_cross/Athaliana_mart_export.txt", sep= "\t", fill=TRUE)
go_ptri <- read.table("/media/sara/easystore/RESULTADOS_DEF/eggnog_mapper/GO_to_cross/Ptrichocarpa_mart_export.txt", sep= "\t", fill=TRUE)
go_vitis <- read.table("/media/sara/easystore/RESULTADOS_DEF/eggnog_mapper/GO_to_cross/vitis_mart_export.txt", sep= "\t", fill=TRUE)
go_slyco <- read.table("/media/sara/easystore/RESULTADOS_DEF/eggnog_mapper/GO_to_cross/Slycopersicum_mart_export.txt", sep= "\t", fill=TRUE)
go_zmays <- read.table("/media/sara/easystore/RESULTADOS_DEF/eggnog_mapper/GO_to_cross/Zmays_mart_export.txt", sep= "\t", fill=TRUE)


# Filter GOs of V.rupestris with the previous lists

fGO_vrp_new_1 <- fGO_vrp_new[fGO_vrp_new$GO %in% go_arab$V3 | fGO_vrp_new$GO %in% go_ptri$V3 | fGO_vrp_new$GO %in% go_vitis$V3 |
                               fGO_vrp_new$GO %in% go_slyco$V3 | fGO_vrp_new$GO %in% go_zmays$V3, ]

fSym_vrp_new_1 <- fSym_vrp_new[fSym_vrp_new$GID %in% fGO_vrp_new_1$GID,]
fChr_vrp_new_1 <- fChr_vrp_new[fChr_vrp_new$GID %in% fGO_vrp_new_1$GID,]

fGO_vrp_new_1 <- data.frame(fGO_vrp_new_1)
fSym_vrp_new_1  <- data.frame(fSym_vrp_new_1)
fChr_vrp_new_1 <- data.frame(fChr_vrp_new_1)


# Filter GOs of V.berlandieri with the previous lists

fGO_vbe_new_1 <- fGO_vbe_new[fGO_vbe_new$GO %in% go_arab$V3 | fGO_vbe_new$GO %in% go_ptri$V3 | fGO_vbe_new$GO %in% go_vitis$V3 |
                               fGO_vbe_new$GO %in% go_slyco$V3 | fGO_vbe_new$GO %in% go_zmays$V3, ]

fSym_vbe_new_1 <- fSym_vbe_new[fSym_vbe_new$GID %in% fGO_vbe_new_1$GID,]
fChr_vbe_new_1 <- fChr_vbe_new[fChr_vbe_new$GID %in% fGO_vbe_new_1$GID,]

fGO_vbe_new_1 <- data.frame(fGO_vbe_new_1)
fSym_vbe_new_1  <- data.frame(fSym_vbe_new_1)
fChr_vbe_new_1 <- data.frame(fChr_vbe_new_1)




# create gene annotation data package of V.berlandieri

makeOrgPackage(gene_info=fSym_vbe_new_1, 
                            go=fGO_vbe_new_1, 
                            chromosome=fChr_vbe_new_1 , 
                            version = "0.1", 
                            maintainer = "sara <sara.pascuale@estudiante.uam.es>",
                            author="sara <sara.pascuale@estudiante.uam.es>",
                            tax_id= "103352",
                            genus= "Vitis", 
                            species = "berlandierigeneNew", 
                            goTable = "go")
                            
# create gene annotation data package of V.rupestris

makeOrgPackage(gene_info=fSym_vrp_new_1, 
                            go=fGO_vrp_new_1, 
                            chromosome=fChr_vrp_new_1, 
                            version = "0.1", 
                            maintainer = "sara <sara.pascuale@estudiante.uam.es>",
                            author="sara <sara.pascuale@estudiante.uam.es>",
                            tax_id= "103352",
                            genus= "portainjerto", 
                            species = "R110", 
                            goTable = "go")
                            
 
# create gene annotation data package of R110

fGO_R110 <- rbind(fGO_vbe_new_1, fGO_vrp_new_1)
fsym_R110 <- rbind(fSym_vbe_new_1, fSym_vrp_new_1)
fChr_R110 <- rbind (fChr_vbe_new_1, fChr_vrp_new_1)


makeOrgPackage(gene_info=fsym_R110, 
                            go=fGO_R110, 
                            chromosome=fChr_R110, 
                            version = "0.1", 
                            maintainer = "sara <sara.pascuale@estudiante.uam.es>",
                            author="sara <sara.pascuale@estudiante.uam.es>",
                            tax_id= "103352",
                            genus= "portainjerto", 
                            species = "R110", 
                            goTable = "go")   
                            
 
