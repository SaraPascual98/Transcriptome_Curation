


# Read blast results between two haplotypes genes

blast_vbe <-read.table("/home/dcarrasco/Resultados/blast/results_dbvrp.outfmt6", sep="\t")
blast_vrp <-read.table("/home/dcarrasco/Resultados/blast/results_dbvbe.outfmt6", sep="\t")

# Get the genes that match with any sequence from the other haplotype
genes_vrp <- blast_vrp$V1
genes_vrp_1 <- blast_vbe$V2
genes_vbe <- blast_vbe$V1
genes_vbe_1 <- blast_vrp$V2

# Read files with total genes present in the transcriptomes

genes_totales_vrp <- read.table("/media/sara/easystore/RESULTADOS_DEF/rvvgo/listas_genes/genes_vrp_new/list_genes.txt", header=TRUE)
genes_totales_vbe <- read.table("/media/sara/easystore/RESULTADOS_DEF/rvvgo/listas_genes/genes_vbe_new/list_genes.txt", header = TRUE)

# Get only the gene name from blast results

genes_vrp <- gsub(".*ID=([^;]+);.*", "\\1", genes_vrp)
genes_vrp <-data.frame(genes_vrp)
genes_vrp_1<- gsub(".*ID=([^;]+);.*", "\\1", genes_vrp_1)
genes_vrp_1 <-data.frame(genes_vrp_1)

genes_vbe<- gsub(".*ID=([^;]+);.*", "\\1", genes_vbe)
genes_vbe <-data.frame(genes_vbe)
genes_vbe_1<- gsub(".*ID=([^;]+);.*", "\\1", genes_vbe_1)
genes_vbe_1 <-data.frame(genes_vbe_1)

# Get the genes from each haplotype transcriome that not are present in blast results

genes_unicos_vrp <- genes_totales_vrp[! genes_totales_vrp$GID %in% genes_vrp$genes_vrp,]
genes_unicos_vrp <- genes_totales_vrp[! genes_totales_vrp$GID %in% genes_vrp_1$genes_vrp_1,]
genes_unicos_vrp <-data.frame(genes_unicos_vrp)

 
genes_unicos_vbe <- genes_totales_vbe[! genes_totales_vbe$GID %in% genes_vbe$genes_vbe,]
genes_unicos_vbe <- genes_totales_vbe[! genes_totales_vbe$GID %in% genes_vbe_1$genes_vbe_1,]
genes_unicos_vbe <-data.frame(genes_unicos_vbe)
colnames(genes_unicos_vbe) <- "GID"
colnames(genes_unicos_vrp) <- "GID"


write.table(genes_unicos_vrp ,"/home/dcarrasco/Resultados/specific_genes/specific_genes_vrp.txt", sep="\t", row.names = FALSE)
write.table(genes_unicos_vbe ,"/home/dcarrasco/Resultados/specific_genes/specific_genes_vbe.txt", sep="\t", row.names = FALSE)
