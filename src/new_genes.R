
#library import

library(dplyr)
library(GenomicRanges)
library(rtracklayer)


# Load the gff3 files

gff_vrp_new<- read.table("/home/dcarrasco/Resultados/transdecoder/VRP/R110_hap_vrp.gff3", header= FALSE, sep = "\t") # read the new transcriptome of Vrp haplotype
gff_vbe_new <- read.table("/home/dcarrasco/Resultados/transdecoder/VBE/R110_hap_vbe.gff3", header= FALSE, sep = "\t") # read the new transcriptome of Vbe haplotype

gff_vrp_ref<- read.table("/home/arodriguez/Data_110R/Vrp/VITVix110R_v1.0.pseudomolecules.hap_Vrp.gff3", header= FALSE, sep = "\t") # read the reference transcriptome of Vrp haplotype
gff_vbe_ref <- read.table("/home/arodriguez/Data_110R/Vbe/VITVix110R_v1.0.pseudomolecules.hap_Vbe.gff3", header= FALSE, sep = "\t") # read the reference transcriptome of Vbe haplotype


# filter gff by genes

genes_ref_vbe <- gff_vbe_ref[gff_vbe_refl$V3 == "gene",]
genes_ref_vrp <- gff_vrp_ref[gff_vrp_ref$V3 == "gene", ]

genes_new_vbe <- gff_vbe_new[gff_vbe_new$V3 == "gene", ]
genes_new_vrp <- gff_vrp_new[gff_vrp_new$V3 == "gene", ]

# make GRAnges object

genes_ref_vbe_gr <- makeGRangesFromDataFrame(genes_ref_vbe)
genes_ref_vrp_gr <- makeGRangesFromDataFrame(genes_ref_vrp)

genes_new_vbe_gr <- makeGRangesFromDataFrame(genes_new_vbe)
genes_new_vrp_gr <- makeGRangesFromDataFrame(genes_new_vrp)

# retain only the genes from new gffs that not overlapping with reference 

no_overlapping_genes_vbe <- subsetByOverlaps(genes_new_vbe_gr, genes_ref_vbe_gr, type=any, invert=TRUE)
no_overlapping_genes_vrp <- subsetByOverlaps(genes_new_vrp_gr, genes_ref_vrp_gr, type=any, invert=TRUE)

# save no overlapping genes as gff file

export.gff3(no_overlapping_genes_vbe, "/home/dcarrasco/Resultados/new_genes/new_genes_vbe.gff3")
export.gff3(no_overlapping_genes_vrp, "/home/dcarrasco/Resultados/new_genes/new_genes_vrp.gff3")
