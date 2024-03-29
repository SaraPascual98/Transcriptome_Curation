#!/bin/sh
#SBATCH -p medium
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=begin #Envía un correo cuando el trabajo inicia
#SBATCH --mail-type=end #Envía un correo cuando el trabajo finaliza
#SBATCH --mail-user=sara.pascuale@estudiante.uam.es #Dirección a la que se envía



source /home/dcarrasco/miniconda3/bin/activate
conda activate transdecoder 


# Definition of variables

TRANSCRIPTOME_VRP=/home/dcarrasco/Resultados/Resultados_ensmablaje/Vrp_merged/110R_v1.1.hap_Vrp.merged.gtf
TRANSCRIPTOME_VBE=/home/dcarrasco/Resultados/Resultados_ensmablaje/Vbe_merged/110R_v1.1.hap_Vbe.merged.gtf 

GENOME_VRP=/home/arodriguez/Data_110R/Vrp/VITVix110R_v1.0.pseudomolecules.hap_Vrp.fasta
GENOME_VBE=/home/arodriguez/Data_110R/Vbe/VITVix110R_v1.0.pseudomolecules.hap_Vbe.fasta

# Get fasta sequence of transcript

/home/dcarrasco/Resultados/Codigo/TransDecoder/TransDecoder/util/gtf_genome_to_cdna_fasta.pl $TRANSCRIPTOME_VBE $GENOME_VRP > /home/dcarrasco/Resultados/transdecoder/VRP/transcripts_VRP.fasta 
/home/dcarrasco/Resultados/Codigo/TransDecoder/TransDecoder/util/gtf_genome_to_cdna_fasta.pl  $TRANSCRIPTOME_VRP $GENOME_VBE > /home/dcarrasco/Resultados/transdecoder/VBE/transcripts_VBE.fasta 

# Get gff3 from gtf file

/home/dcarrasco/Resultados/Codigo/TransDecoder/TransDecoder/util/gtf_to_alignment_gff3.pl $TRANSCRIPTOME_VRP >/home/dcarrasco/Resultados/transdecoder/VRP/transcripts_VRP.gff3
/home/dcarrasco/Resultados/Codigo/TransDecoder/TransDecoder/util/gtf_to_alignment_gff3.pl $TRANSCRIPTOME_VBE > /home/dcarrasco/Resultados/transdecoder/VBE/transcripts_VBE.gff3


# Predicct long ORFs

/home/dcarrasco/Resultados/Codigo/TransDecoder/TransDecoder/util/TransDecoder.LongOrfs -t /home/dcarrasco/Resultados/transdecoder/transcripts_VRP.fasta > /home/dcarrasco/Resultados/transdecoder/VRP/
/home/dcarrasco/Resultados/Codigo/TransDecoder/TransDecoder/util/TransDecoder.LongOrfs -t /home/dcarrasco/Resultados/transdecoder/transcripts_VBE.fasta > /home/dcarrasco/Resultados/transdecoder/VBE/


# Get homology annotation from blast and pfam

blastp -query /home/dcarrasco/Resultados/transdecoder/VBE/longest_orfs.pep  \
    -db uniprot_sprot.fasta  -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 10 > blastp_vbe.outfmt6

blastp -query /home/dcarrasco/Resultados/transdecoder/VRP/longest_orfs.pep  \
    -db uniprot_sprot.fasta  -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 10 > blastp_vrp.outfmt6   
    
    
hmmsearch --cpu 8 -E 1e-10 --domtblout pfam.domtblout /path/to/Pfam-A.hmm /home/dcarrasco/Resultados/transdecoder/VBE/longest_orfs.pep > pfam_vbe.domtblout
hmmsearch --cpu 8 -E 1e-10 --domtblout pfam.domtblout /path/to/Pfam-A.hmm /home/dcarrasco/Resultados/transdecoder/VRP/longest_orfs.pep  > pfam_vrp.domtblout

# Predict cds

/home/dcarrasco/Resultados/Codigo/TransDecoder/TransDecoder/util/TransDecoder.Predict -t transcripts_VRP.fasta --retain_pfam_hits pfam_vrp.domtblout --retain_blastp_hits blastp_vrp.outfmt6> /home/dcarrasco/Resultados/transdecoder/VRP/

/home/dcarrasco/Resultados/Codigo/TransDecoder/TransDecoder/util/TransDecoder.Predict -t transcripts_VBE.fasta --retain_pfam_hits pfam_vbe.domtblout --retain_blastp_hits blastp_vbe.outfmt6> /home/dcarrasco/Resultados/transdecoder/VBE/




# generate a genome-based coding region annotation file

cd /home/dcarrasco/Resultados/transdecoder/VRP

/home/dcarrasco/Resultados/Codigo/TransDecoder/TransDecoder/util/cdna_alignment_orf_to_genome_orf.pl \
     transcripts_VRP.fasta.transdecoder.gff3 \
     transcripts_VRP.gff3 \
     transcripts_VRP.fasta > R110_hap_vrp.gff3
     
 cd /home/dcarrasco/Resultados/transdecoder/VBE    
     
/home/dcarrasco/Resultados/Codigo/TransDecoder/TransDecoder/util/cdna_alignment_orf_to_genome_orf.pl \
     transcripts_VBE.fasta.transdecoder.gff3 \
     transcripts_VBE.gff3 \
     transcripts_VBE.fasta> R110_hap_vbe.gff3
     
     
     
