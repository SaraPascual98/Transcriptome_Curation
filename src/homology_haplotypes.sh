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
conda activate blast


# Make blast databases for genes two haplotype

makeblastdb -in vrp_genes.fasta -parse_seqids -dbtype nucl -out /home/dcarrasco/Resultados/blast/db/vrp_genes
makeblastdb -in vbe_genes.fasta -parse_seqids -dbtype nucl -out /home/dcarrasco/Resultados/blast/db/vbe_genes


# Get similar sequences genes between two haplotypes

blastp -query  vrp_genes.fasta \
    -db  vbe_genes -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-20  > /home/dcarrasco/Resultados/blast/results_dbvbe_vrp.outfmt6
    
    
blastp -query  vbe_genes.fasta \
    -db vrp_genes  -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-20  > /home/dcarrasco/Resultados/blast/results_dbvrp_vbe.outfmt6
