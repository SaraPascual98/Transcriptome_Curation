
#!/bin/sh
#SBATCH -p fast
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=begin #Envía un correo cuando el trabajo inicia
#SBATCH --mail-type=end #Envía un correo cuando el trabajo finaliza
#SBATCH --mail-user=alberto.rodriguez@inia.es #Dirección a la que se envía

source /home/arodriguez/miniconda3/bin/activate

conda activate snps

# Get annotation from gene of R110 v.1.0 transcriptome

python emapper.py -i /home/dcarrasco/Resultados/fasta_genes/VITVix110R_v1.0.pseudomolecules.hap_Vbe_gene.fasta \
--output /home/dcarrasco/Resultados/eggnog-mapper/Vbe_ref \
--cpu 10 \
--itype CDS \
-m diamond

python emapper.py -i /home/dcarrasco/Resultados/fasta_genes/VITVix110R_v1.0.pseudomolecules.hap_Vrp_gene.fasta \
--output /home/dcarrasco/Resultados/eggnog-mapper/Vrp_ref \
--cpu 10 \
--itype CDS \
-m diamond

# Get annotation from gene of R110 v.1.1 transcriptome

python emapper.py -i /home/dcarrasco/Resultados/fasta_genes/VITVix110R_v1.0.pseudomolecules.hap_Vbe_gene.fasta \
--output /home/dcarrasco/Resultados/eggnog-mapper/Vbe_new \
--cpu 10 \
--itype CDS \
-m diamond

python emapper.py -i /home/dcarrasco/Resultados/fasta_genes/VITVix110R_v1.0.pseudomolecules.hap_Vbe_gene.fasta \
--output /home/dcarrasco/Resultados/eggnog-mapper/Vrp_new \
--cpu 10 \
--itype CDS \
-m diamond
