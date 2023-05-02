#!/bin/sh
#SBATCH -p medium
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-98
#SBATCH --mail-type=begin #Envía un correo cuando el trabajo inicia
#SBATCH --mail-type=end #Envía un correo cuando el trabajo finaliza
#SBATCH --mail-user=sara.pascuale@estudiante.uam.es #Dirección a la que se envía


source /home/dcarrasco/miniconda3/bin/activate
conda activate stringtie


# Definition of variables

PARAMS=$(sed -n ${SLURM_ARRAY_TASK_ID}p /home/dcarrasco/Resultados/samples_list.txt)  # List with the names of all samples to assembly

SAMPLES_RES_VRP=/home/dcarrasco/Resultados/Resultados_ensmablaje/vrp # Output directory for V.rupestris assembly
SAMPLES_RES_VBE=/home/dcarrasco/Resultados/Resultados_ensamblaje/vbe # Output directory for V.berlandieri assembly

GFF_VRP=/home/arodriguez/Data_110R/Vrp/VITVix110R_v1.0.pseudomolecules.hap_Vrp.gff3 # V. rupestris annotation reference file in gff3 format
GFF_VBE=/home/arodriguez/Data_110R/Vbe/VITVix110R_v1.0.pseudomolecules.hap_Vbe.gff3 # V.berlandieri annotation reference file in gff3 format

BAM_FILES_VRP=/home/dcarrasco/Resultados/Resultados_alineamiento/vrp # V.rupestris alignment results in BAM format
BAM_FILES_VBE=/home/dcarrasco/Resultados/Resultados_alineamiento/vbe # V.berlandieri alignment results in BAM format

# Creation of the V.rupestris assembly results folders


if [[ ! -f $SAMPLES_RES_VRP/$PARAMS ]]
then
        rm -rf $SAMPLES_RES_VRP/$PARAMS
        mkdir $SAMPLES_RES_VRP/$PARAMS
        echo "####-------------------Creating folders...---------------------####"
fi


# Assembly transcripts of V.rupestris alignment

srun stringtie $BAM_FILE_VRP/$PARAMS/alignment_sorted.bam \
-o $SAMPLES_RES_VRP/$PARAMS/new_annotation.gtf \
-p 4 \
-G $GFF_VRP \
-C $SAMPLES_RES_VRP/$PARAMS/coverage_transcript.gtf \
-A $SAMPLES_RES_VRP/$PARAMS/abundance_genes.gtf



# Creation of the V.berlandieri assembly results folders

if [[ ! -f $SAMPLES_RES_VBE/$PARAMS ]]
then
        rm -rf $SAMPLES_RES_VBE/$PARAMS
        mkdir $SAMPLES_RES_VBE/$PARAMS
        echo "####-------------------Creating folders...---------------------####"
fi


# Assembly transcripts of V.berlandiri alignment

srun stringtie $BAM_FILE_VBE/$PARAMS/alignment_sorted.bam \
-o $SAMPLES_RES_VBE/$PARAMS/new_annotation.gtf \
-p 4 \
-G $GFF_VBE \
-C $SAMPLES_RES_VBE/$PARAMS/coverage_transcript.gtf \
-A $SAMPLES_RES_VBE/$PARAMS/abundance_genes.gtf



### Merge GTF of all samples resulting from de assembly ###

#Define new variables


LIST_GTF_VRP=$(sed -n ${SLURM_ARRAY_TASK_ID}p /home/dcarrasco/Resultados/samples_vrp_path_list.txt) # List with the path of all samples of V.rupestris to merge 
LIST_GTF_VBE=$(sed -n ${SLURM_ARRAY_TASK_ID}p /home/dcarrasco/Resultados/samples_vbe_path_list.txt) # List with the path of all samples of V.berlandieri to merge

ASSEMBLY_RES_VRP=/home/dcarrasco/Resultados/Resultados_ensmablaje/Vrp_merged # Output directory for V.rupestris merge
ASSEMBLY_RES_VBE=/home/dcarrasco/Resultados/Resultados_ensmablaje/Vbe_merged # Output directory for V.berlandieri merge


# Merge assemblies of V.rupestris

srun stringtie $LIST_GTF_VRP \
--merge -G $GFF_VRP \
-o $ASSEMBLY_RES_VRP/110R_v1.1.hap_Vrp.merged.gtf \
-T 10 \
-l Vrp

# Merge assemblies of V.berlandieri

srun stringtie $LIST_GTF_VBE \
--merge -G $GFF_VBE \
-o $ASSEMBLY_RES_VBE/110R_v1.1.hap_Vbe.merged.gtf \
-T 10 \
-l Vbe




