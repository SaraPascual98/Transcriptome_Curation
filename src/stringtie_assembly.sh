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

GFF_VRP=/home/arodriguez/Data_110R/Vrp/VITVix110R_v1.0.pseudomolecules.hap_Vrp.gff3 # Annotation reference file in gff3 format
GFF_VBE=/home/arodriguez/Data_110R/Vbe/VITVix110R_v1.0.pseudomolecules.hap_Vbe.gff3 # Annotation reference file in gff3 format

BAM_FILES_VRP=/home/dcarrasco/Resultados/Resultados_alineamiento/vrp # V.rupestris alignment results in BAM format
BAM_FILES_VBE=/home/dcarrasco/Resultados/Resultados_alineamiento/vbe # V.berlandieri alignment results in BAM format

# Creation of the V.rupestris assembly results folders


if [[ ! -f $SAMPLES_RES_VRP/$PARAMS ]]
then
        rm -rf $SAMPLES_RES_VRP/$PARAMS
        mkdir $SAMPLES_RES_VRP/$PARAMS
        echo "####-------------------Creating folders...---------------------####"
fi


# Assembly transcript of V.rupestris alignment

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


# Assembly transcript of V.berlandiri alignment

srun stringtie $BAM_FILE_VBE/$PARAMS/alignment_sorted.bam \
-o $SAMPLES_RES_VBE/$PARAMS/new_annotation.gtf \
-p 4 \
-G $GFF_VBE \
-C $SAMPLES_RES_VBE/$PARAMS/coverage_transcript.gtf \
-A $SAMPLES_RES_VBE/$PARAMS/abundance_genes.gtf

