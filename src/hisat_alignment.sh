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



module load SAMtools/1.9-GCC-8.2.0-2.31.1 HISAT2/2.1.0-foss-2019a

source /home/dcarrasco/miniconda3/bin/activate
conda activate hisat2

# Creation of the genome index of each haplotype

srun hisat2-build /home/arodriguez/Data_110R/Vrp/VITVix110R_v1.0.pseudomolecules.hap_Vrp.fasta vrp -o /home/arodriguez/Data_110R/Vrp/VITVix110R_v1.0.pseudomolecules.hap_Vrp
srun hisat2-build /home/arodriguez/Data_110R/Vbe/VITVix110R_v1.0.pseudomolecules.hap_Vrp.fasta vbe -o /home/arodriguez/Data_110R/Vbe/VITVix110R_v1.0.pseudomolecules.hap_Vbe

# Definition of variables

PARAMS=$(sed -n ${SLURM_ARRAY_TASK_ID}p /home/dcarrasco/Resultados/samples_list.txt) # List with the names of all samples to align
SAMPLES_RES_VRP=/home/dcarrasco/Resultados/Resultados_alineamiento/vrp  # Output directory for V.rupestris alignment
SAMPLES_RES_VBE=/home/dcarrasco/Resultados/Resultados_alineamiento/vbe # Output directory for V.berlandieri alignment

INDEX_DIR_VRP=/home/arodriguez/Data_110R/Vrp/VITVix110R_v1.0.pseudomolecules.hap_Vrp #Index V.rupestris directory
INDEX_DIR_VBE=/home/arodriguez/Data_110R/Vbe/VITVix110R_v1.0.pseudomolecules.hap_Vbe #Index V.rupestris directory

TRIMMED_DATA_VRP=/home/arodriguez/data_all_experiments/vrp # Directory of trimmed reads of V.rupestris
TRIMMED_DATA_VBE=/home/arodriguez/data_all_experiments/vbe # Directory of trimmed reads of V.berlandieri


# Creation of the V.rupestris alignment results folders

if [[ ! -f $SAMPLES_RES_VRP/$PARAMS ]]
then
        rm -rf $SAMPLES_RES_VRP/$PARAMS
        mkdir $SAMPLES_RES_VRP/$PARAMS
        echo "####-------------------Creating folders...---------------------####"
fi


echo $PARAMS

echo $SAMPLES_RES_VRP/$PARAMS/metrics_.txt
echo $INDEX_DIR_VRP
echo $TRIMMED_DATA_VRP/$PARAMS/output_forward_paired.fq
echo $TRIMMED_DATA_VRP/$PARAMS/output_reverse_paired.fq
echo $SAMPLES_RES_VRP/$PARAMS/alignment.sam


cd $SAMPLES_RES_VRP/$PARAMS

# Alignment of V.rupestris samples against the reference genome

srun hisat2 --dta --mp 2,1 -x $INDEX_DIR_VRP \
-1 $TRIMMED_DATA_VRP/$PARAMS/output_forward_paired.fq \
-2 $TRIMMED_DATA_VRP/$PARAMS/output_reverse_paired.fq \
-S $SAMPLES_RES_VRP/$PARAMS/alignment.sam \
--met-file $SAMPLES_RES_VRP/$PARAMS/metrics.txt \
--min-intronlen 30 \
--novel-splicesite-outfile $SAMPLES_RES_VRP/$PARAMS/novel-splicesites.txt \
--un-gz $SAMPLES_RES_VRP/$PARAMS \
--al-gz $SAMPLES_RES_VRP/$PARAMS \
--un-conc-gz $SAMPLES_RES_VRP/$PARAMS \
--al-conc-gz $SAMPLES_RES_VRP/$PARAMS

conda activate samtools

# convert the SAM format of the V.rupestris alignment to BAM format and sort it

srun samtools view -S -b  $SAMPLES_RES_VRP/$PARAMS/alignment.sam >  $SAMPLES_RES_VRP/$PARAMS/alignment.bam
rm -rf  $SAMPLES_RES_VRP/$PARAMS/alignment.sam
srun samtools sort $SAMPLES_RES_VRP/$PARAMS/alignment.bam -o $SAMPLES_RES_VRP/$PARAMS/alignment_sorted.bam


# Creation of the V.berlandieri alignment results folders

if [[ ! -f $SAMPLES_RES_VBE/$PARAMS ]]
then
        rm -rf $SAMPLES_RES_VBE/$PARAMS
        mkdir $SAMPLES_RES_VBE/$PARAMS
        echo "####-------------------Creating folders...---------------------####"
fi


echo $PARAMS

echo $SAMPLES_RES_VBE/$PARAMS/metrics_.txt
echo $INDEX_DIR_VBE
echo $TRIMMED_DATA_VBE/$PARAMS/output_forward_paired.fq
echo $TRIMMED_DATA_VBE/$PARAMS/output_reverse_paired.fq
echo $SAMPLES_RES_VBE/$PARAMS/alignment.sam



# Alignment of V.berlandieri samples against the reference genome

srun hisat2 --dta --mp 2,1 -x $INDEX_DIR_VBE \
-1 $TRIMMED_DATA_VBE/$PARAMS/output_forward_paired.fq \
-2 $TRIMMED_DATA_VBE/$PARAMS/output_reverse_paired.fq \
-S $SAMPLES_RES_VBE/$PARAMS/alignment.sam \
--met-file $SAMPLES_RES_VBE/$PARAMS/metrics.txt \
--min-intronlen 30 \
--novel-splicesite-outfile $SAMPLES_RES_VBE/$PARAMS/novel-splicesites.txt \
--un-gz $SAMPLES_RES_VBE/$PARAMS \
--al-gz $SAMPLES_RES_VBE/$PARAMS \
--un-conc-gz $SAMPLES_RES_VBE/$PARAMS \
--al-conc-gz $SAMPLES_RES_VBE/$PARAMS



conda activate samtools

# convert the SAM format of the V.berlandieri alignment to BAM format and sort it

srun samtools view -S -b  $SAMPLES_RES_VBE/$PARAMS/alignment.sam >  $SAMPLES_RES_VBE/$PARAMS/alignment.bam
rm -rf  $SAMPLES_RES_VBE/$PARAMS/alignment.sam
srun samtools sort $SAMPLES_RES_VBE/$PARAMS/alignment.bam -o $SAMPLES_RES_VBE/$PARAMS/alignment_sorted.bam

