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


HAPLOTYPE=$(sed -n ${SLURM_ARRAY_TASK_ID}p /home/dcarrasco/Resultados/haplotype.txt) # Haplotype list
PARAMS=$(sed -n ${SLURM_ARRAY_TASK_ID}p /home/dcarrasco/Resultados/samples_list.txt) # List with the names of all samples to align
SAMPLES_RES=/home/dcarrasco/Resultados/Resultados_alineamiento
INDEX_DIR=/home/arodriguez/Data_110R/Vrp/VITVix110R_v1.0.pseudomolecules.hap_Vrp
TRIMMED_DATA=/home/arodriguez/data_all_experiments

if [[ ! -f $SAMPLES_RES/$PARAMS ]]
then
        rm -rf $SAMPLES_RES/$HAPLOTYPE/$PARAMS
        mkdir $SAMPLES_RES/$HAPLOTYPE/$PARAMS
        echo "####-------------------Creating folders...---------------------####"
fi


echo $PARAMS

echo $SAMPLES_RES/$HAPLOTYPE/$PARAMS/metrics_.txt
echo $INDEX_DIR
echo $TRIMMED_DATA/$HAPLOTYPE/$PARAMS/output_forward_paired.fq
echo $TRIMMED_DATA/$HAOLOTYPE/$PARAMS/output_reverse_paired.fq
echo $SAMPLES_RES/$HAPLOTYPE/$PARAMS/alignment.sam


cd $SAMPLES_RES/$HAPLOTYPE/$PARAMS


srun hisat2 --dta --mp 2,1 -x $INDEX_DIR \
-1 $TRIMMED_DATA/$HAPLOTYPE/$PARAMS/output_forward_paired.fq \
-2 $TRIMMED_DATA/$HAPLOTYPE/$PARAMS/output_reverse_paired.fq \
-S $SAMPLES_RES/$HAPLOTYPE/$PARAMS/alignment.sam \
--met-file $SAMPLES_RES/$HAPLOTYPE/$PARAMS/metrics.txt \
--min-intronlen 30 \
--novel-splicesite-outfile $SAMPLES_RES/$HAPLOTYPE/$PARAMS/novel-splicesites.txt \
--un-gz $SAMPLES_RES/$HAPLOTYPE/$PARAMS \
--al-gz $SAMPLES_RES/$HAPLOTYPE/$PARAMS \
--un-conc-gz $SAMPLES_RES/$HAPLOTYPE/$PARAMS \
--al-conc-gz $SAMPLES_RES/$HAPLOTYPE/$PARAMS

conda activate samtools

srun samtools view -S -b  $SAMPLES_RES/$HAPLOTYPE/$PARAMS/alignment.sam >  $SAMPLES_RES/$HAPLOTYPE/$PARAMS/alignment.bam
rm -rf  $SAMPLES_RES/$HAPLOTYPE/$PARAMS/alignment.sam
srun samtools sort $SAMPLES_RES/$HAPLOTYPE/$PARAMS/alignment.bam -o $SAMPLES_RES/$HAPLOTYPE/$PARAMS/alignment_sorted.bam

