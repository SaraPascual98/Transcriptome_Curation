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
conda activate stringtie

# Definition of variables

PARAMS_VRP=$(sed -n ${SLURM_ARRAY_TASK_ID}p /home/dcarrasco/Resultados/samples_vrp_path_list.txt) # List with the path of all samples of V.rupestris to merge 
PARAMS_VRP=$(sed -n ${SLURM_ARRAY_TASK_ID}p /home/dcarrasco/Resultados/samples_vbe_path_list.txt) # List with the path of all samples of V.rupestris to merge

SAMPLES_RES_VRP=/home/dcarrasco/Resultados/Resultados_ensmablaje/Vrp_merged # Output directory for V.rupestris merge
SAMPLES_RES_VBE=/home/dcarrasco/Resultados/Resultados_ensmablaje/Vbe_merged # Output directory for V.berlandieri merge

GFF_VRP=/home/arodriguez/Data_110R/Vrp/VITVix110R_v1.0.pseudomolecules.hap_Vrp.gff3 # V. rupestris annotation reference file in gff3 format
GFF_VBE=/home/arodriguez/Data_110R/Vbe/VITVix110R_v1.0.pseudomolecules.hap_Vbe.gff3 # V.berlandieri annotation reference file in gff3 format


# Merge assemblies of V.rupestris

srun stringtie $PARAMS_VRP \
--merge -G $GFF_VRP \
-o $SAMPLES_RES_VRP/110R_v1.1.hap_Vrp.merged.gtf \
-T 10 \
-l Vrp

# Merge assemblies of V.berlandieri

srun stringtie $PARAMS_VBE \
--merge -G $GFF_VBE \
-o $SAMPLES_RES_VBE/110R_v1.1.hap_Vbe.merged.gtf \
-T 10 \
-l Vbe