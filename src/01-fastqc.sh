#!/bin/bash
#SBATCH --partition=256GBv1 		# partition name
#SBATCH --time=1-00:00 		# hours:minutes runlimit after which job will be killed
#SBATCH --job-name 03-fastqc		# Job name
#SBATCH -o ./log/20230109/03-fastqc-%j.out		# File to which standard out will be written
#SBATCH -e ./log/20230109/03-fastqc-%j.err 		# File to which standard err will be written

DATA_PATH="/project/greencenter/Toprak_lab/shared/TolC-Mutagenesis/data/20230109/illumina"
OUTPUT_DIR="/project/greencenter/Toprak_lab/shared/TolC-Mutagenesis/data/20230109/fastqc"
N_CORES=56

source /home2/s135322/.bashrc
conda activate toprak-ngs

FASTQ=($DATA_PATH/*.fastq.gz)
N_FASTQ=${#FASTQ[@]}

echo $N_FASTQ
for i in $(seq 0 $N_CORES $N_FASTQ); do
    INPUT=""
    #echo "-- $i"
    (( i_end=(($i+$N_CORES-1) < ($N_FASTQ-1)) ? ($i+$N_CORES-1) : ($N_FASTQ-1) ))
    for j in $(seq $i 1 $i_end); do
        INPUT+="${FASTQ[$j]} "
        #INPUT+="$j "
    done
    #echo $INPUT
    fastqc                   \
        -o $OUTPUT_DIR       \
        -t $N_CORES          \
        $INPUT
done