#!/bin/bash
#SBATCH -t  03:00:00 # execution time 3 hours
#SBATCH -c 12 # number of tasks
#SBATCH --mem-per-cpu=2G

module load cesga/2020 gcccore/system kraken2/2.1.2
kraken2-build --threads $SLURM_CPUS_PER_TASK --download-taxonomy --db DB_Kraken2/
kraken2-build --threads $SLURM_CPUS_PER_TASK --download-library viral --db DB_Kraken2/
kraken2-build --build --db DB_Kraken2/
kraken2 --db DB_Kraken2/ --gzip-compressed --paired E1-131222_S1_L001_R1_001.fastq.gz E1-131222_S1_L001_R2_001.fastq.gz --classified-out cseqs#.fq seqs_1.fq seqs_2.fq
gzip seqs_1.fq
gzip seqs_2.fq



