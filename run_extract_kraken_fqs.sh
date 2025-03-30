#!/bin/bash
#SBATCH --mem-per-cpu=2G 	# Memory
#SBATCH -t 00:30:00  	# Run time (d-hh:mm:ss) where days are optional --> 30 min
#SBATCH -c 12		# Cores per task requested

module load cesga/2020 gcccore/system kraken2/2.1.2

python extract_kraken_reads.py -k E1_kraken -s1 E1-131222_S1_L001_R1_001.fastq.gz -s2 E1-131222_S1_L001_R1_001.fastq.gz -t 10239 -o New_R1.fq.gz -o2 New_R2.fq.gz &> pipeline_log.txt