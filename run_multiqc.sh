#!/bin/bash
#SBATCH --mem-per-cpu=2G 	# Memory
#SBATCH -t 00:10:00  	# Run time (d-hh:mm:ss) where days are optional --> 10 min
#SBATCH -c 12		# Cores per task requested

unset DISPLAY

module load cesga/2020 gcccore/system gcc/system fastqc/0.12.1 fastp/0.20.1 multiqc/1.13 openmpi/4.1.1_ft3 ragtag/2.1.0 samtools/1.9 varscan/2.3.9 quast/5.2.0 medusa/20240116 bedtools/2.31.0 unicycler/0.5.0 bwa-mem2/2.2.1 qualimap/2.3 multiqc/1.10.1 freebayes/1.3.7

. multiqc