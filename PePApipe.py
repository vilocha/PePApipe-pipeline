#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, sys, time
import subprocess
import os
import time

## EXAMPLE OF COMMAND LINES TO RUN LOCALLY FOR A BATCH OF SAMPLES:

# while read line; do PePApipe_local.py -s $line -o ./$line -t 12 -r1 ./$line/*R1*fastq.gz -r2 ./$line/*R2*fastq.gz -R Georgia.fasta; done < ./asfv_samples.txt

# multiqc .


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@                                             FUNCTIONS                                                  @@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


####################################################################################################
###                              QUALITY AND TRIMMING                                            ###
####################################################################################################

def createSubdirectories(output):
    print("\n-------------------------------")
    print("### Creating Subdirectories ###")
    print("-------------------------------\n")
    time.sleep(2)
    parent_dir=output
    secondary_dir=os.path.join(parent_dir, "B_mappings_and_new_reads")
    tertiary_dir=os.path.join(parent_dir, "D_variants")
    fastp_dir=os.path.join(parent_dir, "A2_trimming_and_post_trimming_quality")
    dirty_dir=os.path.join(parent_dir, "A1_pre_trimming_quality")
    ragtag_fastp_dir=os.path.join(parent_dir, "C2_denovo2_ragtag")
    medusa_fastp_dir=os.path.join(parent_dir, "C3_denovo3_medusa")
    quast_fastp_dir=os.path.join(parent_dir, "C4_denovo4_quality")
  
    if not os.path.exists(secondary_dir):
      os.mkdir(secondary_dir)
    if not os.path.exists(tertiary_dir):
      os.mkdir(tertiary_dir)
    if not os.path.exists(fastp_dir):
      os.mkdir(fastp_dir)
    if not os.path.exists(dirty_dir):
      os.mkdir(dirty_dir)
    if not os.path.exists(ragtag_fastp_dir):
      os.mkdir(ragtag_fastp_dir)
    if not os.path.exists(medusa_fastp_dir):
      os.mkdir(medusa_fastp_dir)
    if not os.path.exists(quast_fastp_dir):
      os.mkdir(quast_fastp_dir)

def fastqc_pig_virus(dir, sample, threads, reads1, reads2):
    print("\n---------------------------")
    print("### Starting Pre-FastQC ###")
    print("---------------------------\n")
    time.sleep(2)
    parent_dir=dir
    out_dir=os.path.join(parent_dir, "A1_pre_trimming_quality")
    cmd="fastqc -o " + out_dir + " -t " + threads + " -f fastq " + reads1 + " " + reads2
    print(cmd)
    os.system(cmd)

def fastp_pig_virus(dir, sample, threads, reads1, reads2):
    print("\n----------------------")
    print("### Starting Fastp ###")
    print("----------------------\n")
    time.sleep(2)
    parent_dir=dir
    out_dir3=os.path.join(parent_dir, "A2_trimming_and_post_trimming_quality")
    fastp_trimmed_paired_reads1=out_dir3+"/"+sample+"_fastp_trimmed_1P.fastq.gz"
    fastp_trimmed_paired_reads2=out_dir3+"/"+sample+"_fastp_trimmed_2P.fastq.gz"
    fastp_trimmed_unpaired_reads=out_dir3+"/"+sample+"_fastp_trimmed_U.fastq.gz"
    fastp_html_file=out_dir3+"/"+sample+"_fastp_html_file"
    fastp_json_file=out_dir3+"/"+sample+"_fastp_json_file"
    cmd="fastp -i " +reads1+ " -I " +reads2+ " -q 20 -o " +fastp_trimmed_paired_reads1+ " -O " +fastp_trimmed_paired_reads2+ " --unpaired2 " +fastp_trimmed_unpaired_reads+ " -h " +fastp_html_file+ " -j " +fastp_json_file
    print(cmd)
    os.system(cmd)
    '''
    print("\n----------------------------------") 
    print("### Starting Post-Fastp FastQC ###")
    print("----------------------------------\n")
    cmd="fastqc -o " + out_dir3 + " -t " + threads + " -f fastq " + fastp_trimmed_paired_reads1 + " " + fastp_trimmed_paired_reads2
    print(cmd)
    os.system(cmd)
    '''
    

#############################################################################################
###        MAPPING CYCLE TO MAP TRIMMED RAW READS WITH ASFV VIRUS GENOME                  ### 
#############################################################################################

## FIRST OF ALL, IT IS NECESSARY TO INDEX THE VIRUS REFERENCE GENOME WITH THIS COMMAND LINE (IN THE FOLDER WHERE THE VIRUS REFERENCE GENOME IN .FASTA FORMAT IS LOCATED):
# bwa index virus_reference.fasta

def bwa_raw(dir, sample, reads3, reads4, referencev, threads="4"):
    print("\n---------------------------------------------------------------------------------------------------")
    print("###           Starting BWA alignment of trimmed raw reads with virus reference genome           ###")
    print("---------------------------------------------------------------------------------------------------\n")
    time.sleep(2)
    parent_dir=dir
    secondary_dir=os.path.join(parent_dir, "B_mappings_and_new_reads")
    s=" "
    cmd="bwa-mem2 mem -t " + threads + s + referencev + s + reads3 + s + reads4 + " > " + secondary_dir + "/" + sample + "_raw.sam"
    print(cmd)
    os.system(cmd)
    cmd="samtools view -S -b " + secondary_dir + "/" + sample + "_raw.sam -o " + secondary_dir + "/" + sample + "_raw.bam"
    print(cmd)
    os.system(cmd)
    cmd="samtools sort " + secondary_dir + "/" + sample + "_raw.bam -o " + secondary_dir + "/" + sample + "_raw.sort.bam"
    os.system(cmd)
    cmd="samtools index " + secondary_dir + "/" + sample + "_raw.sort.bam"
    os.system(cmd)
    cmd="samtools flagstat " + secondary_dir + "/" + sample + "_raw.sort.bam > " + secondary_dir + "/" + sample + "_raw.txt"
    print(cmd)
    os.system(cmd)     

## IF VIRUS CONTENT IS LOW (<85%) WE MUST CREATE A NEW SET OF FASTQ BY MAPPING TO Taxonomy ID: 10239 (VIRUSES) SO WE CAN WORK WITH VIRAL READS FROM NOW ON:  
    
def new_reads(dir, sample, bam):    
    print("\n----------------------------------------------------------------------------------")
    print("###           Creating and compressing fastq with mapped virus reads           ###")
    print("----------------------------------------------------------------------------------\n")
    time.sleep(2)
    parent_dir=dir
    secondary_dir=os.path.join(parent_dir, "B_mappings_and_new_reads")
    cmd="samtools view -b -F 12 " + secondary_dir + "/" + sample + "_raw.sort.bam > " + secondary_dir + "/" + sample + "_virus_mapped.bam"
    print(cmd)
    os.system(cmd)
    cmd="samtools sort -n " + secondary_dir + "/" + sample + "_virus_mapped.bam -o " + secondary_dir + "/" + sample + "_virus_mapped.sort.bam"
    os.system(cmd)
    cmd="bedtools bamtofastq -i " + secondary_dir + "/" + sample + "_virus_mapped.sort.bam -fq " + secondary_dir + "/" + sample + "_virus_R1.fq -fq2 " + secondary_dir + "/" + sample + "_virus_R2.fq"
    print(cmd)
    os.system(cmd)
    cmd="gzip " + secondary_dir + "/" + sample + "_virus_R1.fq " + secondary_dir + "/" + sample + "_virus_R2.fq "
    print(cmd)
    os.system(cmd)

def bwa_virus(dir, sample, reads5, reads6, referencev, threads="4"):
    print("\n----------------------------------------------------------------------------------------------------")
    print("###           Starting BWA alignment of mapped virus reads with virus reference genome           ###")
    print("----------------------------------------------------------------------------------------------------\n")
    time.sleep(2)
    parent_dir=dir
    secondary_dir=os.path.join(parent_dir, "B_mappings_and_new_reads")
    s=" "
    cmd="bwa-mem2 mem -t " + threads + s + referencev + s + reads5 + s + reads6 + " > " + secondary_dir + "/" + sample + "_virus.sam"
    print(cmd)
    os.system(cmd)
    cmd="samtools view -S -b " + secondary_dir + "/" + sample + "_virus.sam -o " + secondary_dir + "/" + sample + "_virus.bam"
    print(cmd)
    os.system(cmd)
    cmd="samtools sort " + secondary_dir + "/" + sample + "_virus.bam -o " + secondary_dir + "/" + sample + "_virus.sort.bam"
    os.system(cmd)
    cmd="samtools index " + secondary_dir + "/" + sample + "_virus.sort.bam"
    os.system(cmd)
    cmd="samtools flagstat " + secondary_dir + "/" + sample + "_virus.sort.bam > " + secondary_dir + "/" + sample + "_virus.txt"
    print(cmd)
    os.system(cmd)     


#################################################################################################
###                        QUALITY OF ASSEMBLY WITH REFERENCE                                 ###
#################################################################################################

def qualimap_raw(dir, sample):
    print("\n------------------------------------------------------------")
    print("### Starting quality check of BWA alignment of raw reads ###")
    print("------------------------------------------------------------\n")
    time.sleep(2)
    parent_dir=dir
    secondary_dir=os.path.join(parent_dir, "B_mappings_and_new_reads")
    cmd="qualimap bamqc -bam " + secondary_dir + "/" + sample + "_raw.sort.bam"
    print(cmd)
    os.system(cmd)

def qualimap_virus(dir, sample):
    print("\n--------------------------------------------------------------")
    print("### Starting quality check of BWA alignment of virus reads ###")
    print("--------------------------------------------------------------\n")
    time.sleep(2)
    parent_dir=dir
    secondary_dir=os.path.join(parent_dir, "B_mappings_and_new_reads")
    cmd="qualimap bamqc -bam " + secondary_dir + "/" + sample + "_virus.sort.bam"
    print(cmd)
    os.system(cmd)


#################################################################################################
###                                  DE-NOVO ASSEMBLY                                         ###
#################################################################################################
  
def unicycler_raw(dir, sample, reads3, reads4):
    print("\n--------------------------------------------------") 
    print("### Starting Unicycler Post-Fastp on raw reads ###")
    print("--------------------------------------------------\n")
    time.sleep(2)
    parent_dir=dir
    out_dir3=os.path.join(parent_dir, "C1_denovo1_unicycler")
    cmd="unicycler -1 " + reads3 + " -2 " + reads4 + " --mode bold --linear_seqs 1 -o " + out_dir3
    print(cmd)
    os.system(cmd)

def unicycler_virus(dir, sample, reads5, reads6):
    print("\n----------------------------------------------------") 
    print("### Starting Unicycler Post-Fastp on virus reads ###")
    print("----------------------------------------------------\n")
    time.sleep(2)
    parent_dir=dir
    out_dir3=os.path.join(parent_dir, "C1_denovo1_unicycler")
    cmd="unicycler -1 " + reads5 + " -2 " + reads6 + " --mode bold --linear_seqs 1 -o " + out_dir3
    print(cmd)
    os.system(cmd) 

def ragtag_both(dir, sample, referencev, fasta):
    print("\n----------------------------------------------------------------------------------") 
    print("### Implementing RagTag Post Unicycler-Fastp on two stages: correct + scaffold ###")
    print("----------------------------------------------------------------------------------\n")
    time.sleep(2)
    parent_dir=dir
    out_dir4=os.path.join(parent_dir, "C2_denovo2_ragtag")
    ragtag_output_file=out_dir4+"/ragtag.correct.fasta"
    cmd1="ragtag.py correct " + referencev + " ./" + sample + "/C1_denovo1_unicycler/" + fasta + " -o " + out_dir4
    cmd2="ragtag.py scaffold " + referencev + " " + ragtag_output_file + " -o " + out_dir4
    print(cmd1)
    os.system(cmd1)
    print(cmd2)
    os.system(cmd2)
	
def medusa_both(dir, sample):
    print("\n--------------------------------------------")
    print("### Starting MeDuSa Post Unicycler-Fastp ###")
    print("--------------------------------------------\n")
    time.sleep(2)
    parent_dir=dir
    out_dir7=os.path.join(parent_dir, "C3_denovo3_medusa")
    cmd1="cp ./" + sample + "/C2_denovo2_ragtag/ragtag.scaffold.fasta " + out_dir7
    cmd2=" | java -jar $EBROOTMEDUSA/medusa.jar -f ../references_other/ -i" + " ./" + sample + "/C3_denovo3_medusa/ragtag.scaffold.fasta -v -d -gexf -scriptPath $EBROOTMEDUSA/medusa_scripts -o " + out_dir7 + "/ragtag.medusa"
    cmd=cmd1+cmd2
    print(cmd)
    os.system(cmd)
    cmd="awk '/^>/{if(N)exit;++N;} {print;}' " + out_dir7 + "/ragtag.scaffold.fasta > " + out_dir7 + "/ragtag.final.fasta"
    print(cmd)
    os.system(cmd)
    cmd="awk '/^>/{if(N)exit;++N;} {print;}' " + out_dir7 + "/ragtag.medusa > " + out_dir7 + "/medusa.final.fasta"
    print(cmd)
    os.system(cmd)
    cmd1 = "grep '^>' " + out_dir7 + "/medusa.final.fasta > " + out_dir7 + "/medusa.final.fasta.ord | "
    cmd2="tail -n +2 " + out_dir7 + "/medusa.final.fasta | while read lines; do  echo $lines | rev >> " + out_dir7 + "/medusa.final.fasta.ord; done"
    cmd=cmd1+cmd2
    print(cmd)
    os.system(cmd)


#################################################################################################
###                          QUALITY OF DE-NOVO ASSEMBLY                                      ###
#################################################################################################

def quast_raw(dir, sample, scaffold, contigs, reads3, reads4):
    print("\n-----------------------------------------------------------")
    print("### Starting QUAST Post Unicycler-Fastp after raw reads ###")
    print("-----------------------------------------------------------\n")
    time.sleep(2)
    parent_dir=dir
    out_dir7=os.path.join(parent_dir, "C3_denovo3_medusa")
    medusa_output_file=out_dir7+"/medusa.final.fasta"
    quast_fastp_dir=os.path.join(parent_dir, "C4_denovo4_quality")
    if os.path.getsize(medusa_output_file) == 0:
      cmd="quast.py " + " -o " + quast_fastp_dir + " --pe1 " + reads3 + " --pe2 " + reads4 + " " + contigs
      print(cmd)
      os.system(cmd)
    else:
      cmd="quast.py " + " -o " + quast_fastp_dir + " --pe1 " + reads3 + " --pe2 " + reads4 + " " + scaffold
      print(cmd)
      os.system(cmd)

def quast_virus(dir, sample, scaffold, contigs, reads5, reads6):
    print("\n-------------------------------------------------------------")
    print("### Starting QUAST Post Unicycler-Fastp after virus reads ###")
    print("-------------------------------------------------------------\n")
    time.sleep(2)
    parent_dir=dir
    out_dir7=os.path.join(parent_dir, "C3_denovo3_medusa")
    medusa_output_file=out_dir7+"/medusa.final.fasta"
    quast_fastp_dir=os.path.join(parent_dir, "C4_denovo4_quality")
    if os.path.getsize(medusa_output_file) == 0:
      cmd="quast.py " + " -o " + quast_fastp_dir + " --pe1 " + reads5 + " --pe2 " + reads6 + " " + contigs
      print(cmd)
      os.system(cmd)
    else:
      cmd="quast.py " + " -o " + quast_fastp_dir + " --pe1 " + reads5 + " --pe2 " + reads6 + " " + scaffold
      print(cmd)
      os.system(cmd)


#################################################################################################
###                           VARIANT CALLING WITH VARSCAN2                                    ###
#################################################################################################
    
def varscan2_raw(dir, sample, bam, referencev, vcf, threads="4"):
    print("\n-----------------------------------------------------------")
    print("### Starting variant calling of raw reads with VARSCAN2 ###")
    print("-----------------------------------------------------------\n")
    time.sleep(2)
    parent_dir=dir
    tertiary_dir=os.path.join(parent_dir, "D_variants")
    s=" "
    cmd="samtools mpileup -f " + referencev + s + bam + " | java -jar " + "$EBROOTVARSCAN/VarScan.v2.3.9.jar mpileup2cns --min-coverage 20 --p-value 0.05 --min-var-freq 0.05 --min-avg-qual 20 --min-reads2 5 --variants --output-vcf 1 --ploidy 1 >" + s + tertiary_dir + "/" + sample +"_raw_variants.vcf"
    print(cmd)
    os.system(cmd)
        
def varscan2_virus(dir, sample, bam, referencev, vcf, threads="4"):
    print("\n-------------------------------------------------------------")
    print("### Starting variant calling of virus reads with VARSCAN2 ###")
    print("-------------------------------------------------------------\n")
    time.sleep(2)
    parent_dir=dir
    tertiary_dir=os.path.join(parent_dir, "D_variants")
    s=" "
    cmd="samtools mpileup -f " + referencev + s + bam + " | java -jar " + "$EBROOTVARSCAN/VarScan.v2.3.9.jar mpileup2cns --min-coverage 20 --p-value 0.05 --min-var-freq 0.05 --min-avg-qual 20 --min-reads2 5 --variants --output-vcf 1 --ploidy 1 >" + s + tertiary_dir + "/" + sample +"_virus_variants.vcf"
    print(cmd)
    os.system(cmd)
    

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@                                 MAIN                                             @@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

def main():
    #Start timer and initiate pipeline
    start_time = time.time()
    time.sleep(3)
    print("\n#####################################################################################################")
    print("##      'PePApipe': PIPELINE FOR BIOINFORMATIC ANALYSIS OF DNA FROM AFRICAN SWINE FEVER VIRUS      ##")
    print("#####################################################################################################\n")
    time.sleep(3)
    #Initiate the parser
    help_text="Pipeline for NGS analysis of illumina short reads from ASFV samples. This includes the following steps: raw reads quality analysis and trimming, denovo assembly, mapping with reference genome and variant calling"
  
    parser=argparse.ArgumentParser(description=help_text)
  
    required=parser.add_argument_group('required arguments')
  
    required.add_argument("--sample", "-s", help="sample name or id", required=True)
    required.add_argument("--output", "-o", help="output directory", required=True)
    parser.add_argument("--threads", "-t", help="number of threads")
    required.add_argument("--reads1", "-r1", help="fastq file with raw reads1 (paired-end)", required=True)
    required.add_argument("--reads2", "-r2", help="fastq file with raw reads2 (paired-end)", required=True)
    required.add_argument("--referencev", "-R", help="virus reference genome (fasta file)")
      
    args=parser.parse_args()
  
    if args.sample:
      print("sample = %s" % args.sample)
      sample=args.sample

    if args.output:
      print("output dir = %s" % args.output)
      output=args.output

    if args.threads:
      print("threads = %s" % args.threads)
      threads=args.threads
    
    if args.reads1:
      print("reads1 = %s" % args.reads1)
      reads1=args.reads1

    if args.reads2:
      print("reads2 = %s" % args.reads2)
      reads2=args.reads2
	
    if args.referencev:
      print("referencev = %s" % args.referencev)
      referencev=args.referencev
      
    if args.referencep:
      print("referencep = %s" % args.referencep)
      referencep=args.referencep
  
    #create subdirectories
    createSubdirectories(output)
  
    #run prefastqc
    fastqc_pig_virus(output, sample, threads, reads1, reads2)
  
    #run fastp
    fastp_pig_virus(output, sample, threads, reads1, reads2)
  
    ## MAPPING/ALIGNMENT OF RAW READS WITH VIRUS REFERENCE USING BWA ##
    #run bwa_post_fastp_first_cycle_raw
    bwa_raw(output, sample, output+"/A2_trimming_and_post_trimming_quality/"+sample+"_fastp_trimmed_1P.fastq.gz", output+"/A2_trimming_and_post_trimming_quality/"+sample+"_fastp_trimmed_2P.fastq.gz", referencev, threads)
        
    ## CREATION OF FASTQ FILES WITH READS MAPPED WITH VIRUS ##
    #run creation_of_new_fastq_files
    #new_reads(output, sample, output+"/B_mappings_and_new_reads/"+sample+"_virus_mapped.sort.bam")  
    
    ## MAPPING/ALIGNMENT OF VIRUS READS WITH VIRUS REFERENCE USING BWA ##
    #run bwa_post_fastp_second_cycle_virus
    #bwa_virus(output, sample, output+"/B_mappings_and_new_reads/"+sample+"_virus_R1.fq.gz", output+"/B_mappings_and_new_reads/"+sample+"_virus_R2.fq.gz", referencev, threads)
        
    #run unicycler_post_fastp raw reads
    unicycler_raw(output, sample, output+"/A2_trimming_and_post_trimming_quality/"+sample+"_fastp_trimmed_1P.fastq.gz", output+"/A2_trimming_and_post_trimming_quality/"+sample+"_fastp_trimmed_2P.fastq.gz")
    
    #run unicycler_post_fastp virus reads
    #unicycler_virus(output, sample, output+"/B_mappings_and_new_reads/"+sample+"_virus_R1.fq.gz", output+"/B_mappings_and_new_reads/"+sample+"_virus_R2.fq.gz")
        
    #run ragtag_post_unicycler
    ragtag_both(output, sample, referencev, "assembly.fasta")
  
    #run medusa_post_ragtag
    medusa_both(output, sample)
  
    #run quast_post_fastp raw reads
    quast_raw(output, sample, output+"/C3_denovo3_medusa/medusa.final.fasta", output+"/C2_denovo2_ragtag/ragtag.scaffold.fasta", output+"/A2_trimming_and_post_trimming_quality/"+sample+"_fastp_trimmed_1P.fastq.gz", output+"/A2_trimming_and_post_trimming_quality/"+sample+"_fastp_trimmed_2P.fastq.gz")
  
    #run quast_post_fastp virus reads
    #quast_virus(output, sample, output+"/C3_denovo3_medusa/medusa.final.fasta", output+"/C2_denovo2_ragtag/ragtag.scaffold.fasta", output+"/B_mappings_and_new_reads/"+sample+"_virus_R1.fq.gz", output+"/B_mappings_and_new_reads/"+sample+"_virus_R2.fq.gz")
        
    #run qualimap of raw bwa alignment
    qualimap_raw(output, sample)

    #run qualimap of virus bwa alignment
    #qualimap_virus(output, sample)

    #run variant calling with VARSCAN2 on raw bwa alignment
    varscan2_raw(output, sample, output + "/B_mappings_and_new_reads/" + sample + "_raw.sort.bam", referencev, threads)

    #run variant calling with VARSCAN2 on virus bwa alignment
    #varscan2_virus(output, sample, output + "/B_mappings_and_new_reads/" + sample + "_virus.sort.bam", referencev, threads)
     
    #End timer
    end_time = time.time()
    execution_time = end_time - start_time
    executiont_minutes = execution_time/60
    print("\nSample execution time: ",executiont_minutes, "minutes")
    print("\n####################################################################################")
    print("##                      THANK YOU FOR USING 'PePApipe'!!                          ##")
    print("####################################################################################\n")
 
if __name__ == "__main__":

    main()