# PePApipe-pipeline
Complete bioinformatics analysis pipeline for African Swine Fever virus genomes  

## Description
Despite the existence of published generic tools designed to carry out bioinformatics analyses of genomes from different origins, the differences between the genomes may advice to use specifically designed tools, especially in the case of viruses, such as the African Swine Fever Virus (ASFV). 

The problem with viruses in general, and with the ASFV in particular, is that their genome is highly variable, and similarly to what occurs with influenza viruses, it is not unfrequent for these genomes to undergo point mutations of small indels (insertions and deletions) although other times these changes can be large homologous recombinations. An additional peculiarity of ASFV genomes is that both left and right terminus regions of its DNA are very variable, often containing inverted terminal repeats which further complicate both the assembly and annotation of these genomes. 
It is therefore very important to genetically characterize these changes if we want to have an accurate representation of the true structure of the genome. 

Although several tools exist than can be used to perform assembly and analysis of DNA genomes of viruses, such as ASFV, the authors of this work have not found to date a single tool able to compile all steps involved in this analysis in a complete, simple and accessible way. ‘PePApipe’ is fit for this purpose, and implements the work in a rapid, complete, efficient and clear manner. 

This pipeline has been developed bearing in mind that end-users are mostly laboratory people with limited bioinformatics skills.


## Pipeline Workflow (Diagram)
This pipeline, designed and programmed in Python, can be run using basic Linux commands to launch a Bash script with a Slurm protocol and can be executed on multi-sample batches. The work is sequentially performed over three main areas of genomic analysis: quality control and pre-processing of raw reads, denovo genome assembly and multiple variant calling (Figure 1). Moreover, it may be run in a locally or remotely by accessing a computing server. A basic knowledge of Linux programming language is necessary both to install the programmes and to run the pipeline.

Starting from raw data (paired .FASTQ files) obtained from short-read sequencing platforms (Illumina), ‘PePApipe’ implements in a sequential manner all crucial steps by using 14 software tools adequately built into a single pipeline.


![image](https://github.com/user-attachments/assets/f84b44a6-66a9-473c-9bd6-0ece897f9146)

Figure 1. Overall bioinformatic analysis flow followed with PePApipe (in-house built python pipeline specifically designed for ASF viruses). 

## Installation Steps


## Example Usage


## Expected Outputs


## License & Citation

