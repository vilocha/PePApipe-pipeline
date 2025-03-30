# PePApipe-pipeline
Complete bioinformatics analysis pipeline for African Swine Fever Virus (ASFV) genomes  

## Description
Despite the existence of published generic tools designed to carry out bioinformatics analyses of genomes from different origins, the differences between the genomes may advice to use specifically designed 
tools, especially in the case of viruses, such as the African Swine Fever Virus (ASFV). 

The problem with viruses in general, and with the ASFV in particular, is that their genome is highly variable, and similarly to what occurs with influenza viruses, it is not unfrequent for these genomes to undergo point mutations of small indels (insertions and deletions) although other times these changes can be large homologous recombinations. An additional peculiarity of ASFV genomes is that both left and right terminus regions of its DNA are very variable, often containing inverted terminal repeats which further complicate both the assembly and annotation of these genomes. 
It is therefore very important to genetically characterize these changes if we want to have an accurate representation of the true structure of the genome. 

Although several tools exist than can be used to perform assembly and analysis of DNA genomes of viruses, such as ASFV, the authors of this work have not found to date a single tool able to compile all steps involved in this analysis in a complete, simple and accessible way. ‘PePApipe’ is fit for this purpose, and implements the work in a rapid, complete, efficient and clear manner. 

This pipeline has been developed bearing in mind that end-users are mostly laboratory people with limited bioinformatics skills.


## Pipeline Workflow (Diagram)
This pipeline, designed and programmed in Python, can be run using basic Linux commands to launch a Bash script with a Slurm protocol and can be executed on multi-sample batches. The work is sequentially performed over three main areas of genomic analysis: quality control and pre-processing of raw reads, denovo genome assembly and multiple variant calling (Figure 1). Moreover, it may be run in a locally or remotely by accessing a computing server. A basic knowledge of Linux programming language is necessary both to install the programmes and to run the pipeline.

Starting from raw data (paired .FASTQ files) obtained from short-read sequencing platforms (Illumina), ‘PePApipe’ implements in a sequential manner all crucial steps by using 14 software tools adequately built into a single pipeline.

PIPELINE DIAGRAM
![image](https://github.com/user-attachments/assets/f84b44a6-66a9-473c-9bd6-0ece897f9146)

Figure 1. Overall bioinformatic analysis flow followed with PePApipe (in-house built Python pipeline specifically designed for ASF viruses). 

## Installation Steps
The following steps must be implemented in the order suggested to successfully install all necessary programmes to be able to run this pipeline. The actual code lines used within the pipeline are clearly shown in the pipeline script itself, which can be downloaded from this page. These code lines are already built within the pipeline and will be executed at once when the pipeline is run.

According to needs, the different modules in the pipeline can be activated/ deactivated by just uncommenting/ commenting them in the main section of the pipeline script, respectively.

The input or starting point of the pipeline is always the set of raw reads which are the outputs of the sequencing platform (short reads from Illumina, although long reads from Oxford Nanopore Technologies can also be used) organized in two files: .FASTQ R1, or forward reads, and .FASTQ R2, or reverse reads. The output files of the different steps are specified on the installation webpages.


STEPS:

QUALITY CONTROL OF READS

This is the very first step which is necessary to ensure a good quality of the two .FASTQ files provided by the sequencing platform. The quality check of raw reads is carried out with the tool FastQC. As a general rule, the per-base sequence quality (Phred score) should always be ideally over 28, and under any circumstances it should be below 20.

Installation: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Once the quality control step is complete and one is satisfied that the quality of sequences is adequate to continue the processing, the reads are ‘cleaned’ with the tool fastp. This step is necessary to ensure that adapters, low quality reads and duplicate reads are removed off the .FASTQ files. In addition, fast produces a post-trimming quality analysis report. 
The adapters were included during the sequencing process but are not needed for the bioinformatics analysis. Only reads with a minimum phred quality score of 20 and a minimum length of 15 base pairs are kept after applying fastp, but these parameters can be changed in the code line according to user needs. 
A further quality control using FastQC can be performed for comparison purposes although this is not necessary because fastp already provides a post-trimming quality report, as mentioned above.

Installation: https://github.com/OpenGene/fastp 


VIRAL READS CONTENT CHECK

In this step the reads that passed the previous filters (‘clean reads’) are mapped against the virus reference genome in .FASTA format chosen for convenience in each particular case using the software tools BWA-mem2 and Samtools. 
This step is probably the most important of all steps included in this part of the flow, and of the pipeline for that matter. The viral genome reference used for convenience by the authors of the pipeline was ASFV Georgia 2007/1 (Genbank assembly: GCA_003047755.2).
The objective of this mapping is two-fold: on the one hand it allows to quantify the amount of reads within the original .FASTQ files actually mapping with an ASFV (and thus to also quantify the amount of reads not mapping with an ASFV), and on the other hand, to create a set of reference-mapped contigs (.BAM file) as a first step in the search of genetic variants (SNPs and indels) in the genome of our virus problem compared to the genome of the virus reference used.

Installation: https://github.com/bwa-mem2/bwa-mem2; https://www.htslib.org/download/ 


QUALITY CHECK OF READS MAPPED TO REFERENCE GENOME OF CHOICE

One of the tools most commonly used for this purpose is Qualimap 2. The rationale behind is to assess the quality of the mapping of trimmed raw reads to virus reference genome which was performed to calculate the amount of viral reads present. The outputs can be verified in .PDF and .HTML formats. 

Installation: http://qualimap.conesalab.org/ 


QUANTITATIVE ASSESSMENT OF VIRAL READS (THEORETICAL CUT-OFF POINT: 75%)

It is important to carry out all bioinformatics analyses only on viral reads free of contamination with reads from other origins (host, environment or handling personnel) or, if this is not possible, at least on a as high as possible number of viral reads. This is because the results we may obtain if the analyses are performed on a high number of non-viral reads may be spurious and therefore not representative of the virus we are trying to characterize.
This is a point in the pipeline where a decision has to be made as to the quality level we want to apply to our analysis in relation to the available reads (Figure 1). A theoretical threshold (cut-off point) of 75% has been established by the authors in order to decide whether to continue the analyses with the original .FASTQ files already filtered (‘clean reads’) or whether to carry out yet another filtering step to get rid of non-viral reads. If the latter is decided, a new set of .FASTQ files must be created by mapping the ‘clean reads’ against all published reliable reads belonging to the Superkingdom Viruses (NCBI: txid10239).
Apart from the methods used to isolate and characterize the viral DNA, the method used to prepare the libraries (with or without capture) is the single most influential step on the amount of foreign (non-viral) DNA found in our sample. If a non-capture library preparation method is used, there will be a high number of non-viral reads present in the .FASTQ files, thus these may be larger in size although not necessarily better in quality.
Hence, this decision threshold is arbitrary and may vary depending on the initial .FASTQ files and whether the library has been prepared using capture or not. Libraries prepared without capture are likely to contain percentages of viral reads lower than 75%.      


REMOVAL OF READS NOT BELONGING TO VIRUSES

To carry out the removal of non-viral reads a loop with an additional Python script using the tool Kraken2 is necessary. This extra step will create a new set of .FASTQ files containing only viral reads. 
As already mentioned, this step can be performed whenever the chosen threshold of viral reads is not reached (in our case 75%) or it may also be systematically run every time, independently of the amount of non-viral reads present in the original ‘cleaned’ set of .FASTQ files. 
If there is availability of a substantial amount of samples, a benchmark exercise can be carried out to refine the decision threshold or cut-off percentage point, making it lower or higher than 75% according to each specific situation.

Installation:

The programme used for this purpose is called Kraken2 and the procedure consists of two steps: building of the virus database and extraction of the viral reads. Instructions about these procedures are detailed on this website: https://ccb.jhu.edu/software/krakentools/index.shtml?t=extractreads. For this step, purposely-built BASH scripts have been designed and are downloadable through this webpage.


DE NOVO VIRAL ASSEMBLY

This step is mainly performed with the software tool Unicycler, and ultimately refined with the programmes RagTag (using Unicycler outputs as inputs) and MeDuSa (using RagTag outputs as inputs). 
This suite uses several dependency programmes in order to obtain a full viral de novo genome assembly starting from the .FASTQ files containing our already ‘cleaned’ reads, be these the original .FASTQ files or the new set created using Kraken2. 
Although Unicycler was originally built for bacteria, ideally combining both short (Illumina) and long reads (Oxford Nanopore Technologies), it performs equally well using only viral genomes short reads. 
A most positive aspect about this type of assembly is that a viral genome reference is not used at all, thus no bias from this source can be introduced, allowing for a genuine assembly of the putative whole genome of our virus problem. This step is the most important within the de novo assembly part of the flow, but unfortunately it is also the most consuming in terms of time and resources.      
There are particularities in relation to the three tools used in the de novo assembly. In order to obtain a proper assembly of our virus whole genome it is important to specify that what we are looking for is a linear sequence, in other words, that the expected number of linear (non-circular) sequences in the underlying sequence being assembled is just one, and so these requirements must be included in the code in the Unicycler section.

Installation: https://github.com/rrwick/Unicycler 

After the first assembly step with Unicycler is concluded, a correction of the obtained contigs converting these into longer contigs is performed with the RagTag tool. This step is necessary to guide the de novo assembly performed by Unicycler into plausible contigs belonging to the ASFV virus, using for that purpose the chosen viral genome reference in .FASTA format (in our case ASFV Georgia 2007/1; Genbank assembly: GCA_003047755.2). 
It is important to note that even though a viral reference is used, not a single fragment of it is incorporated into the assembly obtained with Unicycler, so the only activity performed in this step is that of contig correction and scaffolding. RagTag comprises four software tools (correct, scaffold, patch and merge) but only the first two, correct and scaffold, are used in the pipeline.      

Installation: https://github.com/malonge/RagTag

After the correction and scaffolding performed by RagTag, a third tool called MeDuSa is incorporated as part of the de novo assembly process. MeDuSa performs a scaffold correction in a similar way to RagTag, but instead of using only one viral reference it uses as many well curated references as we are able to provide. Up to 10 references belonging to ASFV virus Genotype II have been included in this pipeline, which in alphabetical order are: Arm_07_CBM_c2_LR812933.1, ASFV_Georgia_2007_1_FR682468.2.fasta, ASFV_HU_2018_MN715134.1, ASFV_LT14_1490_MK628478.fasta, Belgium_2018_1_LR536725.1, Estonia_2014_LS478113.1, Korea_YC1_2019_ON075797.1, POL_2015_Podlaskie_MH681419.1, Tanzania_Rukwa_2017_1_LR813622.1 and Ulyanovsk_19_WB_5699_MW306192.1  
The higher the number of references used, the better the result of this step will be, always making sure the references are akin to the ASFV Genotype we are studying. Briefly, the contigs obtained from the previous step are further assembled into larger scaffolds with MeDuSa. If the quality of the original reads is high this step generally results in one final contig, which is the putative consensus sequence of the virus under study. On the contrary, poor quality sequences may result in a large number of contigs, and if the number of contigs exceeds 3-5 it would be advisable to reject those sequences.

Installation: https://github.com/combogenomics/medusa 


QUALITY ASSESSMENT OF THE DE NOVO ASSEMBLY

This is a quality analysis that can be performed over all steps in the de novo assembly loop with the intention of investigating the overall quality of the process. The tool included in the pipeline for this task is Quast (QUality ASsessment Tool) and the outputs are available on different formats: .TXT, .PDF, .TSV, .HTML, etc.

Installation: https://github.com/ablab/quast 


ANNOTATION OF THE CONSENSUS SEQUENCE WITH REFERENCE GENOME OF CHOICE

This is not an essential step in the pipeline but rather a complementary analysis that can be carried out alongside with it. The annotation of the de novo consensus sequence is an important task to identify and tag the different genes which have been assembled together. The best tool released so far for this purpose is a universal tool for general annotation that works with viruses, called GATU (General Annotation Transfer Utility), although many viral genes remain unannotated after its use. Two files are necessary to carry out the annotation: a reference sequence in .GBK format and the problem sequence in .FASTA or .GBK formats. 
The ASFV poses a problem that may render the use of this tool inappropriate, and that is the high genome variability of these viruses. Because of this, a full manual curation exercise is always necessary to arrange the genes according to references in a left to right terminus direction (5’ end to 3’ end) and to clean off all hypothetical genes whose probabilities of belonging to our virus problem are very low (the so called ‘false predictions’). The output of GATU is on a .TSV format and this can be easily transformed into Excel files.    

Installation: https://4virology.net/virology-ca-tools/gatu/ 


VISUALIZATION OF MAPPED READS AGAINST VIRUS REFERENCE GENOME OF CHOICE

This is also an optional step that is not an integral part of the pipeline. The software tool most universally used for visualizing genome mappings is IGV (Integrative Genomics Viewer). 
Two genomes are necessary to feed the tool: a reference genome in .FASTA format and the genome under study, which is always provided as a .BAM file accompanied by its index (.BAI file).   

Installation: https://igv.org/ 


SEARCH FOR VIRAL GENOME VARIANTS (VARIANT CALLING)

Several tools are used to reveal genetic variants (SNPs and indels) present in the virus genome under study when compared to a viral reference genome. 
In general, there are two main types of tools: those that use a Bayesian inference and those that use a deterministic/numerical approach. The tool included in this pipeline is called VarScan 2 and is one of the second type (numerical approach). 

In order to obtain the most plausible record of variants with this tool the following five parameters have been defined with the specified values:

-	Minimum number of reads with change to call a variant: 5
-	Minimum depth (coverage=number of reads in total supporting that variant): 20
-	P-value: 0.001
-	Minimum average quality: 30
-	Minimum variant frequency in that position: 0.5
-	Ploidy: 1 (haploid)

This tool may be complemented with other tools in order to obtain more robust/reliable results. We complemented VarScan 2 with Qiagen CLC Genomics Workbench (https://digitalinsights.qiagen.com/products-overview/discovery-insights-portfolio/analysis-and-visualization/qiagen-clc-genomics-workbench/). 

Installation: https://varscan.sourceforge.net/ 


OVERALL QUALITY ASSESSMENT OF THE ANALYSIS

We performed an overall quality assessment of the analysis using the tool MultiQC. This tool evaluates all inputs and outputs in the working directory and performs a general summary assessment in the same folder where all folders used in the analysis are located.  

Installation: https://github.com/MultiQC/MultiQC


## Example Usage

Once the installation of the tools is complete it has to be decided whether to analyse one sample or a batch of samples. The pipeline (File 1), as a Python extension file, is then run using either the script for a single sample (File 2) or the multisample script (File 3), both in bash. If we decide to analyse a batch of samples, a text file with the list of samples (File 4) must also be provided. 

If we opt for generating a new set of .FASTQ files that only contain viral reads, there is an additional Python extension file (File 5) which is launched with a script in bash (File 6), but first we must build the Kraken2 database with yet another script in bash (File 7).

Finally, the last bash script is needed to carry out an overall analysis using MultiQC (File 8). 

This pipeline was originally built and designed to be launched using a slurm protocol. The launching command is either 'sbatch run_pipe_single.sh' or 'sbatch run_pipe_batch.sh' and if we want to follow the execution on the screen the following command must be typed on the same terminal (do not open a new terminal for this): 'tail -F pipeline_log.txt'. After typing this line, an immediate message will appear telling us that "the pipeline_log.txt has not been found", but once the job enters a node and the execution commences the message will be that "the pipeline_log.txt was found" and from that moment we can start seeing the progress on the terminal screen.

In summary, all files used, their running codes and the orden of execution are listed below:

A) Pipeline in Python language 'PePApipe.py' (File 1) to be run on a single sample with bash script 'run_pipe_single.sh' (File 2) or on multiple samples with bash script 'run_pipe_batch.sh' (File 3) in which case a 'samples.txt' file (File 4) must be provided

Code lines in script 'run_pipe_single.sh':

unset DISPLAY   #This command is necessary to avoid any files from any of the tools from opening during execution of the pipeline, which will abort the pipeline execution

module load ... (branches, modules and versions to be loaded in each particular server)

python -u ./PePApipe.py -s SAMPLE -o ./SAMPLE -t 12 -r1 ./SAMPLE/*R1*fastq.gz -r2 ./SAMPLE/*R2*fastq.gz -R ../virus_reference/virus_reference.fasta &> pipeline_log.txt

Code lines in script 'run_pipe_batch.sh':

unset DISPLAY   #This command is necessary to avoid any files from any of the tools from opening during execution of the pipeline, which will abort the pipeline execution

module load ... (branches, modules and versions to be loaded in each particular server)

while read -r line; do
python -u ./PePApipe.py -s "$line" -o "./$line" -t 12 -r1 "./$line/*R1*fastq.gz" -r2 "./$line/*R2*fastq.gz" -R ../virus_reference/virus_reference.fasta &>> pipeline_log.txt
done < ./samples.txt 

B) Building of Kraken2 database with script 'run_krakenDB_build' (File7) 



C) Python extension file 'extract_kraken_reads.py' for extracting the viral reads (File 5) using script 'run_extract_kraken_fqs' (File 6)     

D) 

The eigth files that can be downloaded and adapted to each particular case to successfully run PePApipe are thus:

NB: If we want to run the pipeline locally, the running codes which are embedded within the bash scripts must be typed straigth into the terminal, by-passing the slurm loop. 

In addition, please note that in the analysis working directory where we want to execute the pipeline there must be the eigth files mentioned above (two Python files, five bash scripts and a text file with the list of samples) and as many folders named with as many different samples we have. Inside each of the sample folders we must have the two reads corresponding to that sample, Forward and Reverse, either from Illumina (short reads) or Oxford Nanopore Technologies (long reads). There will be eight output folders (see Section of Expected Outputs below) which will be created within each of the sample folders.

Finally, the two sets of genome references, which can be called 'virus_reference' (needed for bwa-mem2) and 'other_references' (needed for MeDuSa) must be located one hierarchy above of the directory where the sample folders and scripts are located, i.e., at the same level as our analysis working directory. The references for MeDuSa can be used by the tool direct from .FASTA format files but the viral genome reference used by bwa-mem2 must be indexed first in order to be available for the tool with this command: bwa index virus_reference.fasta.    

## Expected Outputs

The pipeline produces all files necessary to interpret the results, classifying them into folders, in an approximate time which depends on the size of the .FASTQ files. The time used may span from 10 minutes with files of around 0.5 Mb through to 30 minutes with files of 1 Mb in size or to longer times if high capacity sequencers are used. Several intermediate files are also created in the process, which can be used as inputs for further or parallel analyses, as well as to check for possible errors during the analysis. Among these complementary post-analysis steps are viral annotation (GATU) or genome visualization (IGV).  

All steps are amenable to user control by means of switching on/off the necessary sections in each particular case. There are several steps where the quality of the processes can be checked to make sure the progress of the results is in line with our expectations. 
The pipeline can be easily adapted to viruses other than ASFV by changing the parameters relevant to the new virus and running the specific pipeline sections accordingly.

The eigth folders that are created when the pipeline has finished running along with the information retrievable for the interpretation of results are listed below:


The most important outputs are the consensus sequence of the genome of the virus we are studying and the table of variants for that particular virus against the genome reference chosen in the first place. Secondary outputs are reference mapping parameters, de novo assembly parameters and quality assessment summaries. As already mentioned, the findings obtained may be further investigated using IGV and GATU.

## License & Citation

LICENSE

This pipeline is publicly licensed to be used under the GNU GENERAL PUBLIC LICENSE Version 3, 29th June 2007

CITATION

Once published the title and author list of this work will be as follow:

Title:
'PePApipe': a complete bioinformatics analysis pipeline for African Swine Fever Virus genomes

Authors:
Vicente Lopez Chavarrias 1, Irene Aldea Ramos 1, Jovita Fernandez Pinero 1 
1.	Centro de Investigación en Sanidad Animal (CISA), Instituto Nacional de Investigación y Tecnología Agraria y Alimentaria (INIA), Consejo Superior de Investigaciones Científicas (CSIC), Valdeolmos, Madrid, Spain
