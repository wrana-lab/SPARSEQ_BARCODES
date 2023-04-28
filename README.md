# SPARSEQ_BARCODES
Sub-version of SPARseq pipeline designed to handle specialized paired-end data with 5bp barcodes for tracking well-to-well spread of reads.  Outputs the same files as the normal SPARseq pipeline with additional outputs related to the barcode contents per sample. 


Part 1: Barcode analysis

Requirements:

OSX or Linux system

Python3.9 [https://www.python.org/]

Bowtie [http://bowtie-bio.sourceforge.net/index.shtml]

HTSeq [https://htseq.readthedocs.io/en/master/]

Cutadapt [https://cutadapt.readthedocs.io/en/stable/]

Descriptions

sparseq_pipelineNS_V5.sh

Wrapper script to run cutadapt for trimming of 5bp barcode from R1 files, then run normal processing of trimmed SPARseq R1 files with sparseq_analysisNS_V5.1.py. 

sparseq_analysisNS_V5.1.py

Script for analyzing top reads of R1 files after the 5bp barcode was trimmed. Basically the same as our normal sparseq script but with a few numbers adjusted for the Srbdv2 barcodes. This allows us to ensure consistency with non-barcode data.

sparseq_pipelineNS_V6.sh

Wrapper script to do trimming of non-barcode bases from R1 and R2 fastq files. It uses cutadapt to trim last 139bp of R1 and last 3bp of 8bp R2. It executes sparseq_analysisNS_V6.1.py and sparseq_analysisNS_V7.1.py which reformat the bowtie outputs and check barcode pairs to measure well identity.
 
sparseq_analysisNS_V6.1.py

Script to reformat outputs from R1 and R2 bowtie mapping from the wrapper script. 

sparseq_analysisNS_V7.1.py

Script to process pairs of reads and match up R1 and R2 barcodes to check well identity. 





Part 2: sequence analysis of barcoded data.

Extension of SPARSeq barcode pipeline to detect possible quasispecies with paired end reads.

Requirements:

OSX or Linux system

Python3.9 [https://www.python.org/]

Clustal Omega running locally [http://www.clustal.org/omega/]

R v4.2.2 [https://cran.r-project.org/]

Workflow

quasispecies_wrapper_v3_barcode.sh
Wrapper script for the QS analysis pipeline. Unzips fastq files, writes list of filenames, executes sequence pileup script, executes clustal omega locally for alignment and conversion to aligned fasta, executes fasta cut/trim script, then re-zips fastq files and reports run time. 

sparseq_srbdv2.py
Initial processing script for raw fastq files. It processes pairs of fastq files by checking for matching read IDs, checking that the R1 and R2 barcodes match the sample's well in the 384w plate, then selecting, binning, and counting Srbdv2 sequences. Next it uses the kneed library to fit a knee plot model to determine an appropriate count cutoff value for the sample, which is used to filter the sequences within the sample. If no knee point is found it uses a count cutoff value of 0.5% (based on total Srbdv2 counts). Sequences passing the cutoffs are written to a .fa file.

Clustal Omega is used locally to align sequences that passed the cutoffs and convert to an aligned .fasta file.

fasta_cutter_barcode_v3.py

This script takes the aligned .fasta file and organizes the sequences, and outputs a file that has the base for each position for each sequence. This output is used later to determine the contents of each position in the Srbdv2 sequence.

quasispecies_barcode.R
This script takes aligned sequences and translates them, outputting the sequences with a reference sequence for easy viewing in snapgene.

Run final output from r through clustal omega locally again to get finalized alignments of top aggregated sequences. This can be opened in snapGene for visualization.
