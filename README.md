# SPARSEQ_BARCODES


This is the pipeline for testing barcoded clinical samples to assess barcode contamination across wells. It is a sub-version of SPARseq pipeline designed to handle specialized paired-end data with 5bp barcodes for tracking well-to-well spread of reads. Outputs the same files as the normal SPARseq pipeline with additional outputs related to the barcode contents per sample.


Requirements:
Python3 [https://www.python.org/] [os, argparse, re, openpyxl, pandas, datetime]
Bowtie 0.12.7 [http://bowtie-bio.sourceforge.net/index.shtml]
HTSeq 0.13.5 [https://htseq.readthedocs.io/en/master/]
Cutadapt 3.4 [https://cutadapt.readthedocs.io/en/stable/]
Tested on MacOS 12.7.5

The run folder and resource folder must be subfolders of the same parent folder:

Run folder structure

	/runclXX
		-/R1andR2_files
			-fastq files

Resource folder structure

	/resources/
		-/bowtie_index ** this is the same set of index files as in the main sparseq pipeline.
				-bowtie index files
				-reference fa file
				-reference GTF file
		-/bowtie_index_R1_barcode 
				-bowtie R1 index files
				-reference R1 fa file
				-reference R1 GTF file
		-/bowtie_index_R2_barcode
				-bowtie R2 index files
				-reference R2 fa file
				-reference R2 GTF file
		-codon_aa_table.csv
		-BCrowmatch.csv
		-BCcolmatch.csv
		-sparseq_pipelineNS_V5.sh
		-sparseq_pipelineNS_V6.sh
		-sparseq_analysisNS_V5.1.py
		-sparseq_analysisNS_V6.1.py
		-sparseq_analysisNS_V7.1.py

Explanations of files:

bowtie_index contains the same index files as the main clinical sparseq pipeline. The R1 and R2 barcode references are based on the row and column barcodes, and are used to quantify the barcodes. 

sparseq_pipelineNS_V5.sh: Wrapper script to run cutadapt for trimming of 5bp barcode from R1 files, then run normal processing of trimmed SPARseq R1 files with sparseq_analysisNS_V5.1.py.

sparseq_analysisNS_V5.1.py: Script for analyzing top reads of R1 files after the 5bp barcode was trimmed. Similar to our normal sparseq script but with a few numbers adjusted for the barcodes. This allows us to ensure consistency with non-barcode data.

sparseq_pipelineNS_V6.sh: Wrapper script to do trimming of non-barcode bases from R1 and R2 fastq files. It uses cutadapt to trim last 139bp of R1 and last 3bp of 8bp R2. It executes sparseq_analysisNS_V6.1.py and sparseq_analysisNS_V7.1.py which reformat the bowtie outputs and check barcode pairs to measure well identity.

sparseq_analysisNS_V6.1.py: Script to reformat outputs from R1 and R2 bowtie mapping from the wrapper script.

sparseq_analysisNS_V7.1.py: Script to process pairs of reads and match up R1 and R2 barcodes to check well identity.

codon_aa_table.csv contains nucleotide to amino acid conversion info.

BCrowmatch.csv and BCcolmatch.csv contain tables of barcodes used per row and column.



Initial one-time steps:

Build normal clinical bowtie index from within the bowtie_index folder with 
```
bowtie-build clinical_V2_sparsq_V1p2_amplicons.fa clinical_V2_sparsq_V1p2_amplicons
```
Build R1 barcode bowtie index from within the bowtie_index_R1_barcode folder with
```
bowtie-build sparsq_R1_barcode_amplicons.fa sparsq_R1_barcode_amplicons
```
Build R2 barcode bowtie index from within the bowtie_index_R2_barcode folder with
```
Bowtie-build sparsq_R2_barcode_amplicons.fa sparsq_R2_barcode_amplicon
```


Instructions for one run:
```
1. Create run folder runclX and subdirectory R1andR2_files; place fastq files into R1andR2_files
2. Run $bash sparseq_pipelineNS_V5.sh path/to/run/folder/runclXX path/to/resource/folder
		> Wrapper script performs the following: 
			>Copies R1 fastq files to R1_files directory
			>unzips fastqs
			>runs cutadapt to trim barcodes
			>runs bowtie and htseq as in clinical pipeline
			>runs sparseq_analysisNS_V5.1.py for usual clinical processing
				>Creates samcounts directory then runs bowtie and HTseq
				>Creates results directory then run Clinical_V2_sparseq_analysis_V1.3.2 which outputs several key results tables which are compared to the outputs from the matching clinical run to confirm consistency
					-runclXX_AllVariantDetails.xlsx
					-runclXX_CountTable.xlsx
					-runclXX_SelectedVariants.xlsx
					-sparseq_report_runclXX.xlsx
					-please note that as in the clinical pipeline variant calling is manually completed by a supervisor based on review of the 4 outputs 
3. Run $bash sparseq_pipelineNS_V6.sh path/to/run/folder/runclXX path/to/resource/folder
		> Wrapper script performs the following:
			>unzips R1 and R2 fastqs
			>runs cut adapt to trim everything except the barcodes
			>bowtie and HTSeq to quantify the R1 and R2 barcodes per file
			>runs sparseq_analysisNS_V6.1.py
			>runs sparseq_analysisNS_V7.1.py
			*****In this portion of the analysis some reads can be skipped with the warning "read pair error" which is normal.
			>Key outputs from this step include
				>runclXX_R2BC_Table_Reordered.xlsx which has counts of R2 barcodes per sample/well
				>runclXX_R1BC_Table.xlsx which has counts of R1 barcodes per sample/well
				>runclXX_BC_PairedAnalysis.xlsx which has combined counts of R1 and R2 barcode pairs per sample
```



