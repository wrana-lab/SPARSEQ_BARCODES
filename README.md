# SPARSEQ_BARCODES
Sub-version of SPARseq pipeline designed to handle specialized paired-end data with 5bp barcodes for tracking well-to-well spread of reads.  Outputs the same files as the normal SPARseq pipeline with additional outputs related to the barcode contents per sample. 


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

