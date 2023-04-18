#/bin/bash

###This is the wrapper for the tiny UMI versions of sparseq.

# Dependencies:
#	1. Bowtie 1 (http://bowtie-bio.sourceforge.net/index.shtml)
#	2. HTSeq* (https://htseq.readthedocs.io/en/release_0.11.1/index.html)
#	3. Python3 [os, argparse, re, openpyxl, pandas, datetime]*
#		a. sparseq_analysis_v1.3.X.py
#	*To install pkgs
# 		$ module load conda
#		$ pip install --user [package]

run_folder=$1 #./runclX/
resource_folder=$2 #./resources/SPAR_SEQ/

#module load conda
echo "***** Beginning SPAR-Seq Pipeline "$(date)"*****"

echo "	*Unzipping FASTQ files...*"
gunzip $run_folder/R1andR2_files/*.gz

##execute cutadapt to chop off LAST 139 bp of R1
#and last 3 of 8bp R2
#this is for NextSeq


echo "  *Running cutadapt on R1s..."

for f in $run_folder/R1andR2_files/*_R1_001.fastq; do cutadapt -u -139 -o $f.trimmed.R1.fastq $f; done
##execute cutadapt to chop off LAST 1 bp of R2
echo "  *Running cutadapt on R2s..."
for f in $run_folder/R1andR2_files/*_R2_001.fastq; do cutadapt -u -3 -o $f.trimmed.R2.fastq $f; done

echo "	*Running bowtie alignment and counting reads with HTseq for R1 UMIs...*"
mkdir -p $run_folder/R1andR2_files/samcounts

#may need to verify if *R1*.sam and *R2*.sam will work
# -v 0 for 0 mismatch; -k1 for top 1 alignment, -m 1 to suppress read if more than 1 alignment exists
for f in $run_folder/R1andR2_files/*.trimmed.R1.fastq; do bowtie $resource_folder/bowtie_index_R1_UMI/sparsq_R1_UMI_amplicons $f $f.sam --best --norc -v 0 -k 1 -m 1 -S; done && for x in $run_folder/R1andR2_files/*trimmed.R1.fastq.sam; do python3 -m HTSeq.scripts.count -f sam -t CDS $x $resource_folder/bowtie_index_R1_UMI/sparsq_R1_UMI_amplicons.gtf > $x.R1count.txt; done;
echo "	*Running bowtie alignment and counting reads with HTseq for R2 UMIs...*"
#for f in $run_folder/R1andR2_files/*.trimmed.R2.fastq; do bowtie $resource_folder/bowtie_index_R2_UMI/sparsq_R2_UMI_amplicons $f $f.sam --best -v 0 -k 1 -m 1 -S; done && for x in $run_folder/R1andR2_files/*R2*.sam; do python3 -m HTSeq.scripts.count -f sam -t CDS $x $resource_folder/bowtie_index_R2_UMI/sparsq_R2_UMI_amplicons.gtf > $x.R2count.txt; done;
for f in $run_folder/R1andR2_files/*.trimmed.R2.fastq; do bowtie $resource_folder/bowtie_index_R2_UMI/sparsq_R2_UMI_amplicons $f $f.sam --best --norc -v 0 -k 1 -m 1 -S; done && for x in $run_folder/R1andR2_files/*trimmed.R2.fastq.sam; do python3 -m HTSeq.scripts.count -f sam -t CDS -s reverse $x $resource_folder/bowtie_index_R2_UMI/sparsq_R2_UMI_amplicons.gtf > $x.R2count.txt; done;


mv $run_folder/R1andR2_files/*count.txt $run_folder/R1andR2_files/samcounts

echo "	*Running analysis of mapped UMIs...*"
mkdir -p $run_folder/results
python3 $resource_folder/sparseq_analysisNS_V6.1.py $run_folder $resource_folder

echo "	*Remove trimmed fastq files*"
rm $run_folder/R1andR2_files/*trimmed.R1.fastq
rm $run_folder/R1andR2_files/*trimmed.R2.fastq

#run paired analysis of fastq files
echo "  #Running paired fastq UMI analysis...*"
python3 $resource_folder/sparseq_analysisNS_V7.1.py $run_folder $resource_folder

echo "	*Zipping FASTQ files...*"
gzip $run_folder/R1andR2_files/*.fastq

echo "	*Remove sam files*"
rm $run_folder/R1andR2_files/*.sam

echo "*****"$(date)" SPAR-Seq UMI Pipeline is now complete. *****"
