#/bin/bash

# Dependencies:
#	1. Bowtie 1 (http://bowtie-bio.sourceforge.net/index.shtml)
#	2. HTSeq* (https://htseq.readthedocs.io/en/release_0.11.1/index.html)
#	3. Python3 [os, argparse, re, openpyxl, pandas, datetime]*
#		a. sparseq_analysis_v1.3.X.py
#	*To install pkgs
# 		$ module load conda
#		$ pip install --user [package]


#NextSeq read cycles can't be less than index cycles so R1 has 144 and R2 has 8 bp
run_folder=$1 #./runclX/
resource_folder=$2 #./resources/SPAR_SEQ/

#module load conda
echo "***** Beginning SPAR-Seq Pipeline "$(date)"*****"

#copy R1s only to R1_files
echo "  *Copying R1 files...*"
mkdir -p $run_folder/R1_files
cp $run_folder/R1andR2_files/*_R1_001.fastq.gz $run_folder/R1_files/

echo "	*Unzipping FASTQ files...*"
gunzip $run_folder/R1_files/*.gz

##execute cutadapt to chop of first 5 bp of R1
echo "  *Running cutadapt..."
for f in $run_folder/R1_files/*.fastq; do cutadapt -u 5 -o $f.trimmed.fastq $f; done

echo "	*Running bowtie alignment and counting reads with HTseq...*"
mkdir -p $run_folder/R1_files/samcounts
for f in $run_folder/R1_files/*.trimmed.fastq; do bowtie $resource_folder/bowtie_index_V1p2/sparsq_V1p2_amplicons $f $f.sam --best -v 3 -k 1 -m 1 -S; done && for x in $run_folder/R1_files/*.sam; do python3 -m HTSeq.scripts.count -f sam -t CDS $x $resource_folder/bowtie_index_V1p2/sparsq_V1p2_amplicons.gtf > $x.count.txt; done;

#for f in $run_folder/R1_files/*.fastq; do bowtie $resource_folder/bowtie_index_V1p2/bw1_sparsq_V1p2_amplicons_index $f $f.sam --best -v 3 -k 1 -m 1 -S; done && for x in $run_folder/R1_files/*.sam; do python3 -m HTSeq.scripts.count -f sam -t CDS $x $resource_folder/bowtie_index_V1p2/sparsq_V1p2_amplicons.gtf > $run_folder/R1_files/samcounts/$x.count.txt; done;
#$f $run_folder/R1_files/$f.sam --best -v 3 -k 1 -m 1 -S; done && for x in $run_folder/R1_files/*.sam; \
#note from L - bowtie refuses to work if full path is given to the sam file name
mv $run_folder/R1_files/*.count.txt $run_folder/R1_files/samcounts

echo "	*Running analysis...*"
mkdir -p $run_folder/results
python3 $resource_folder/sparseq_analysisNS_V5.1.py $run_folder $resource_folder

echo "	*Remove trimmed fastq files*"
rm $run_folder/R1_files/*trimmed.fastq

echo "	*Zipping FASTQ files...*"
gzip $run_folder/R1_files/*.fastq

echo "	*Remove sam files*" #will prompt user
rm $run_folder/R1_files/*.sam

echo "*****"$(date)" SPAR-Seq Pipeline is now complete. *****"
