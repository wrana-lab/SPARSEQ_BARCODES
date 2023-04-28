#/bin/bash
#Originally written by Lauren C, Feb 2022
#updated Dec 2022 for further use

#quick wrapper script for quasispecies Analysis workflow

start=`date +%s`

#get SE_fastq_files
echo "Unzipping fastq files..."
gunzip *.gz

ls *.fastq > SE_fastq_files.txt

##run first python scripts
echo "Running sequence pileup script..."
python3 sparseq_srbdv2.py

#take those outputs and put them in clustal
# --force means it will force overwrite
echo "Aligning Srbd..."
clustalo -i filtered_sv2_sequences.fa  -o filtered_sv2_sequences.fasta --force

#run second python script with the aligned fasta files
echo "Running fasta cutter script..."
python3 fasta_cutter_barcode_v3.py

#Zip fastq files back up
echo "Zipping fastq files..."
gzip *.fastq

end=`date +%s`

#report run time
runtime=$((end-start))
echo "Runtime was $runtime seconds"
