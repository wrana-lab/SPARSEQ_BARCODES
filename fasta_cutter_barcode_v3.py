#!/usr/bin/python

##script for cutting up fasta files from the quasispecies script
#after they have been aligned in clustal

import re

#set up output file, write headers
fileout_Srbdv2 = open("srbdv2_aligned_list.txt", "w")
fileout_Srbdv2.write("sample, sequence\n")
srbdv2_sequences = list()

#due to the way I set up the previous python script to have seqs of specific lengths
#the line structure for these  files should be consistent
#if the structure breaks, and the sequences are too long or too short,
#it means that the sequences are not being filtered properly and low quality sequences
#are being allowed through that are causing gaps in the alignments

#Srbdv2 processing
firstcount = 0
secondcount = 1
thirdcount = 2
#
print("Processing Srbdv2...")
with open("filtered_sv2_sequences.fasta") as fasta:
    alllines = fasta.readlines()
    #check num of lines in file
    number_of_lines = len(alllines)
    #limitline = number_of_lines/3
    while firstcount < number_of_lines :
        firstline = alllines[firstcount].strip()
        secondline = alllines[secondcount].strip()
        thirdline = alllines[thirdcount].strip()

        #remove sequence and keep sample ID and count
        firstline = re.sub(">", "", firstline)
        firstline_trimmed = re.sub("\..*$", "", firstline)
        sequence = secondline + thirdline #+ fourthline
        srbdv2_sequences.append(sequence)
        fileout_Srbdv2.write(firstline_trimmed + "," + sequence + "\n")
        firstcount = firstcount+3
        secondcount = secondcount+3
        thirdcount = thirdcount+3

fasta.close()

print("Done.")
fasta.close()
fileout_Srbdv2.close()

###Part 2: Crunching
#crunch the numbers of the contents of the sequences per base
#this goes into r to make bar graphs
#basically col 1 is a position in the sequence and col 2 is A/T/C/G/-
print("Now crunching bases...")

#set up file, write headers
fileout2_srbdv2 = open("srbdv2_crunched_bases.txt", "w")
fileout2_srbdv2.write("position, base\n")

#get length of sequences (ie the nchar)
srbdv2_seq_length = len(srbdv2_sequences[0])

##get length of each list
num_of_seqs_srbdv2 = len(srbdv2_sequences)

counter = 1
for sequ in srbdv2_sequences :
    if counter%100==0 :
        print("Srbdv2: sequence " + str(counter) + " out of " + str(num_of_seqs_srbdv2) )
    for nchar in range(srbdv2_seq_length) :
        base = sequ[nchar]
        pos = nchar+1
        fileout2_srbdv2.write( str(pos) + "," + base + "\n")
    counter = counter+1

#close file
fileout2_srbdv2.close()
