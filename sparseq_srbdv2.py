#!/usr/bin/python

#initial script for quasispecies analysis of paired end sparseq barcode data

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt

import os, sys
import argparse
import re
import openpyxl
from openpyxl import Workbook
import pandas as pd
from kneed import KneeLocator

###This version is for R1s and R2s with barcodes --
#The purpose is to get paired reads and get exact match for row and column barcode from the matching pair of R1 and R2
#and run the initial QS analysis on them - sorting reads.
Srbdv2 = "ATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTT" #WT

def main():
    Srbdv2_var_dic = dict()
    Srbdv2_sum = open("Srbdv2_summary.csv", "w")
    Srbdv2_sum.write("sample" + "," + "sequence" + "," + "count" + "\n")
    fileout_Srbdv2 = open("Srbdv2_sequence_list.fa", "w")
    fileout_Srbdv2.write(">Refseq_Srbdv2_WT" + "\n" + Srbdv2 + "\n")
    knee_tracker = open("kneepoint_tracker.csv", "w")
    knee_tracker.write("sampleID" + "," + "kneepoint" + "\n")
    filtered_sequences = open("filtered_sv2_sequences.fa", "w")
    filtered_sequences.write(">Refseq_Srbdv2_WT" + "\n" + Srbdv2 + "\n")
    seqpercenttablefile = open("seq_percents.csv", "w")

    sampleID_listR1 = [] #used to check for duplicate sampleIDs
    input_file_setR1 = []
    sampleID_listR2 = [] #used to check for duplicate sampleIDs
    input_file_setR2 = []

#get R1 AND R2 as original fastq files
    with PdfPages('sv2_knees.pdf') as pdf:
        for fastq in os.listdir("./"):
            #edited this to use the trimmed fastq
            if "R1_001.fastq" in fastq:
                matches = re.match(r'.*?_(.*?)_.*?', fastq)
                sampleID_listR1.append(matches.group(1))
                input_file_setR1.append(fastq)
            elif "R2_001.fastq" in fastq:
                matches = re.match(r'.*?_(.*?)_.*?', fastq)
                sampleID_listR2.append(matches.group(1))
                input_file_setR2.append(fastq)

        input_file_setR1.sort()
        input_file_setR2.sort()

        #create list of duplicate samples
        sIDR1 = set()
        sampleIDR1_dups = {x for x in sampleID_listR1 if x in sIDR1 or (sIDR1.add(x) or False)}
        #since we process as pairs, only need r1

        #get the row and column dicts
        row_dic = {}
        bcrowtable = './BCrowmatch.csv'
        with open(bcrowtable) as rowlines:
            for line in rowlines:
                info = line.strip()
                infol = info.split(",")
                row_dic[(infol[0])] = infol[1]
        rowlines.close()

        col_dic = {}
        coltable = './BCcolmatch.csv'
        with open(coltable) as collines:
            for line in collines:
                info = line.strip()
                infol = info.split(",")
                col_dic[(infol[0])] = infol[1]
        collines.close()

        counterone = 0
        while counterone < len(input_file_setR1):
        #for itemnum in len(input_file_setR1):
            print("Processing file ", counterone+1)

            #it+=1
            r1file = input_file_setR1[counterone]
            r2file = input_file_setR2[counterone]

            countfilepathR1 = "./"+r1file
            countfilepathR2 = "./"+r2file

            #pasre sample ID with regex
            matches = re.match(r'(.*?)_(.*?)_S.{1,4}_.*?_(.{2,3})_.*?_.*', r1file)
            date = matches.group(1)
            sampleID = matches.group(2)
            sampleID = re.sub("-V1-2", "", sampleID)
            well = matches.group(3)

            if sampleID in sampleIDR1_dups:
                sampleID = sampleID+"_"+well
            count_dataR1 = [sampleID]

            #open pair of fastqs
            it=3 #start at it = 3
            with open(countfilepathR1) as f1, open(countfilepathR2) as f2:
                #start at 0, get ID lines, add 1, get seq lines, then add 2 and get next ID lines again
                long_read_dict_Srbdv2 = dict()
                for r1line,r2line in zip(f1,f2):
                    it+=1 #add one to it immediately

                    if it%4 == 0:
                        #get ID line (first line, 5th line...) check if ID lines match, and remove " 1" and " 2"
                        #you do this to make sure it's a proper pair of reads
                        idlineR1 = re.sub(" 1", "", r1line)
                        idlineR2 = re.sub(" 2", "", r2line)

                    #if ID lines match AND it's a 2nd, 6th etc line  get sequences on next line
                    if it%4 == 1 and idlineR1 == idlineR2:
                        read1 = r1line[:5] #the 5 is not included, it is 0:4
                        read2 = r2line[:5]
                        #match reads to barcode list
                        #read 1 is row (letter) and r2 is cols (num)
                        if read1 in row_dic and read2 in col_dic:
                            givenrow = row_dic[read1]
                            givencol = col_dic[read2]
                            givenwell = givencol + givenrow

                            if givenwell == well:
                            #if correct well pair get seq and
                                #this should be seq
                                #seq = r1line[:1].strip()
                                seq = r1line.strip()

                                if len(seq)>100 :

                                    if "ACCTTTTGAGA" in seq : ##Srbdv2
                                    #usually 27 - -22 but here we have additional 5 bp for barcode
                                    #the only part of rev primer involved here is ACCAACCATAC so 11
                                        seq = seq[32:] #
                                        #remove leading AA to clean up alignments
                                        seq = re.sub("^AAA", "A", seq)
                                        seq = re.sub("^AA", "A", seq)
                                        #hard trim so it is 101 long
                                        seq = seq[:96]
                                        if seq not in long_read_dict_Srbdv2 : long_read_dict_Srbdv2[seq] = 1
                                        else : long_read_dict_Srbdv2[seq] = 1+long_read_dict_Srbdv2[seq]

            Srbdv2_ = list()
            Srbdv2_c_ = 0
            #Srbdv2 >> trimmed version above is 101
            long_read_t_Srbdv2 = long_read_dict_Srbdv2.items()
            for key,val in long_read_t_Srbdv2 :
                #append to dict
                #have to have low level cutoffs of counts and lengths within reason
                if val>100 and len(key)<105 and len(key)>85:
                    Srbdv2_.append((val, key))
                    Srbdv2_c_ = Srbdv2_c_+ val
                    #Srbdv2_c_ is total counts of sample sv2, including dropped top seq

            Srbdv2_.sort(reverse=True)
            if Srbdv2_: #drop top count sequence
                Srbdv2_.pop(0)

            #run knee plot testing
            if Srbdv2_:
                #if there is anything in Srbdv2_ after dropping top
                #convert dict to pandas frame
                seqs_merged = pd.DataFrame.from_records(Srbdv2_)
                #fix colnames
                seqs_merged.columns = ['count', 'sequence']
                #get percent within sample
                seqs_merged['percent_count'] = seqs_merged['count'] / Srbdv2_c_ * 100
                rangelength = len(seqs_merged)
                #add dummy variable
                seqs_merged['seq_number'] = range(1,rangelength+1)
                seqs_merged['sample'] = sampleID

                if len(seqs_merged) > 2: #if long enough to run kneedle
                    #https://github.com/arvkevi/kneed; https://kneed.readthedocs.io/en/stable/parameters.html#s
                    kneedle1 = KneeLocator(seqs_merged['seq_number'], seqs_merged["percent_count"], S=2.0, curve="convex", direction="decreasing")
                    #get y axis value
                    if kneedle1.knee_y:
                        #plot actual data
                        plt.figure()
                        plt.suptitle(sampleID)
                        plt.plot(seqs_merged['seq_number'], seqs_merged["percent_count"], marker='o', linestyle = 'dotted')
                        plt.ylabel('percent of total counts')
                        plt.xlabel('sequence #')
                        pdf.savefig()
                        #now plot kneeplot
                        kneepoint = round(kneedle1.knee_y, 2)
                        kneedle1.plot_knee()
                        pdf.savefig()
                        if kneepoint > 0.5:
                            filtered_seqs_merged = seqs_merged[seqs_merged["percent_count"] > kneepoint]
                            knee_tracker.write(sampleID + "," + str(kneepoint) + "\n")
                            filtered_seqs_merged.to_csv('seq_percents.csv',mode='a', header=False)
                        else:
                            filtered_seqs_merged = seqs_merged[seqs_merged["percent_count"] > 0.5]
                            knee_tracker.write(sampleID + "," + "kneepoint too low: used 0.5% cutoff" + "\n")
                            filtered_seqs_merged.to_csv('seq_percents.csv',mode='a', header=False)

                        filt_dict = pd.Series(filtered_seqs_merged['count'].values,index=filtered_seqs_merged['sequence']).to_dict()
                        filt_dict_items = filt_dict.items()
                        for key,val in filt_dict_items :
                            filtered_sequences.write(">" + sampleID + "_" + str(val) + "." + key + "." + "\n" + key + "\n")

                    else:
                        #if knee not found just apply cutoff
                        filtered_seqs_merged = seqs_merged[seqs_merged["percent_count"] > 0.5]
                        filt_dict = pd.Series(filtered_seqs_merged['count'].values,index=filtered_seqs_merged['sequence']).to_dict()
                        filt_dict_items = filt_dict.items()
                        for key,val in filt_dict_items :
                            filtered_sequences.write(">" + sampleID + "_" + str(val) + "." + key + "." + "\n" + key + "\n")
                        knee_tracker.write(sampleID + "," + "could not fit model: used 0.5% cutoff" + "\n")
                        filtered_seqs_merged.to_csv('seq_percents.csv',mode='a', header=False)

                elif len(seqs_merged) > 1: #if less than 3 but more than 1 add sequences
                    knee_tracker.write(sampleID + "," + "only 2 sequences so excluded from further analysis" + "\n")

                else:
                    knee_tracker.write(sampleID + "," + "only 1 sequence so excluded from further analysis" + "\n")
                    Srbdv2_ = list()

            #this second out is the complete set of sequences.
            for val,item in Srbdv2_ :
                fileout_Srbdv2.write(">" + sampleID + "-" + str(val) + "." + item + "." + "\n" + item + "\n")
                Srbdv2_sum.write(sampleID + "," + item + "," + str(val) + "\n")
            counterone+=1


    filtered_sequences.close()
    knee_tracker.close()
    Srbdv2_sum.close()
    fileout_Srbdv2.close()
    seqpercenttablefile.close()
    print("\t\t> R1 and R2 barcode QS analysis now complete...")

main()
