#!/usr/bin/python

import os, sys
import argparse
import re
import openpyxl
from openpyxl import Workbook
import pandas as pd

###This version is for R1s and R2s with UMIs --
#The purpose is to get paired reads and get exact match for row and column UMI from the matching pair of R1 and R2

def main():

    parser = argparse.ArgumentParser(description = 'This outputs an alignment file for each amplicon and \'important_vars_detailed_v3.txt\'. This important_vars file is what quickly tells us if there is a variant')
    parser.add_argument('run_folder', help = 'Full path of run folder with fastq list file')
    parser.add_argument('resource_folder', help = 'Full path of SPAR Seq resource folder')
    #parser.add_argument('--start', help='Stage to start in [select_vars|all_top|bowtie_count|variant_agg] Default: select_vars')
    #parser.add_argument('--intermediates', help = 'Keep Intermediates [Y|N] DEFAULT: N') #Optional argument: --intermediates option
    args = parser.parse_args()  #Needed for args to be recognized

    args.run_folder = os.path.abspath(args.run_folder) #abspath
    args.resource_folder = os.path.abspath(args.resource_folder)

    #parsing script version
    versionre = re.match(r'.*/.[^/]*(V.+)\.py', sys.argv[0])
    version = versionre.group(1)

    #runID
    runIDre = re.match(r'.*(run.[^/]*).*', sys.argv[1])
    runID = runIDre.group(1)

    #script
    scriptre = re.match(r'.*/(.[^/]+)\.py', sys.argv[0])
    script = scriptre.group(1)

    sampleID_listR1 = [] #used to check for duplicate sampleIDs
    input_file_setR1 = []
    sampleID_listR2 = [] #used to check for duplicate sampleIDs
    input_file_setR2 = []

#get R1 AND R2 as original fastq files
    for fastq in os.listdir(args.run_folder+"/R1andR2_files/"):
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
    umirowtable = args.resource_folder+'/UMIrowmatch.csv'
    with open(umirowtable) as rowlines:
        for line in rowlines:
            info = line.strip()
            infol = info.split(",")
            row_dic[(infol[0])] = infol[1]
    rowlines.close()

    col_dic = {}
    coltable = args.resource_folder+'/UMIcolmatch.csv'
    with open(coltable) as collines:
        for line in collines:
            info = line.strip()
            infol = info.split(",")
            col_dic[(infol[0])] = infol[1]
    collines.close()

    #open output dump file
    bcoutput_file = open(args.run_folder+"/results/UMI_analysis_"+version+".txt", "w") #fileout
    bcoutput_file.write("sample, real_well, bc_well, count\n")

    counts_wbname = args.run_folder+"/results/"+runID+"_UMI_PairedAnalysis.xlsx"
    counts_wb = Workbook()
    ws_counts = counts_wb.active
    ws_counts.title = "UMIs_"+runID
    ws_counts.append(["sampleID", "Actual_Well", "BC_Well", "Count"])


    counterone = 0
    while counterone < len(input_file_setR1):
    #for itemnum in len(input_file_setR1):

        #it+=1
        r1file = input_file_setR1[counterone]
        r2file = input_file_setR2[counterone]

        countfilepathR1 = args.run_folder+"/R1andR2_files/"+r1file
        countfilepathR2 = args.run_folder+"/R1andR2_files/"+r2file

        #print("processing file pair " + str(counterone+1) + "\n" + r1file + "\n" + r2file)

        #Using regex to parse sampleID; include well to account for duplicates e.g. 6062021_H2O-multi_S41_R2C3_23M_S41_R1_001.fastq
        #matches = re.match(r'(.*?)_(.*?)_S.{1,4}_.*?_(.{2,3})_.*', filename)
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
            welldict = {}
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
                    #what if val is not in dict
                    if read1 in row_dic and read2 in col_dic:
                        givenrow = row_dic[read1]
                        givencol = col_dic[read2]
                        givenwell = givencol + givenrow
                    else:
                        givenwell = "invalid"

                    if givenwell not in welldict : welldict[givenwell] = 1
                    else : welldict[givenwell] = 1+welldict[givenwell]



                if it%4 == 1 and idlineR1 != idlineR2:
                    print("Read pair error in files " + r1file + "and" + r2file)

        #after the files for this sample are done looping you can unload the dict into a file
        for x in welldict:
            linetowrite = sampleID + "," + well + "," + x + "," + str(welldict[x]) + "\n"
            bcoutput_file.write(linetowrite)
            ws_counts.append([sampleID, well, x, welldict[x]])

        f1.close()
        f2.close()
        counterone+=1


    bcoutput_file.close()
    counts_wb.save(filename = counts_wbname) #save wb
    print("\t\t\t> R1 and R2 UMI analysis now complete...")


main()




#ACTGAATATAAACCAACTGTTCGTTTAGCAACTCGTGATGATGTTCCACGTGCAGTTCGTACTTTAGCAGCAGCATTTGCAGATTATCCAGCAACTCGTCATACTGTTGATCCAGATCGTCATATTGAACGTGTTACTGAATTACAAGAATTATTTTTAACTCGTGTTGGTTTAGATATTGGTAAAGTTTGGGTTGCAGATGATGGTGCAGCAGTTGCAGTTTGGACTACTCCAGAATCTGTTGAAGCAGGTGCAGTTTTTGCAGAAATTGGTCCACGTATGGCAGAATTATCTGGTTCTCGTTTAGCAGCACAACAACAAATGGAAGGTTTATTAGCACCACATCGTCCAAAAGAACCAGCATGGTTTTTAGCAACTGTTGGTGTTTCTCCAGATCATCAAGGTAAAGGTTTAGGTTCTGCAGTTGTTTTACCAGGTGTTGAAGCAGCAGAACGTGCAGGTGTTCCAGCATTTTTAGAAACTTCTGCACCACGTAATTTACCATTTTATGAACGTTTAGGTTTTACTGTTACTGCAGATGTTGAAGTTCCAGAAGGTCCACGTACTTGGTGTATGACTCGTAAACCAGGTGCATAA
