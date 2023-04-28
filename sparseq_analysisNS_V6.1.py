#!/usr/bin/python

import os, sys
import argparse
import re
import openpyxl
from openpyxl import Workbook
import pandas as pd

###This version is for R1s and R2s with barcodes --
#all this does is reformat bowtie outputs for R1 and R2 barcodes

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

    for fastq in os.listdir(args.run_folder+"/R1andR2_files/samcounts"):
        #edited this to use the trimmed fastq
        if "fastq.trimmed.R1.fastq.sam.R1count" in fastq:
            matches = re.match(r'.*?_(.*?)_.*?', fastq)
            sampleID_listR1.append(matches.group(1))
            input_file_setR1.append(fastq)
        elif "fastq.trimmed.R2.fastq.sam.R2count" in fastq:
            matches = re.match(r'.*?_(.*?)_.*?', fastq)
            sampleID_listR2.append(matches.group(1))
            input_file_setR2.append(fastq)

    input_file_setR1.sort()
    input_file_setR2.sort()

    #create list of duplicate samples
    sIDR1 = set()
    sampleIDR1_dups = {x for x in sampleID_listR1 if x in sIDR1 or (sIDR1.add(x) or False)}
    sIDR2 = set()
    sampleIDR2_dups = {x for x in sampleID_listR2 if x in sIDR2 or (sIDR2.add(x) or False)}


    #defining outputs
    bctable_fileR1 = open(args.run_folder+"/results/R1BC_count_table_"+version+".txt", "w") #fileout
    bctable_fileR1.write("samples")
    bctable_fileR2 = open(args.run_folder+"/results/R2BC_count_table_"+version+".txt", "w") #fileout
    bctable_fileR2.write("samples")


    #generating workbook in parallel
    #.active sets active sheet in excel file
    counts_wbnameR1 = args.run_folder+"/results/"+runID+"_R1BC_Table.xlsx"
    counts_wbR1 = Workbook()
    ws_countsR1 = counts_wbR1.active
    ws_countsR1.title = "R1BC_count_"+version

    #R2s
    counts_wbnameR2 = args.run_folder+"/results/"+runID+"_R2BC_Table.xlsx"
    counts_wbR2 = Workbook()
    ws_countsR2 = counts_wbR2.active
    ws_countsR2.title = "R2BC_count_"+version

    it=0
    remove_data = ["__no_feature", "__ambiguous", "__too_low_aQual","__alignment_not_unique"]
    count_header = ["sample","well"]

    for countfile in input_file_setR1:

        it+=1
        countfilepath = args.run_folder+"/R1andR2_files/samcounts/"+countfile

        #Using regex to parse sampleID; include well to account for duplicates e.g. 6062021_H2O-multi_S41_R2C3_23M_S41_R1_001.fastq
        matches = re.match(r'(.*?)_(.*?)_S.{1,4}_.*?_(.{2,3})_.*?_.*', countfile)
        #matches = re.match(r'(.*?)_(.*?)_S.{1,4}_.*?_(.{2,3})_.*', countfile)
        date = matches.group(1)
        sampleID = matches.group(2)
        well = matches.group(3)

        if sampleID in sampleIDR1_dups:
            sampleID = sampleID+"_"+well

        count_dataR1 = [sampleID]
        with open(countfilepath) as f1:

            header= ""
            data = sampleID
            data+=","+well

            for line in f1:
                info = line.strip()
                infol = info.split()

                if infol[0] not in remove_data:

                    if it==1: #if first iteration add header line
                        header+=","+infol[0]
                        count_header.append(infol[0])

                    data+=","+infol[1]
                    count_dataR1.append(int(infol[1]))


            if it == 1:
                bctable_fileR1.write(",well"+header+"\n")
                header+=",filename\n"
                count_header.extend(["filename"])
                ws_countsR1.append(count_header)

        bctable_fileR1.write(data+","+countfile+"\n")
        data+=","+countfile+"\n"
        count_dataR1.extend([countfile])
        #stick the well after the dang sample name
        count_dataR1.insert(1,well)
        ws_countsR1.append(count_dataR1)

    #repeat for R2
    it=0
    remove_data = ["__no_feature", "__ambiguous", "__too_low_aQual","__alignment_not_unique"]
    count_header = ["sample","well"] #excel

    for countfile in input_file_setR2:

        it+=1
        countfilepath = args.run_folder+"/R1andR2_files/samcounts/"+countfile

        #Using regex to parse sampleID; include well to account for duplicates e.g. 6062021_H2O-multi_S41_R2C3_23M_S41_R1_001.fastq
        matches = re.match(r'(.*?)_(.*?)_S.{1,4}_.*?_.*?_(.{2,3})_.*', countfile)
        #matches = re.match(r'(.*?)_(.*?)_S.{1,4}_.*?_(.{2,3})_.*', countfile)
        date = matches.group(1)
        sampleID = matches.group(2)
        well = matches.group(3)

        if sampleID in sampleIDR2_dups:
            sampleID = sampleID+"_"+well

        count_dataR2 = [sampleID]
        with open(countfilepath) as f1:

            header= ""
            data = sampleID
            data+=","+well

            for line in f1:
                info = line.strip()
                infol = info.split()

                if infol[0] not in remove_data:

                    if it==1: #if first iteration add header line
                        header+=","+infol[0]
                        count_header.append(infol[0])

                    data+=","+infol[1]
                    count_dataR2.append(int(infol[1]))

            if it == 1:
                bctable_fileR2.write(",well"+header+"\n")
                header+=",filename\n"
                count_header.extend(["filename"])
                ws_countsR2.append(count_header)

        bctable_fileR2.write(data+","+countfile+"\n")
        data+=","+countfile+"\n"
        count_dataR2.extend([countfile])
        #stick the well after the dang sample name
        count_dataR2.insert(1,well)
        ws_countsR2.append(count_dataR2)

    f1.close()

    bctable_fileR1.close()
    bctable_fileR2.close()
    counts_wbR1.save(filename = counts_wbnameR1) #save wb
    counts_wbR2.save(filename = counts_wbnameR2) #save wb

    #open bc2 counts csv files

    bc2_table = pd.read_excel(counts_wbnameR2, header=0)

    #reorder bc2_table
    bc2_table = bc2_table[['sample', 'well',	'1', '2', '3',	'4', '5',	'6',	'7',	'8',	'9', '10', '11', '12',	'13',	'14',	'15',	'16',	'17',	'18',	'19',	'20',	'21',	'22',	'23',	'24', '__not_aligned', 'filename']]

    counts_wbnameR2_reordered = args.run_folder+"/results/"+runID+"_R2BC_Table_Reordered.xlsx"

    bc2_table.to_excel(counts_wbnameR2_reordered,  sheet_name='Reordered', index=False)


    print("\t\t\t> R1 and R2 BC bowtie_count now complete...")


main()
