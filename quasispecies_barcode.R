#feb 2022
#lauren c
#a simple script to set up some visualization and analysis of the quasispecies analysis outputs

library(dplyr)
#BiocManager::install("Biostrings")
library(Biostrings) #used for doing translations within R

#setwd("...")

### translate top sequences from python workflow to see if there are variants or AA changes
#need to select sequences, remove - gaps, get correct reading frame, do any other trimming required

#Reference sequences are used to write in output files for ease of viewing alignments
#usually used WT seq for WT samples and alpha samples, then delta for delta samples, omicron for omicron samples, etc
sv2_wt<-"ATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTT"
sv2_wt_AA<-"IYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGV"

# sv2_delta<-"ATCTATCAGGCCGGTAGCAAACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTT"
# sv2_delta_AA<-"IYQAGSKPCNGVEGFNCYFPLQSYGFQPTNGV"

# sv2_omicron<-"ATCTATCAGGCCGGTAACAAACCTTGTAATGGTGTTGCAGGTTTTAATTGTTACTTTCCTTTACGATCATATAGTTTCCGACCCACTTATGGTGTT"
# sv2_omicron_AA<-"IYQAGNKPCNGVAGFNCYFPLRSYSFRPTYGV"


#import aligned list from fasta_cutter_barcode_v3.py
srbdv2_aligned_list<-read.csv("srbdv2_aligned_list.txt")
head(srbdv2_aligned_list)
#drop refseq
srbdv2_aligned_list<-srbdv2_aligned_list[2:nrow(srbdv2_aligned_list),]

#count occurences of sequences
srbdv2_counted_seqs<-as.data.frame(srbdv2_aligned_list %>% dplyr::count(sequence))
srbdv2_counted_seqs<-srbdv2_counted_seqs[order(srbdv2_counted_seqs$n, decreasing = T),]

#adjust edges of sequences
srbdv2_counted_seqs$sequence<-gsub("-", "", srbdv2_counted_seqs$sequence)

# depending on version of pipeline may need to skip first starting -AA and final C to keep sequences trimmed and lined up
srbdv2_counted_seqs$sequence<-gsub("^AAA", "A", srbdv2_counted_seqs$sequence)
srbdv2_counted_seqs$sequence<-gsub("^AA", "A", srbdv2_counted_seqs$sequence)
srbdv2_counted_seqs$sequence<-gsub("C$", "", srbdv2_counted_seqs$sequence)


#take top sequences - here we settled on keeping all sequences
srbdv2_top<-subset(srbdv2_counted_seqs, srbdv2_counted_seqs$n > 0)

srbdv2_dnastring<-DNAStringSet(srbdv2_top$sequence)
#convert to AAs
srbdv2_aas<-as.data.frame(translate(srbdv2_dnastring))
#bind
srbdv2_translated<-cbind(srbdv2_aas,srbdv2_top)
colnames(srbdv2_translated)<-c("aas", "sequence", "numSamples")


#organize and write out finalized fa files
#write.table(">Refseq_srbdv2_omicron", file = "delta_srbdv2_topsequences.fa", row.names = FALSE, append = FALSE, col.names = FALSE, quote = F)
#write.table(sv2_omicron, file = "delta_srbdv2_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
#write.table(">Refseq_srbdv2_delta", file = "delta_srbdv2_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
#write.table(sv2_delta, file = "delta_srbdv2_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(">Refseq_srbdv2_WT", file = "alpha_srbdv2_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(sv2_wt, file = "alpha_srbdv2_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
for(i in 1:nrow(srbdv2_translated)){
  seq<-srbdv2_translated[i,2]
  numsamples<-srbdv2_translated[i,3]
  lines<-paste0(">",numsamples,"_samples")
  write.table(lines, file = "alpha_srbdv2_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
  write.table(seq, file = "alpha_srbdv2_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
}

#AAs
#write.table(">Refseq_srbdv2_omicron", file = "delta_srbdv2_topAAsequences.fa", row.names = FALSE, append = FALSE, col.names = FALSE, quote = F)
#write.table(sv2_omicron_AA, file = "delta_srbdv2_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
#write.table(">Refseq_srbdv2_delta", file = "delta_srbdv2_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
#write.table(sv2_delta_AA, file = "delta_srbdv2_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(">Refseq_srbdv2_WT", file = "alpha_srbdv2_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(sv2_wt_AA, file = "alpha_srbdv2_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)

for(i in 1:nrow(srbdv2_translated)){
  seq<-srbdv2_translated[i,1]
  numsamples<-srbdv2_translated[i,3]
  lines<-paste0(">",numsamples,"_samples")
  write.table(lines, file = "alpha_srbdv2_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
  write.table(seq, file = "alpha_srbdv2_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
}

#next would run clustal omega locally again to get the aligned fasta files for visualizing in snapgene
# eg clustalo -i omicron_spbs_topsequences.fa -o omicron_spbs_topsequences.fasta
