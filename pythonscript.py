#!/usr/local/bin/python3

import os
import sys
import subprocess
import shutil

# Split the DNA into coding and non-coding parts, write each of the sub-sequences into individual files AND two separate files.

# First import local DNA seq.

subprocess.call("cp /localdisk/data/BPSM/Lecture11/plain_genomic_seq.txt .", shell = True)
print("\n".join(os.listdir())) #see the files that now exist
localDNAseq = open("plain_genomic_seq.txt").read().rstrip() #need to take away the new line at the end

---------

#Now import remote DNA: the bellow is done in linux shell or use os.system("")

wget -qO AJ223353.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AJ223353&strand=1&rettype=fasta&retmode=text"

grep -v ">" > AJ223353_noheader.fasta

#Alternatively, we can use esearch (part of edirect software package)

esearch -db nucleotide -query "AJ223353[accession]" | efetch -db nucleotide -format fasta | grep -v ">" > AJ223353_noheader.fasta

# We can also access the file from GenBank

wget -qO AJ223353.genbank "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AJ223353&strand=1&rettype=gb&retmode=text"

#alternatively, can get the sequence from GenBank using edirect

esearch -db nucleotide -query "AJ223353[accession]" | efetch -db nucleotide -format gb > AJ223353.genbank

#Using genbank, we can extract the segment that is coding

cat AJ223353.genbank | grep "CDS" 

------- 

#back to python

remotefile = open("AJ223353_noheader.fasta").read().upper()
localfile = open("plain_genomic_seq.txt").read().upper()

#convert the remote file into a single line

remotefile_singleline = remotefile.replace("\n","")
remotefile_singleline_ready = remotefile_singleline

local_seq_singleline = local_seq.replace("\n","")

# Remove non-DNA characters: First check what needs to be replaced

remotefile_singleline_anythingleft = remotefile_singleline.replace("G","").replace("A","").replace("T","").replace("C","")
remotefile_singleline_anythingleft

local_seq_singleline_anythingleft = local_seq_singleline.replace("G","").replace("A","").replace("T","").replace("C","")
local_seq_singleline_anythingleft

# Now replace the characters in local DNA

local_seq_singleline_trueDNA = local_seq_singleline.replace("X", "").replace("S", "").replace("K", "").replace("L", "")
local_seq_singleline_trueDNA

# Alternatively, use: set(list(local_seq.rstrip()))

# Get the remote coding/non-coding sequences

remote_noncoding01 = remotefile_singleline_ready[0:28]
remote_exon01 = remotefile_singleline_ready[28:409]
remote_noncoding02 = remotefile_singleline_ready[409:]

# Get the local coding/non-coding sequences
local_exon01 = local_seq_singleline_ready[0:63]
local_intron01 = local_seq_singleline_ready[63:90]
local_exon02 = local_seq_singleline_ready[90:]

# Now we write the remote sequences to new file

remote_noncoding01_out = open("remote_noncoding01.fasta", "w")
remote_noncoding01_out.write(">AJ223353_noncoding01_length" + str(len(remote_noncoding01)) + "\n") #write the fasta header
remote_noncoding01_out.write(remote_noncoding01)
remote_noncoding01_out.close()
print(open("remote_noncoding01.fasta").read())

# Write the local sequences to new file

local_exon01_out = open("local_exon01.fasta", "w")
local_exon01_out.write(">LocalSeq_exon01_length" + str(len(local_exon01)) + "\n")
local_exon01_out.write(local_exon01)
local_exon01_out.close()
print(open("local_exon01.fasta").read())

#Do this for intron 1

local_intron01_out = open("local_intron01.fasta", "w")
local_intron01_out.write(">LocalSeq_intron01_length" + str(len(local_intron01)) + "\n")
local_intron01_out.write(local_intron01)
local_intron01_out.close()
print(open("local_intron01.fasta").read())

#Now for exon 2

local_exon02_out = open("local_exon02.fasta", "w")
local_exon02_out.write(">LocalSeq_exon02_length" + str(len(local_exon02)) + "\n")
local_exon02_out.write(local_exon02)
local_exon02_out.close()
print(open("local_exon02.fasta").read())

#Now make a file with all the exons
exons_out = open("All_exons.fasta", "w")
exons_out.write(">AJ223353_exon01_length" + str(len(remote_exon01)) + "\n" + remote_exon01 + "\n")
exons_out.write(">LocalSeq_exon01_length" + str(len(local_exon01)) + "\n" + local_exon01 + "\n")
exons_out.write(">LocalSeq_exon02_length" + str(len(local_exon02)) + "\n" + local_exon02)
exons_out.close()
print(open("All_exons.fasta").read())

#Now make a file with all the introns

introns_out = open("All_noncodings.fasta", "w")
introns_out.write(">AJ223353_noncoding01_length" + str(len(remote_noncoding01)) + "\n" + remote_noncoding01 + "\n")
introns_out.write(">AJ223353_noncoding02_length" + str(len(remote_noncoding02)) + "\n" + remote_noncoding02 + "\n")
introns_out.write(">LocalSeq_intron01_length" + str(len(local_intron01)) + "\n" + local_intron01)
introns_out.close()
print(open("All_noncodings.fasta").read())

#Now run blast

blastx -db xxx -query fastafile -outfmt 7 > someoutfile

blastx -db xxx -query fastafile -outfmt 7 | head -n6 #Parse to only get the top result (the first 5 lines of a BLAST output are comments) 

blastx -db $db -query remote_exon01.fasta -outfmt 7 | head -n6 | tail -n1 > remote_exon01.tophit

