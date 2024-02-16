#!/usr/bin/env python3
# # -*- coding: utf-8 -*-

"""
Created on Mon Feb 20 16:14:07 2020

@author: gottince
"""
#----------------------------------#
#    DESCRIPTION
#----------------------------------#

# Read table gene model file from annotation transfert output and test for canonical model
# Retrieve gene ID for wich non-canonical model is detected

#----------------------------------#
#    PARAMETRES
#----------------------------------#
# - File name of input genomic fasta
# - File name of input Table gene file 
# - File name of output Table gene file (with corrected intron)
import sys
import argparse
import csv
from Bio import SeqIO
from Bio.Seq import translate, reverse_complement


parser = argparse.ArgumentParser()

parser.add_argument("-f","--fasta", type=str, help="Exact path of genomic fasta file")
parser.add_argument("-t","--table", type=str, help="Exact path of feature table file")
parser.add_argument("-o","--output", type=str, help="Exact path of output file")

args = parser.parse_args()


#----------------------------------#
#            FUNCTIONS
#----------------------------------#

def isCanonicalIntron_forward(DNA_dict,Chr,startIntron,stopIntron) :
    "Check if donnor and acceptor splice sites are canonical in forward strand"
    donnor=DNA_dict[Chr][startIntron:startIntron+2].seq.upper()
    acceptor=DNA_dict[Chr][stopIntron-3:stopIntron-1].seq.upper() 
    return (donnor=="GT" and acceptor=="AG") or (donnor=="GC" and acceptor=="AG")
 

def isCanonicalIntron_reverse(DNA_dict,Chr,startIntron,stopIntron) :
    "Check if donnor and acceptor splice sites are canonical in reverse strand"
    donnor=DNA_dict[Chr][stopIntron-3:stopIntron-1].seq.upper() 
    acceptor=DNA_dict[Chr][startIntron:startIntron+2].seq.upper()
    return (donnor == "AC" and acceptor == "CT") or (donnor == "GC" and acceptor == "CT") 


def noStart(codon) :
    if (codon.upper()=="ATG") :
        return False
    return True


def noStop(codon) :
    if (codon.upper() in ("TAA", "TAG", "TGA")):
        return False
    return True




#----------------------------------#
#              MAIN
#----------------------------------#
outfile=open(args.output,"w")
outfile.write("Chr\tProtID\tNOstart\tNOstop\tFS\tunexpectedSplicingSite\tStopInFrame\tLengthPb\tValidity\n")

## Read genome
chr_dict = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

## read table with one gene per line / one line per gene
gff=open(args.table, mode='r')
gff_reader = csv.reader(gff, delimiter=';')

for row in gff_reader :
    start=""
    stop=""
    
    frameshift=False
    notCanonic=False
    NCintron=False
    myProt=row[0]
    strand=row[1]
    #Chr='_'.join(row[0].split("_")[1:-1]) ## Chromosome Id
    Chr=row[0].split("_")[0] ## pour Nip
    ## VR add stopInFrame and CDS length
    exons = []
    for i in range(2, len(row) - 1, 2):
        # index in Python starts at 0 hence -1
        start, end = int(row[i]) - 1, int(row[i + 1]) - 1
        exons.append(str(chr_dict[Chr][start:(end + 1)].seq))

    cds = ''.join(exons)
    if strand == "-":
        cds = str(reverse_complement(cds))
    extra_nucl=len(cds)%3
    lengthPb = (extra_nucl!=0)
    if lengthPb != 0:
        cds = cds[:-extra_nucl]

    prot = translate(cds, to_stop=False)
    stopInFrame = '*' in prot[:-1]
    ## end VR

    if(strand=="+") :
        ## Check start and stop
        start=chr_dict[Chr][int(row[2])-1:int(row[2])+2].seq
        stop=chr_dict[Chr][int(row[-1])-3:int(row[-1])].seq
        ## Check intron/frameshift
        ## if frameshift or (intron and not canonical intron)
        for i in range(3,len(row)-2,2) :
            if(int(row[i+1])<=int(row[i])+25) :
                frameshift=True
            elif(int(row[i+1])>int(row[i])+25 and not isCanonicalIntron_forward(chr_dict,Chr,int(row[i]),int(row[i+1]))) :
                NCintron=True

    else :
        ## Check start and stop
        stop=chr_dict[Chr][int(row[2])-1:int(row[2])+2].reverse_complement().seq
        start=chr_dict[Chr][int(row[-1])-3:int(row[-1])].reverse_complement().seq
        ## Check intron/frameshift
        ## if frameshift or (intron and not canonical intron)
        for i in range(3,len(row)-2,2) :
            if(int(row[i+1])<=int(row[i])+25):
                frameshift=True
            elif(int(row[i+1])>int(row[i])+25 and not isCanonicalIntron_reverse(chr_dict,Chr,int(row[i]),int(row[i+1]))) :
                NCintron=True
   
    non_canonique =noStart(start) or noStop(stop) or frameshift or NCintron or stopInFrame or lengthPb;
    annotation_status = "notValid" if non_canonique else "Valid"
    line_info=[Chr,myProt,str(noStart(start)),str(noStop(stop)), str(frameshift),str(NCintron),str(stopInFrame),str(lengthPb),annotation_status];
    outfile.write("\t".join(line_info)+ "\n")
    
gff.close()    
outfile.close()
