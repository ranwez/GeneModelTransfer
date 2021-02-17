# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 15:15:02 2020

@author: gottince

@description: script permettant d'extraire des sequences proteiques d'un
fasta genomique a partir d'un gff.
"""
import argparse
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

#----------------------------------#
#    PARAMETRES
#----------------------------------#
# - File name of genomic fasta
# - output filename
# - output sequence type dna/rna/protein

parser = argparse.ArgumentParser()

parser.add_argument("-f","--fasta", type=str, help="Exact path of genomic fasta file")
parser.add_argument("-g","--gff", type=str, help="Exact path of gff file")
parser.add_argument("-o","--output", type=str, help="Output file name")
parser.add_argument("-t","--type", type=str, choices=["gene","cdna","prot","cds","exon","frameshift"], help="type of extracting sequence : <gene> extract genomique dna sequence with intron, <cdna> extract reconstruct coding nucleotide sequence with frameshift, <prot> extract reconstruct protein sequence, <exons> extract isolated coding exons, <cds> extract reconstruct coding nucleotide sequence without frameshift.")

args = parser.parse_args()

#----------------------------------#
#            FUNCTIONS
#----------------------------------#

# 1. extracting complete gene
def extract_gene(fasta,gff) :
    gff_reader = csv.reader(gff, delimiter='\t')
    
    sid=""
    allseq=[]
    for row in gff_reader :
        if(row[2]=="gene") :
            tmp=row[8].split(';')
            sid=tmp[0][3:]
            subseq=chr_dict[row[0]][int(row[3])-1:int(row[4])]
            if(row[6]=="-") :
                dna=SeqRecord(subseq.reverse_complement().seq, id=sid, description="")
                allseq.append(dna)
            else :
                dna=SeqRecord(subseq.seq, id=sid, description="")
                allseq.append(dna)

    return allseq
            
# 2. extracting coding sequence : prot and cds
def extract_coding(fasta,gff,typeseq) :
    gff_reader = csv.reader(gff, delimiter='\t')
    
    dna= ""
    sid=""
    allseq=[]
    
    for row in gff_reader :
        if(row[2]=="gene") :
            if(len(dna)>0) :
            # Export sequence
                dnaSeq = Seq(dna, IUPAC.ambiguous_dna)
                if(typeseq=="prot") :
                    dnaRec = SeqRecord(dnaSeq.translate(), id=sid, description="")
                    allseq.append(dnaRec)               
                else :
                    dnaRec=SeqRecord(dnaSeq, id=sid, description="")
                    allseq.append(dnaRec)
        
            dna = ""
            sid=row[8][3:]
        elif(row[2]=="CDS" or row[2]=="cds") :
        # Extract sequence
            subseq=chr_dict[row[0]][int(row[3])-1:int(row[4])]
            if(row[6]=="-") :
                #dna=dna+str(subseq.reverse_complement().seq)
                dna=str(subseq.reverse_complement().seq)+dna
            else :
                dna=dna+str(subseq.seq)

    ## dernier gène
    dnaSeq = Seq(dna, IUPAC.ambiguous_dna)
    if(typeseq=="prot") :
        dnaRec = SeqRecord(dnaSeq.translate(), id=sid, description="")
        allseq.append(dnaRec)               
    elif(typeseq=="rna") :
        dnaRec = SeqRecord(dnaSeq.transcribe(), id=sid, description="")
        allseq.append(dnaRec)
    else :
        dnaRec=SeqRecord(dnaSeq, id=sid, description="")
        allseq.append(dnaRec)

    return allseq

# 3. extracting individual exons
def extract_exons(fasta,gff) :
    gff_reader = csv.reader(gff, delimiter='\t')
    
    sid=""
    dna=""
    allseq=[]
    for row in gff_reader :
        if(row[2]=="CDS" or row[2]=="cds") :
            tmp=row[8].split(";") #store first field corresponding to ID=...
            sid=tmp[0][3:] #remove "ID=" from sequence id
            subseq=chr_dict[row[0]][int(row[3])-1:int(row[4])]
            if(row[6]=="-") :
                dna=subseq.reverse_complement().seq
            else :
                dna=subseq.seq
            cds=SeqRecord(dna, id=sid, description="")
            allseq.append(cds)
    
    return allseq
    
# 4. extracting coding sequence : cdna
def extract_cdna(fasta,gff) :
    gff_reader = csv.reader(gff, delimiter='\t')
    
    dna= ""
    sid=""
    allseq=[]
    
    for row in gff_reader :
        if(row[2]=="gene") :
            if(len(dna)>0) :
            # Export sequence
                dnaSeq = Seq(dna, IUPAC.ambiguous_dna)
                dnaRec=SeqRecord(dnaSeq, id=sid, description="")
                allseq.append(dnaRec)
        
            dna = ""
            sid=row[8][3:]
            lastStop=0
        elif(row[2]=="CDS" or row[2]=="cds") :
        # Extract sequence
            # check short intron
            if((int(row[3])-1)<=(lastStop+15)) :
                # integration of the short intron in the sequence
                subseq=chr_dict[row[0]][lastStop:int(row[4])]
            else :
                subseq=chr_dict[row[0]][int(row[3])-1:int(row[4])]
                
            lastStop=int(row[4])
            if(row[6]=="-") :
                dna=str(subseq.reverse_complement().seq)+dna
            else :
                dna=dna+str(subseq.seq)

    ## dernier gène
    dnaSeq = Seq(dna, IUPAC.ambiguous_dna)
    dnaRec=SeqRecord(dnaSeq, id=sid, description="")
    allseq.append(dnaRec)

    return allseq

#5. extracting cdna with frameshift completing with "!"
def extract_frameshift(fasta,gff) :
    gff_reader = csv.reader(gff, delimiter='\t')
    
    dna= ""
    sid=""
    allseq=[]
    lastStop=0
    lastFrame=0
    
    for row in gff_reader :
        if(row[2]=="gene") :
            if(len(dna)>0) :
            # Export sequence
                dnaSeq = Seq(dna, IUPAC.ambiguous_dna)
                dnaRec=SeqRecord(dnaSeq, id=sid, description="")
                allseq.append(dnaRec)
        
            dna = ""
            sid=row[8][3:]
            lastStop=0
            lastFrame=0
        elif(row[2]=="CDS" or row[2]=="cds") :
        # Extract sequence
            # check short intron
            if((int(row[3])-1)<=(lastStop+15)) :
                #check if one or two "!" is needed
                if ((int(row[3])-1-lastStop)%3 == 1):
                    comp="!!"
                elif ((int(row[3])-1-lastStop)%3 == 2) :
                    comp="!"
                # integration of the short intron in the sequence
                subseq=comp+chr_dict[row[0]][lastStop:int(row[4])]
            else :
                subseq=chr_dict[row[0]][int(row[3])-1:int(row[4])]
                
            lastStop=int(row[4])
            if(row[6]=="-") :
                dna=str(subseq.reverse_complement().seq)+dna
            else :
                dna=dna+str(subseq.seq)
            lastFrame=len(dna)%3

    ## dernier gène
    dnaSeq = Seq(dna, IUPAC.ambiguous_dna)
    dnaRec=SeqRecord(dnaSeq, id=sid, description="")
    allseq.append(dnaRec)

    return allseq



#----------------------------------#
#              MAIN
#----------------------------------#

# 1. Importing Fasta Genome 
#============================================= 
chr_dict = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

# 2. Reading GFF and processing
#=============================================
gff=open(args.gff, mode='r')


if(args.type=="gene") :
    #full length gene
    seq_list=extract_gene(chr_dict,gff)
    
elif(args.type=="exon") :
    #indivudual exon
    seq_list=extract_exons(chr_dict,gff)
    
elif(args.type=="cdna") :
    #complete coding sequence with insertions and frameshifts
    seq_list=extract_cdna(chr_dict,gff)
    
elif(args.type=="frameshift") :
    #exact cDNA sequence with frameshift completed with "!"
    seq_list=extract_frameshift(chr_dict,gff)

else :
    #coding sequence as given by CDS feature (without frameshift to allow true protein sequence to be extracted)
    seq_list=extract_coding(chr_dict,gff,args.type)

gff.close()

# 3. Export sequence
#=======================================
SeqIO.write(seq_list, args.output, "fasta")


