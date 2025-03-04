# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 15:15:02 2020

@author: gottince

@description: script permettant d'extraire des sequences proteiques et
nucleotidiques d'un fasta genomique a partir d'un gff.
"""
print("Starting Extract_sequences_from_genome.py", flush=True)
import argparse
import csv
import tempfile
import os
from Bio import SeqIO
from Bio.Seq import translate, reverse_complement
#from CANDIDATE_LOCI.gff_utils import parse_gff, sort_gff
from sort_gff import sort_and_write_gff

#----------------------------------#
#    PARAMETRES
#----------------------------------#
# - File name of genomic fasta
# - output filename
# - output sequence type dna/rna/protein

def positive_int(value):
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError(f"{value} is an invalid positive int value")
    return ivalue

parser = argparse.ArgumentParser()

parser.add_argument("-f","--fasta", type=str, help="Exact path of genomic fasta file")
parser.add_argument("-g","--gff", type=str, help="Exact path of gff file")
parser.add_argument("-o","--output", type=str, help="Output file name")
parser.add_argument("-t","--type", type=str, choices=["gene","cdna","prot","FScdna","FSprot","exon"], help="type of extracting sequence : <gene> extract genomique dna sequence with intron, <cdna> extract reconstruct coding nucleotide sequence with frameshift, <prot> extract reconstruct protein sequence, <exons> extract isolated coding exons, <cds> extract reconstruct coding nucleotide sequence without frameshift.")
parser.add_argument("-m", "--margin", type=positive_int, default=0, help="Margin size (positive integer) to include flanking regions around gene extracted sequences. Default is 0.")
args = parser.parse_args()

#----------------------------------#
#            FUNCTIONS
#----------------------------------#

def should_skip_line(row):
    line = "".join(row)
    return len(line.strip()) == 0 or line.strip().startswith('#')


# 1. extracting complete gene
def extract_gene(fasta,gff,margin=0) :
    gff_reader = csv.reader(gff, delimiter='\t')
    print("--- [extract_gene] gff read", flush=True)
    sid=""
    allseq=[]
    for row in gff_reader :
        if should_skip_line(row):
            continue
        if(row[2]=="gene") :
            tmp=row[8].split(';')
            sid=tmp[0][3:]

            begin = int(row[3]) - 1 - margin
            begin = max(0, begin)  # Ensure begin is not negative
            end = int(row[4]) + margin
            end = min(end, len(str(chr_dict[row[0]].seq)))  # Ensure end does not exceed sequence length

            subseq = chr_dict[row[0]][begin:end]  # Extract sequence with margin
            #subseq=chr_dict[row[0]][int(row[3])-1:int(row[4])]
            if(row[6]=="-") :
                dna=reverse_complement(str(subseq.seq))
                allseq.append((sid,dna))
            else :
                dna=str(subseq.seq)
                allseq.append((sid,dna))

    return allseq


# 2. extracting coding sequence : prot and cdna without FS recoded
## the protein sequence should end at the first stop codon
def extract_coding(fasta,gff,typeseq) :
    gff_reader = csv.reader(gff, delimiter='\t')

    dna= ""
    sid=""
    allseq=[]

    for row in gff_reader :
        if should_skip_line(row):
            continue
        if(row[2]=="gene") :
            if(len(dna)>0) :
            # Export sequence
                if(typeseq=="prot") :
                    #dnaRec = SeqRecord(Seq(dna).translate(to_stop=True)+"*", id=sid, description="")
                    prot=translate(dna, to_stop=True)+"*"
                    allseq.append((sid,prot))
                else :
                    allseq.append((sid,dna))

            dna = ""
            tmp=row[8].split(';')
            sid=tmp[0][3:]
        elif(row[2]=="exon") :
        # Extract sequence
            subseq=chr_dict[row[0]][int(row[3])-1:int(row[4])]
            if(row[6]=="-") :
                #dna=dna+str(subseq.reverse_complement().seq)
                dna=str(subseq.reverse_complement().seq)+dna
            else :
                dna=dna+str(subseq.seq)

    ## dernier gene
    if(typeseq=="prot") :
        #dnaRec = SeqRecord(Seq(dna).translate(to_stop=True)+"*", id=sid, description="")
        prot=translate(dna, to_stop=True)+"*"
        allseq.append((sid, prot))
    else :
        allseq.append((sid,dna))

    return allseq



# 3. extracting individual cds fragment
def extract_exons(fasta,gff) :
    gff_reader = csv.reader(gff, delimiter='\t')

    sid=""
    dna=""
    allseq=[]
    for row in gff_reader :
        if should_skip_line(row):
            continue
        if(row[2] == "gene"):
            tmp=row[8].split(';')
            gid=tmp[0][3:]
        if(row[2]=="CDS" or row[2]=="cds") :
            tmp=row[8].split(";") #store first field corresponding to ID=...
            sid=tmp[0][3:] #remove "ID=" from sequence id
            subseq=chr_dict[row[0]][int(row[3])-1:int(row[4])]
            if(row[6]=="-") :
                dna=str(subseq.reverse_complement().seq)
            else :
                dna=str(subseq.seq)
            allseq.append((gid+"_"+sid,dna))

    return allseq


#4. extracting cdna or prot with frameshift completing with "!"
def extract_frameshift(fasta, gff, typeseq) :
    gff_reader = csv.reader(gff, delimiter='\t')

    dna = ""
    sid = ""
    allseq = []  ## list de tuple (id, seq)
    lastStop = 0

    for row in gff_reader :
        if should_skip_line(row):
            continue
        if(row[2]=="gene") :
            if(len(dna)>0) :
            # Export sequence
                if(typeseq=="FSprot"):
                    if(len(dna) %3 != 0):
                        print("Warning: Partial codon, len(sequence) not a multiple of three:")
                        print(sid)
                        print(dna)
                    ## change "!" to "X" to allow auto translation
                    prot = translate(dna.replace("!","N"))
                    allseq.append((sid,prot))
                else:
                    allseq.append((sid,dna))

            dna = ""
            tmp = row[8].split(';')
            sid = tmp[0][3:]
            lastStop = 0
        elif(row[2]=="CDS" or row[2]=="cds") :
        # Extract sequence
            # mark short intron with frameshift (could lead to 1 or 2 ! in AA translation AT! !!G ATG or ATG !!! ATG)
            if((int(row[3])-1)<=(lastStop+15) and len(dna)>0 ) :
                comp="!!!"
            else:
                comp=""

            subseq=chr_dict[row[0]][int(row[3])-1:int(row[4])].upper()
            lastStop=int(row[4])
            if(row[6]=="-") :
                dna=reverse_complement(str(subseq.seq))+comp+dna
            else :
                dna=dna+comp+str(subseq.seq)
            #print (dna)

    ## dernier gene
    if(typeseq=="FSprot"):
        ## change "!" to "N" to allow auto translation
        prot = translate(dna.replace("!","N"))
        allseq.append((sid,prot))
    else:
        allseq.append((sid,dna))

    return allseq

def write_fasta(mylist, myfile):
    if myfile:  # Check if filename is provided
        fastafile = open(myfile, mode='w')
        for seq in mylist:
            line = ">"+seq[0]+"\n"+seq[1]+"\n"
            fastafile.write(line)
        fastafile.close()
    else:  # If filename is not provided, print to console
        for seq in mylist:
            print(">"+seq[0])
            print(seq[1])
    print(f"--- output fasta file written to: {myfile}", flush=True)


#----------------------------------#
#              MAIN
#----------------------------------#

# 1. Importing Fasta Genome
#=============================================
chr_dict = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

# 2. Reading GFF and processing
#=============================================


with tempfile.NamedTemporaryFile(suffix=".gff", delete=False) as temp_file:
    temp_filename = temp_file.name  # Get the path of the temp file
    try:
        # ensure the gff is sorted so that features of the same gene are grouped
        sort_and_write_gff(args.gff, temp_filename)
        
        gff=open(temp_filename, mode='r')

        seq_list=[]

        if(args.type=="gene") :
            #full length gene
            seq_list=extract_gene(chr_dict,gff,args.margin)

        elif(args.type=="exon") :
            #indivudual exon
            seq_list=extract_exons(chr_dict,gff)

        elif(args.type=="cdna" or args.type=="prot") :
            #complete coding sequence with insertions and frameshifts
            seq_list=extract_coding(chr_dict, gff, args.type)

        elif(args.type=="FScdna" or args.type=="FSprot") :
            #exact cDNA or protein sequence with frameshift completed with "!"
            seq_list=extract_frameshift(chr_dict, gff, args.type)

        else :
            print("argument error : unknown type -t "+args.type)

        gff.close()
        

    finally:
        # Ensure cleanup
        if os.path.exists(temp_filename):
            os.remove(temp_filename)  # Delete the temp file



# 3. Export sequence
#=======================================
write_fasta(seq_list, args.output)
