# -*- coding: utf-8 -*-
#!/usr/bin/python3
"""
Created on Mon Feb 17 11:38:07 2020

@author: gottince
"""

#----------------------------------#
#    DESCRIPTION
#----------------------------------#


# Read GFF file from exonerate output and test for canonical intron splice site (GT-AG sequence)
# Correct intron position if a canonical site (GT-AG) is found at 1 or 2 bases from current position
# Check for gene to be of length *3. If not, add nucleotide at last position to correct the length


#----------------------------------#
#    PARAMETRES
#----------------------------------#
# - File name of input genomic fasta
# - File name of input Table gene file 
# - File name of output Table gene file (with corrected intron)

import argparse
import csv
import GFFclass
from Bio import SeqIO



parser = argparse.ArgumentParser()

parser.add_argument("-f","--fasta", type=str, help="Exact path of genomic fasta file")
parser.add_argument("-g","--gff", type=str, help="Exact path of feature gff file")

args = parser.parse_args()


#----------------------------------#
#            FUNCTIONS
#----------------------------------#

def importGFF(myfile) :

    gff = open(myfile, mode = "r")
    gff_reader = csv.reader(gff,delimiter = "\t")

    genes=[]
    gnCount=0

    for row in gff_reader :
        if(len(row)==9) :
            if(row[2]=="gene") :
                feat=""
                ident=""
                L=row[8].split(";")
                for ft in L :
                    if("ID=" in ft) :
                        #ft.replace("ID=","")
                        ident=ft.replace("ID=","")
                    if("origin=" in ft) :
                        feat=ft.replace("origin=","Origin:")
                #ident=str(row[1]+"_"+row[0]+"_"+row[3])
                genes.append(GFFclass.GeneFeatures(ident,row[0],int(row[3]),int(row[4]),row[6],row[1]))
                if(len(feat)>0) :
                    genes[gnCount].add_feature(feat," / ")
                gnCount+=1
                #genes[gnCount].add_feature("pred:MappingCDS",";")
            elif(row[2]=="CDS") :
                rk=genes[gnCount-1].get_nbCDS()
                genes[gnCount-1].add_CDS((rk+1),int(row[3]),int(row[4]))
                genes[gnCount-1].set_nbCDS(rk+1)
                size=genes[gnCount-1].get_len()+(int(row[4])-int(row[3])+1)
                genes[gnCount-1].set_len(size)
            #else : ##mRNA or other
            ## feature si section comment
            L=row[8].split(";")
            for ft in L :
                if("comment=" in ft) :
                    ft.replace("comment=","",1)
                    genes[gnCount-1].add_feature(ft," / ")
    
    
    gff.close()
    
    return genes

def isCanonical_forward(DNA_dict,Chr,startIntron,stopIntron) :
    "Check if donnor and acceptor splice sites are canonical in forward strand"
    donnor=DNA_dict[Chr][startIntron:startIntron+2] ##GT
    acceptor=DNA_dict[Chr][stopIntron-3:stopIntron-1] ##AG
    #print("donnor=",donnor.seq," acceptor=",acceptor.seq)
    #if((donnor.seq=="GT" or donnor.seq=="GC") and acceptor.seq=="AG") :
    if(donnor.seq=="GT" and acceptor.seq=="AG") :
        return True
    else :
        return False
    
def findCanonical_forward(DNA_dict,Chr,startIntron,stopIntron) :
    "Look for canonical splice site at -2, -1, +1 or +2 from actual intron position in forward strand."
    if(isCanonical_forward(DNA_dict,Chr,startIntron-2,stopIntron-2)) :
        return [-2,-2]
    elif(isCanonical_forward(DNA_dict,Chr,startIntron-2,stopIntron+1)) :
        return [-2,1]
    elif(isCanonical_forward(DNA_dict,Chr,startIntron-1,stopIntron-1)) :
        return [-1,-1]
    elif(isCanonical_forward(DNA_dict,Chr,startIntron-1,stopIntron+2)) :
        return [-1,+2]
    elif (isCanonical_forward(DNA_dict,Chr,startIntron+1,stopIntron+1)):
        return [1,1]
    elif(isCanonical_forward(DNA_dict,Chr,startIntron+1,stopIntron-2)) :
        return [1,-2]
    elif(isCanonical_forward(DNA_dict,Chr,startIntron+2,stopIntron+2)) :
        return [2,2]
    elif(isCanonical_forward(DNA_dict,Chr,startIntron+2,stopIntron-1)) :
        return [2,-1]
    else :
        return 0


def isCanonical_reverse(DNA_dict,Chr,startIntron,stopIntron) :
    "Check if donnor and acceptor splice sites are canonical in reverse strand"
    donnor=DNA_dict[Chr][stopIntron-3:stopIntron-1] ##AC
    acceptor=DNA_dict[Chr][startIntron:startIntron+2] ##CT
    #print("donnor=",donnor.seq," acceptor=",acceptor.seq)
    #if((donnor.seq=="AC" or donnor.seq=="CG") and acceptor.seq=="CT") :
    if(donnor.seq=="AC" and acceptor.seq=="CT") :
        return True
    else :
        return False    

def findCanonical_reverse(DNA_dict,Chr,startIntron,stopIntron) :
    "Look for canonical splice site at -2, -1, +1 or +2 from actual intron position in reverse strand."
    if(isCanonical_reverse(DNA_dict,Chr,startIntron-2,stopIntron-2)) :
        return [-2,-2]
    elif(isCanonical_reverse(DNA_dict,Chr,startIntron-2,stopIntron+1)) :
        return [-2,1]
    elif(isCanonical_reverse(DNA_dict,Chr,startIntron-1,stopIntron-1)) :
        return [-1,-1]
    elif(isCanonical_reverse(DNA_dict,Chr,startIntron-1,stopIntron+2)) :
        return [-1,+2]
    elif (isCanonical_reverse(DNA_dict,Chr,startIntron+1,stopIntron+1)):
        return [1,1]
    elif(isCanonical_reverse(DNA_dict,Chr,startIntron+1,stopIntron-2)) :
        return [1,-2]
    elif(isCanonical_reverse(DNA_dict,Chr,startIntron+2,stopIntron+2)) :
        return [2,2]
    elif(isCanonical_reverse(DNA_dict,Chr,startIntron+2,stopIntron-1)) :
        return [2,-1]
    else :
        return 0


def isStop(DNA_dict,Chr,stop, strand) :
    if(strand=="+") :
        seq=DNA_dict[Chr][stop-3:stop].seq
        #print(seq)
        if(seq=="TAA" or seq=="TAG" or seq=="TGA"):
            return True
        else :
            return False
    else :
        seq=DNA_dict[Chr][stop-1:stop+2].seq
        #print(seq)
        if(seq=="TTA" or seq=="CTA" or seq=="TCA"):
            return True
        else :
            return False
        

def findStop_forward(DNA_dict,Chr,posInit,distanceMax) :
    for i in range(1,distanceMax+1) :
        pos=posInit+i*3
        if(isStop(DNA_dict,Chr,pos,"+")) :
               return pos
    return 0

def findStop_reverse(DNA_dict,Chr,posInit,distanceMax) :
    for i in range(1,distanceMax+1) :
        pos=posInit-i*3
        if(isStop(DNA_dict,Chr,pos,"-")) :
               return pos
    return 0

def isStart(DNA_dict,Chr,start,strand):
    if(strand=="+") :
        seq=DNA_dict[Chr][start-1:start+2].seq
        #print(seq)
        if(seq=="ATG"):
            return True
        else :
            return False
    else :
        seq=DNA_dict[Chr][start-3:start].seq
        #print(seq)
        if(seq=="CAT"):
            return True
        else :
            return False

def findStart_forward(DNA_dict,Chr,posInit,distanceMax) :
    for i in range(1,distanceMax+1) :
        pos=posInit-i*3
        if(isStart(DNA_dict,Chr,pos,"+")) :
               return pos
    return 0

def findStart_reverse(DNA_dict,Chr,posInit,distanceMax) :
    for i in range(1,distanceMax+1) :
        pos=posInit+i*3
        if(isStart(DNA_dict,Chr,pos,"-")) :
               return pos
    return 0


#----------------------------------#
#              MAIN
#----------------------------------#

## Read genome
chr_dict = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

## read gff
#gff=open(args.gff, mode='r')
#gff_reader = csv.reader(gff, delimiter='\t')

## importing data
myGenes=importGFF(args.gff)

for ign in range(len(myGenes)) : ## for each gene
    toCheck=False
    correction=False
    modifstop=False
    modifstart=False
    myProt=myGenes[ign].id
    strand=myGenes[ign].strand
    #print(myProt+" ; "+strand)
    #Chr='_'.join(row[0].split("_")[1:-1]) ## Chromosome Id
    Chr=myGenes[ign].chr
    ns1=0
    ns2=0
    

    ##Checking for CDS length (should be a *3)
    if(not myGenes[ign].length%3==0) :    
        modifstop=True
        #print("->pb len ;"+str(myGenes[ign].length)+" "+str(myGenes[ign].length%3))
        #Change stop
        if(strand=="+") :
            #Os=myGenes[ign].get_stop()
            myGenes[ign].set_stop(myGenes[ign].get_stop()+(3-myGenes[ign].length%3)) ## add missing nucleotides
            myGenes[ign].cds[myGenes[ign].nbCDS].set_stop(myGenes[ign].cds[myGenes[ign].nbCDS].get_stop()+(3-myGenes[ign].length%3))
        else :
            #Os=myGenes[ign].get_start()
            myGenes[ign].set_start(myGenes[ign].get_start()-(3-myGenes[ign].length%3))
            myGenes[ign].cds[1].set_start(myGenes[ign].cds[1].get_start()-(3-myGenes[ign].length%3))
            #print(str("     modif stop ; oldstop "+str(Os)+" ; newstop "+str(myGenes[ign].get_start())))

    ##Checking for start and stop
    if(strand=="+") :
        #start
        if(not isStart(chr_dict,Chr,myGenes[ign].start,strand)):
            #Ost=myGenes[ign].get_start()
            #print("--pb start ;")
            ns1=findStart_forward(chr_dict,Chr,myGenes[ign].start,20)
            if(ns1!=0) :
                modifstart=True
                myGenes[ign].set_start(ns1)
                myGenes[ign].cds[1].set_start(ns1) 
            #print(str("     modif start ; oldstart "+str(Ost)+" ; newstop "+str(ns1)))
        #stop
        #Os=myGenes[ign].get_stop()
        if(not isStop(chr_dict,Chr,myGenes[ign].stop,strand)):
            #Os=myGenes[ign].get_stop()
            #print("--pb stop ;")
            ns2=findStop_forward(chr_dict,Chr,myGenes[ign].stop,20)
            if(ns2!=0) :
                modifstop=True
                myGenes[ign].set_stop(ns2)
                myGenes[ign].cds[myGenes[ign].nbCDS].set_stop(ns2) 
            #print(str("     modif stop ; oldstop "+str(Os)+" ; newstop "+str(ns2)))
    else :
        #start
        if(not isStart(chr_dict,Chr,myGenes[ign].stop,strand)):
            #Ost=myGenes[ign].get_stop()
            #print("--pb start ;")
            ns1=findStart_reverse(chr_dict,Chr,myGenes[ign].stop,20)
            if(ns1!=0) :
                modifstart=True
                myGenes[ign].set_stop(ns1)
                myGenes[ign].cds[myGenes[ign].nbCDS].set_stop(ns1)
            #print(str("     modif start ; oldstart "+str(Ost)+" ; newstop "+str(ns1)))
        #stop
        if(not isStop(chr_dict,Chr,myGenes[ign].start,strand)):
            #Os=myGenes[ign].get_start()
            #print("--pb stop ;")
            ns2=findStop_reverse(chr_dict,Chr,myGenes[ign].start,20)
            if(ns2!=0) :
                modifstop=True
                myGenes[ign].set_start(ns2)
                myGenes[ign].cds[1].set_start(ns2)
            #print(str("     modif stop ; oldstop "+str(Os)+" ; newstop "+str(ns2)))


    ##Checking for intron
    if(myGenes[ign].nbCDS>=2) :
        if(strand=="+") :
            #for i in range(3,len(row)-2,2) :
            for icds in range(1,myGenes[ign].nbCDS):
                ## si pas frameshift et pas canonique
                #if(int(row[i+1])>int(row[i])+10 and not isCanonical_forward(chr_dict,Chr,int(row[i]),int(row[i+1]))) :
                if( (myGenes[ign].cds[icds+1].start>myGenes[ign].cds[icds].stop+20) and not isCanonical_forward(chr_dict,Chr,myGenes[ign].cds[icds].stop,myGenes[ign].cds[icds+1].start) ) :
                    #decal=findCanonical_forward(chr_dict,Chr,int(row[i]),int(row[i+1]))
                    decal=findCanonical_forward(chr_dict,Chr,myGenes[ign].cds[icds].stop,myGenes[ign].cds[icds+1].start)
                    if(decal!=0) : #if(decal==0):
                        #toCheck=True
                        ## non canonical splicing site and not corrected --> to check
                        ## keep old positions
                        #newLine=str(newLine+";"+str(row[i])+";"+str(row[i+1]))
                    #else :
                        ## non canonical splicing site and corrected
                        ## Change position
                        correction=True
                        myGenes[ign].cds[icds].set_stop(myGenes[ign].cds[icds].stop+decal[0])
                        myGenes[ign].cds[icds+1].set_start(myGenes[ign].cds[icds+1].start+decal[1])
                        #newLine=str(newLine+";"+str(nstart)+";"+str(nstop))
                #else :
                    ## write old
                    #newLine=str(newLine+";"+str(row[i])+";"+str(row[i+1]))
            ## last position
            #newLine=str(newLine+";"+row[i+2])
        else :
            #...CT--->AC... forward
            #...GA<---TG... reverse        
            #for i in range(3,len(row)-2,2) :
            for icds in range(1,myGenes[ign].nbCDS):
                #if(int(row[i+1])>int(row[i])+10 and not isCanonical_reverse(chr_dict,Chr,int(row[i]),int(row[i+1]))) :
                if( (myGenes[ign].cds[icds+1].start>myGenes[ign].cds[icds].stop+20) and not isCanonical_reverse(chr_dict,Chr,myGenes[ign].cds[icds].stop,myGenes[ign].cds[icds+1].start) ) :
                    #decal=findCanonical_reverse(chr_dict,Chr,int(row[i]),int(row[i+1]))
                    decal=findCanonical_reverse(chr_dict,Chr,myGenes[ign].cds[icds].stop,myGenes[ign].cds[icds+1].start)
                    if(decal==0) :
                        toCheck=True
                        ## site non canonique et pas corrigé --> a voir
                        ## write old
                        #newLine=str(newLine+";"+str(row[i])+";"+str(row[i+1]))
                    else :
                        ## site non canonique et corrigé
                        ## write new
                        correction=True
                        myGenes[ign].cds[icds].set_stop(myGenes[ign].cds[icds].stop+decal[0])
                        myGenes[ign].cds[icds+1].set_start(myGenes[ign].cds[icds+1].start+decal[1])
                        #nstart=int(row[i])+decal[0]
                        #nstop=int(row[i+1])+decal[1]
                        #newLine=str(newLine+";"+str(nstart)+";"+str(nstop))
                #else :
                    ##write old
                    #newLine=str(newLine+";"+str(row[i])+";"+str(row[i+1]))
            ## last position
            #newLine=str(newLine+";"+row[i+2])
        
    #if(toCheck) :
        #print(newLine," toCheck")
        #myGenes[ign].add_feature("exo_corr:non_canonical_intron"," / ")
    if(correction) :
        #print(newLine," corrected")
        myGenes[ign].add_feature("exo_corr:corrected_intron"," / ")
    if(modifstop) :
        myGenes[ign].add_feature("exo_corr:modif_stop"," / ")
    if(modifstart):
        myGenes[ign].add_feature("exo_corr:modif_start"," / ")    
    if(not correction and not modifstop and not modifstart) :
        myGenes[ign].add_feature("exo_corr:NA"," / ")
    #else :
        #print(newLine)
    

    ##export
    myGenes[ign].stdexport()


## Fin 
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    