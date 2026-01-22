"""
Created on Tue Feb  1 11:45:02 2022

@author: gottince

@description: script permettant dde formatter un GFF.
"""
import sys
import os
import argparse
import GFFclass
import csv

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from CANDIDATE_LOCI.gff_utils import parse_and_sort_gff

#----------------------------------#
#    PARAMETRES
#----------------------------------#
# - File name of gff
# - output filename

parser = argparse.ArgumentParser()

parser.add_argument("-g","--gff", type=str, help="Exact path of gff file")
parser.add_argument("-o","--output", type=str, help="Output gff file name")

args = parser.parse_args()

#----------------------------------#
#            FUNCTIONS
#----------------------------------#

def importGFF(gff) :
    gff_reader = csv.reader(gff,delimiter = "\t")

    genes=[]
    gnCount=0

    for row in gff_reader :        
        if len(row)==9 :            
            if(row[2]=="gene") :
                gnCount+=1
                ident=""
                L=row[8].split(";")                
                for ft in L :
                    if "ID=" in ft :
                        ident=ft.replace("ID=","")
                genes.append(GFFclass.GeneFeatures(ident,row[0],int(row[3]),int(row[4]),row[6],row[1]))                                
                genes[gnCount-1].add_feature(row[8],";")               

                
            elif row[2]=="CDS" :
                rk=genes[gnCount-1].get_nbCDS()
                genes[gnCount-1].add_CDS((rk+1),int(row[3]),int(row[4]))
                genes[gnCount-1].set_nbCDS(rk+1)

                
            elif row[2]=="exon" :  
                rk=genes[gnCount-1].get_nbExon()
                genes[gnCount-1].add_Exon((rk+1),int(row[3]),int(row[4]))
                genes[gnCount-1].set_nbExon(rk+1)
       
    return genes




#----------------------------------#
#              MAIN
#----------------------------------#

gff=parse_and_sort_gff(args.gff)
myGenes=importGFF(gff) 
#myGenes=importGFF("C:/Users/gottince/Documents/DATA/ORYZA/GELOC/release_2/NIPPONBARE/Oryza_Nipponbare_IRGSP-1.0_LRR-CR__20210715.gff")


for ign in range(len(myGenes)) : ## for each gene
    if not myGenes[ign].eval_features() :
        myGenes[ign].predict_exon()
    myGenes[ign].predict_sequence_alteration()
    myGenes[ign].export(args.output) 
 
# When no genes are found create an empty file to respect implicit contract
if len(myGenes) == 0:
    with open(args.output, 'w') as f:
        pass  # 'pass' simply does nothing, leaving the file empty
       
