#!/usr/bin/env python3
# # -*- coding: utf-8 -*-

"""
Created on Fri Jan 28 13:14:07 2022

@author: gottince
"""
#----------------------------------#
#    DESCRIPTION
#----------------------------------#

# read alignment results and extract candidate regions

#----------------------------------#
#    PARAMETRES
#----------------------------------#
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-g","--gff_file", type=str, help="Exact path of gff file")
parser.add_argument("-t","--table", type=str, help="Exact path of alignment res table")

args = parser.parse_args()


#----------------------------------#
#    SCRIPT
#----------------------------------#


# dict with max intron len
MAX_INTRON = {}

with open(args.gff_file, 'r') as gff :
	for line in gff :
		L=line.split("\t")
		L[3:5]=map(int, L[3:5])
		if (L[2]=="gene"):
			com=L[8].split(";")
			ID=com[0].replace("ID=","")
			newgene=True
		else:
			if (L[2]=="CDS"):
				if newgene :
					stop=int(L[4])
					newgene=False
					MAX_INTRON[ID]=0
				else:
					if(L[3]-stop>MAX_INTRON[ID]):
						MAX_INTRON[ID]=L[3]-stop
						stop=L[4]


# concat consecutive hit
with open(args.table) as res :
	newline=True
	for line in res :
		L=line.replace("\n","").split("\t")
		L[2:13]=map(float, L[2:13])
		if(L[6]<L[7]):
			strand="+"
		else:
			strand="-"
		
		limIntron=5000
		if(MAX_INTRON[L[0]]>limIntron):
			limIntron=MAX_INTRON[L[0]]+500
		
		if newline :
			savedStrand=strand
			savedL=L
			newline=False
		else:
			# compare actual res line (L) with previous one (savedL) and check for merging
			#print(savedL)
			#print(L)
			if(L[0]==savedL[0] and L[1]==savedL[1] and strand==savedStrand and ((strand=="+" and L[5]>savedL[5] and savedL[5]<L[4]+100) or (strand=="-" and L[4]<savedL[4] and L[5]<savedL[4]+100)) and (L[6]-savedL[7]<limIntron or L[7]-savedL[6]<limIntron)):
				## merging
				L[3]=L[3]+savedL[3];
				if(L[4]>savedL[4]):
					L[4]=savedL[4]
				if(L[5]<savedL[5]):
					L[5]=savedL[5]
					
				if(strand=="+"):
					if(L[6]>savedL[6]):
						L[6]=savedL[6]
					if(L[7]<savedL[7]):
						L[7]=savedL[7]
				else:
					if(L[6]<savedL[6]):
						L[6]=savedL[6]
					if(L[7]>savedL[7]):
						L[7]=savedL[7]
					
				L[8]=L[8]+savedL[8]
				L[10]=L[10]+savedL[10]
				L[12]=L[12]+savedL[12]
				L[9]=(L[8]/L[3])*100
				savedL=L
				
			else:
				## not merging
				# check previous res for writting
				if((savedL[5]-savedL[4]+1)/savedL[2]>=0.6 and savedL[9]>=50.0):
					print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(savedL[0],savedL[1],str(savedL[2]),str(savedL[3]),str(savedL[4]),str(savedL[5]),str(savedL[6]),str(savedL[7]),str(savedL[8]),str(savedL[9]),str(savedL[10]),str(savedL[11]),str(savedL[12])))
				
				#update
				savedStrand=strand
				savedL=L
	##last res
	print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(savedL[0],savedL[1],str(savedL[2]),str(savedL[3]),str(savedL[4]),str(savedL[5]),str(savedL[6]),str(savedL[7]),str(savedL[8]),str(savedL[9]),str(savedL[10]),str(savedL[11]),str(savedL[12])))




