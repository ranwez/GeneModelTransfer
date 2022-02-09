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

parser.add_argument("-t","--table", type=str, help="Exact path of alignment res table")
parser.add_argument("-o","--outgff", type=str, help="Exact path of output regions in gff")

args = parser.parse_args()


#----------------------------------#
#    FUNCTIONS
#----------------------------------#

def import_res(resDict, table, ths_inf=0, ths_sup=100) :
	with open(table, "r") as res :
		ln=0
		for line in res :
			ln=ln+1
			L=line.replace("\n","").split("\t")
			L[2:13]=map(float, L[2:13])
			if L[9]>=ths_inf and L[9]<=ths_sup :
				resDict[ln]=L


def select_best_query_per_region(resDict) :
	newline=True
	L=[]
	savedL=[]
	cle=""
	savedKey=""
	for cle in sorted(resDict):
		L=resDict[cle]
		del resDict[cle]
		if newline :
			savedL=L
			savedKey=cle
			newline=False
		else:
			if L[1]==savedL[1] and (L[6]<savedL[7] or savedL[6]>L[7]) :
				if savedL[12]<L[12] :
					#si score actuel est meilleur que le précédent
					#ont ne garde pas la ligne precedente
					savedL=L
					savedKey=cle 
			else :
				resDict[savedKey]=savedL
				savedL=L
				savedKey=cle
	resDict[savedKey]=savedL

def max_boundaries_dict(resDict) :
	"function that extract 5' and 3' limit of each regions according to the previous and next one"
	"to avoid overlapping"
	BOUND = {}
	START = {}
	STOP = {}
	chr=""
	savedKey = -1
	for cle in resDict :
		L=resDict[cle]
		if L[6]<L[7] : #strand +
			start=L[6]
			stop=L[7]
		else :
			start=L[7]
			stop=L[6]
		
		#Overlap STOP precedent
		if L[1]==chr and STOP[savedKey]>start : ## si meme chromosome et stop precedent chevauche la région actuelle
			STOP[savedKey]=start-1

		## START actuel
		if L[4]>10 : # si plus de 10 base manquante en 5'
			START[cle]=start-3000
		else :
			START[cle]=start-300
		
		# Overlap START actuel
		if savedKey>=0 and START[cle]<STOP[savedKey] :
			START[cle]=STOP[savedKey]+1
		
		##STOP actuel
		if L[5]<L[2]-10 : # si plus de 10 base manquante en 3'
			STOP[cle]=stop+3000
		else : 
			STOP[cle]=stop+300
		
		## SAVE
		BOUND[cle]=(START[cle],STOP[cle])
		savedKey=cle
		chr=L[1]
	
	START.clear()
	STOP.clear()
	return BOUND

def print_gff_regions(resDict, limitDict, outfile) :
	with open(outfile, "w") as gff :
		for cle in sorted(resDict) :
			L=resDict[cle]
			B=limitDict[cle]
			if L[6]<L[7] : #strand +
				ident="{}_{:0>8d}".format(L[1], int(B[0]))
				line=L[1]+"\tLRRtransfer\tgene\t"+str(int(B[0]))+"\t"+str(int(B[1]))+"\t.\t+\t.\tID="+ident+";origin="+L[0]+"\n"
				gff.write(line)
				print('{}\t{}\t{}'.format(ident,L[0],"+"))
			else : # starnd -
				ident="{}_{:0>8d}".format(L[1], int(B[1]))
				line=L[1]+"\tLRRtransfer\tgene\t"+str(int(B[0]))+"\t"+str(int(B[1]))+"\t.\t-\t.\tID="+ident+";origin="+L[0]+"\n"
				gff.write(line)
				print('{}\t{}\t{}'.format(ident,L[0],"-"))


#----------------------------------#
#    MAIN
#----------------------------------#

# import result table
RES={}

import_res(RES, args.table)
select_best_query_per_region(RES)

#Limites des zones en 5' et 3'
LIMITS = max_boundaries_dict(RES)

print_gff_regions(RES, LIMITS, args.outgff)
