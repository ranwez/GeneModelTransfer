#!/usr/bin/env python3
# # -*- coding: utf-8 -*-

"""
Created on Fri Jan 28 13:14:07 2022

@author: gottince
"""
# ----------------------------------#
#    DESCRIPTION
# ----------------------------------#

# read alignment results and extract candidate regions

# ----------------------------------#
#    PARAMETRES
# ----------------------------------#
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-g", "--gff_file", type=str,
                    help="Exact path of gff file")
parser.add_argument("-t", "--table", type=str,
                    help="Exact path of alignment res table")
# qseqid sseqid qlen length qstart qend sstart send nident pident gapopen evalue bitscore

args = parser.parse_args()


# ----------------------------------#
#    FUNCTIONS
# ----------------------------------#

def import_res(resDict, table, ths_inf, ths_sup):
    with open(table, "r") as res:
        ln = 0
        for line in res:
            ln = ln+1
            L = line.replace("\n", "").split("\t")
            L[2:13] = map(float, L[2:13])
            if L[9] >= ths_inf and L[9] <= ths_sup:
                resDict[ln] = L


def max_intron_dict(gff_file):
    MAX_INTRON = {}

    with open(gff_file, 'r') as gff:
        for line in gff:
            if (line.startswith("#")):
                continue
            L = line.strip().split("\t")
            L[3:5] = map(int, L[3:5])
            if (L[2] == "gene"):
                com = L[8].split(";")
                ID = com[0].replace("ID=", "")
                newgene = True
            else:
                if (L[2] == "CDS"):
                    if newgene:
                        stop = int(L[4])
                        newgene = False
                        MAX_INTRON[ID] = 0
                    else:
                        if (L[3]-stop > MAX_INTRON[ID]):
                            MAX_INTRON[ID] = L[3]-stop
                            stop = L[4]

    return MAX_INTRON


def concat_consecutive_hit(resDict):
    newline = True
    L = []
    savedL = []
    cle = ""
    savedKey = ""
    for cle in sorted(resDict):
        L = resDict[cle]
        del resDict[cle]
        if (L[6] < L[7]):
            strand = "+"
        else:
            strand = "-"

        limIntron = 5000
        if (MAX_INTRON[L[0]] > limIntron):
            limIntron = MAX_INTRON[L[0]]+500

        if newline:
            savedStrand = strand
            savedL = L
            savedKey = cle
            newline = False
        else:
            # compare actual res line (L) with previous one (savedL) and check for merging
            # print(savedL)
            # print(L)
            if (L[0] == savedL[0] and L[1] == savedL[1] and strand == savedStrand and ((strand == "+" and L[5] > savedL[5] and savedL[5] < L[4]+100) or (strand == "-" and L[4] < savedL[4] and L[5] < savedL[4]+100)) and (L[6]-savedL[7] < limIntron or L[7]-savedL[6] < limIntron)):
                # merging
                L[3] = L[3]+savedL[3]
                if (L[4] > savedL[4]):
                    L[4] = savedL[4]
                if (L[5] < savedL[5]):
                    L[5] = savedL[5]

                if (strand == "+"):
                    if (L[6] > savedL[6]):
                        L[6] = savedL[6]
                    if (L[7] < savedL[7]):
                        L[7] = savedL[7]
                else:
                    if (L[6] < savedL[6]):
                        L[6] = savedL[6]
                    if (L[7] > savedL[7]):
                        L[7] = savedL[7]

                L[8] = L[8]+savedL[8]
                L[10] = L[10]+savedL[10]
                L[12] = L[12]+savedL[12]
                L[9] = (L[8]/L[3])*100
                savedL = L

            else:
                # not merging
                # check previous res for writting
                if ((savedL[5]-savedL[4]+1.0)/savedL[2] >= 0.6 and savedL[9] >= 50.0):
                    resDict[savedKey] = savedL
                # update
                savedStrand = strand
                savedL = L
                savedKey = cle
    # last res
    resDict[savedKey] = savedL


def remove_overlap(resDict):
    "update resDict, eliminate hits inside other hits"
    newline = True
    L = []
    for cle in sorted(resDict):
        L = resDict[cle]
        # print(L)
        if (L[6] < L[7]):
            strand = "+"
        else:
            strand = "-"
        if newline:
            savedStrand = strand
            savedL = L
            savedKey = cle
            newline = False
        else:
            if L[0] != savedL[0] or L[1] != savedL[1] or strand != savedStrand or (L[6] < L[7] and L[7] > savedL[7]) or (L[6] > L[7] and L[6] > savedL[6]):
                # print("keepline")
                savedStrand = strand
                savedL = L
                savedKey = cle
            else:
                # print("del line")
                del resDict[cle]


def print_res(resDict):
    for cle in sorted(resDict):
        L = resDict[cle]
        # print(str(cle))
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(L[0], L[1], str(L[2]), str(L[3]), str(
            L[4]), str(L[5]), str(L[6]), str(L[7]), str(L[8]), str(L[9]), str(L[10]), str(L[11]), str(L[12])))


# ----------------------------------#
#    MAIN
# ----------------------------------#


# dict with max intron len per query
MAX_INTRON = max_intron_dict(args.gff_file)

# import result table
RES = {}

# process
import_res(RES, args.table, 65, 100)
concat_consecutive_hit(RES)
import_res(RES, args.table, 45, 65)
remove_overlap(RES)
concat_consecutive_hit(RES)
print_res(RES)
