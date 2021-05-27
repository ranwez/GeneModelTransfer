# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 09:53:25 2020

@author: gottince
"""

# GFF features representation

class GeneFeatures :
    
    def __init__(self,Identifier,Chromosome,Start,Stop,Strand,Mode) :
        self.id=Identifier
        self.chr=Chromosome
        self.start=Start
        self.stop=Stop
        self.length=0
        self.strand=Strand
        self.mode=Mode
        self.nbCDS=0
        self.cds={}
        self.feature=""
        
    def set_nbCDS(self,nb) :
        self.nbCDS=nb
        
    def get_nbCDS(self) :
        return self.nbCDS

    def set_len(self,ln) :
        self.length=ln

    def get_len(self) :
        return self.length

    def set_stop(self,st) :
        self.stop=st
    
    def get_stop(self) :
        return self.stop

    def set_start(self,st) :
        self.start=st

    def get_start(self) :
        return self.start
    
    def add_feature(self,feat,delim) :
        if(self.feature=="") :
            self.feature=str(feat)
        else :
            self.feature=delim.join([self.feature,feat])
    
    def add_CDS(self,rank,start,stop) :
        self.cds[rank]=CDS(rank,start,stop)
        
    def export(self,filename) :
        with open(filename,mode='a') as File :
            #gene
            line="\t".join([self.chr,self.mode,"gene",str(self.start),str(self.stop),".",self.strand,".",str("ID="+self.id+";comment="+self.feature)])
            File.write(line)
            #mRNA
            line="\t".join([self.chr,self.mode,"mRNA",str(self.start),str(self.stop),".",self.strand,".",str("ID="+self.id+"_mrna;Parent="+self.id)])
            File.write(line)
            #CDS
            for i in range(1,self.nbCDS+1) :
                line="\t".join([self.chr,self.mode,"CDS",str(self.cds[i].start),str(self.cds[i].stop),".",self.strand,".",str("Parent="+self.id+"_mrna")])
                File.write(line)

    def stdexport(self) : ##print
        #gene
        line="\t".join([self.chr,self.mode,"gene",str(self.start),str(self.stop),".",self.strand,".",str("ID="+self.id+";comment="+self.feature)])
        print(line)
        #mRNA
        line="\t".join([self.chr,self.mode,"mRNA",str(self.start),str(self.stop),".",self.strand,".",str("ID="+self.id+"_mrna;Parent="+self.id)])
        print(line)
        #CDS
        for i in range(1,self.nbCDS+1) :
            line="\t".join([self.chr,self.mode,"CDS",str(self.cds[i].start),str(self.cds[i].stop),".",self.strand,".",str("Parent="+self.id+"_mrna")])
            print(line)
        
        
        
class CDS :
    
    def __init__(self,Rank,Start,Stop) :
        self.rank=Rank
        self.start=Start
        self.stop=Stop

    def set_stop(self,st) :
        self.stop=st

    def get_stop(self) :
        return self.stop

    def set_start(self,st) :
        self.start=st

    def get_start(self) :
        return self.start
