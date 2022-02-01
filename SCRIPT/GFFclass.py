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
        self.cds={} ## dict of SeqFeature
        self.nbExon=0
        self.Exon={} ## dict of SeqFeature
        self.feature=""
        
    def set_nbCDS(self,nb) :
        self.nbCDS=nb
        
    def get_nbCDS(self) :
        return self.nbCDS
    
    def set_nbExon(self,nb) :
        self.nbExon=nb
        
    def get_nbExon(self) :
        return self.nbExon

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
        ## Allow to modify information from the 9th column
        if(self.feature=="") :
            self.feature=str(feat)
        else :
            self.feature=delim.join([self.feature,feat])
    
    def add_CDS(self,rank,start,stop) :
        self.cds[rank]=SeqFeature(rank,start,stop)
        
    def add_Exon(self,rank,start,stop) :
        self.Exon[rank]=SeqFeature(rank,start,stop)
        
    def concat_CDS(self, indexExon1, indexExon2) :
        # fusion de deux exon si même phase et pas de stop entre
        stop = self.cds[indexExon2].get_stop()
        self.cds[indexExon1].set_stop(stop)
        # remove exon2, reorder exon
        if(indexExon2<self.nbCDS) :
            for rk in range(indexExon2+1,self.nbCDS+1) :
                self.cds[rk-1]=self.cds[rk]
        self.cds.pop(self.nbCDS)
        self.set_nbCDS(len(self.cds))
        
        
    def eval_features(self) :
        # return True if Exons and CDS are matching, False otherwise
        correspStart = False
        correspStop = False
        
        if len(self.cds)==0 or len(self.Exon)==0 :
            return False
        else :
            for exon in self.Exon :
                for cds in self.cds :
                    if exon.get_start() == cds.get_start() :
                        correspStart=True
                    if exon.get_stop() == cds.get_stop() :
                        correspStop=True
                if not correspStart or not correspStop :
                    return False
        return True
                    

    def predict_exon(self) :
        # Determine Exon features from CDS (concatenate CDS separated by frameshift)
        self.add_Exon(rank=1, start=self.cds[1].get_start(), stop=self.cds[1].get_stop())
        rkexon=1
        for rkcds in range(2,len(self.cds)+1) :
            if self.cds[rkcds].get_start()<self.Exon[rkexon].get_stop()+30 :
                #si CDS 2 trop proche de cds 1 -> FS
                # alors update exon précédent
                self.Exon[rkexon].set_stop(self.cds[rkcds].get_stop())
            else :
                # sinon ajout nouvel exon
                rkexon+=1
                self.add_Exon(rank=rkexon, start=self.cds[rkcds].get_start(), stop=self.cds[rkcds].get_stop())
        self.set_nbExon(len(self.Exon))
       
        
    def export(self,filename) :
        with open(filename,mode='a') as File :
            #gene
            line="\t".join([self.chr,self.mode,"gene",str(self.start),str(self.stop),".",self.strand,".",str(self.feature)])
            File.write(line+"\n")
            #mRNA
            line="\t".join([self.chr,self.mode,"mRNA",str(self.start),str(self.stop),".",self.strand,".",str("ID="+self.id+"_mrna;Parent="+self.id)])
            File.write(line+"\n")
            #exon
            for i in range(1,len(self.Exon)+1) :
                j=i
                if (self.strand=="-") :
                    j=len(self.Exon)-i+1
                line="\t".join([self.chr,self.mode,"exon",str(self.Exon[i].start),str(self.Exon[i].stop),".",self.strand,".",str("ID="+self.id+":exon_"+str(j)+";Parent="+self.id+"_mrna")])
                File.write(line+"\n")   
            #CDS
            for i in range(1,self.nbCDS+1) :
                j=i
                if self.strand=="-" :
                    j=len(self.cds)-i+1
                line="\t".join([self.chr,self.mode,"CDS",str(self.cds[i].start),str(self.cds[i].stop),".",self.strand,".",str("ID="+self.id+":cds_"+str(j)+";Parent="+self.id+"_mrna")])
                File.write(line+"\n")



    def stdexport(self) : ##print
        #gene
        line="\t".join([self.chr,self.mode,"gene",str(self.start),str(self.stop),".",self.strand,".",str("ID="+self.id+";"+self.feature)])
        print(line)
        #mRNA
        line="\t".join([self.chr,self.mode,"mRNA",str(self.start),str(self.stop),".",self.strand,".",str("ID="+self.id+"_mrna;Parent="+self.id)])
        print(line)
        #exon
        for i in range(1,self.nbExon+1) :
            line="\t".join([self.chr,self.mode,"exon",str(self.Exon[i].start),str(self.Exon[i].stop),".",self.strand,".",str("ID="+self.id+":exon"+str(i)+";Parent="+self.id+"_mrna")])
            print(line)               
        #CDS
        for i in range(1,self.nbCDS+1) :
            line="\t".join([self.chr,self.mode,"CDS",str(self.cds[i].start),str(self.cds[i].stop),".",self.strand,".",str("Parent="+self.id+"_mrna")])
            print(line)
        
        

        
class SeqFeature :
    
    def __init__(self,Rank,Start,Stop) :
        self.rank=Rank
        self.start=Start
        self.stop=Stop
    
    def __str__(self) :
       return ("Feature number "+str(self.rank)+" with start="+str(self.start)+" and stop="+str(self.stop))

    def set_stop(self,st) :
        self.stop=st

    def get_stop(self) :
        return self.stop

    def set_start(self,st) :
        self.start=st

    def get_start(self) :
        return self.start

    def get_length(self) :
        return self.stop-self.start+1
