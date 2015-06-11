#!/usr/bin/env python2.7

"""
Sam Bryson
10 June 2015
sammasam@gmail.com

This Class provides code for the creation of a Sipros database object

    - this object holds sequences and annotations for a proteomics database
        - fasta sequences (amino acid, nucleotide)
        - taxonomic annotation in D|P|C|O|F|G|S style
        - perform analysis and create figures
        
    - methods allow for the summary of database
        - at sequence characteristics level (length, cds completeness, GC%, etc.)
        - at hierarchical taxonomic and functional levels
        
"""

import io
import sys
import matplotlib.pyplot as plt
import numpy as np

class SiprosDatabase():
    
    def __init__(self, name):
        self.name = name
        self.seq_count = 0
        self.seq_lengths = []
        self.faa = {}
        self.fnt = {}
        self.taxonomy = {}
        self.bactNOG = {}
        self.KO = {}
        self.go = {}
        self.SEED = {}
        
        
#------------------------------------------------------------------#        
#                   Initialization Methods                         #
#------------------------------------------------------------------#

    def load_faa (self, fasta_file):
        print("Loading fasta file: " + fasta_file)
        inF = io.open(fasta_file)
        seq = ""
        seq_id = ""
        for line in inF:
            line = line.strip()
            if line[0] == ">":
                if seq:
                    self.seq_count += 1
                    seq_length = len(seq)
                    self.seq_lengths.append(seq_length)
                    self.faa[seq_id] = seq
                    self.taxonomy[seq_id] = ('na','na','na','na','na','na','na')
                    self.bactNOG[seq_id] = ('na','na','na')
                    seq = ""
                seq_id = line[1:]
            else:
                seq += line
        if seq:
            self.seq_count += 1
            seq_length = len(seq)
            self.seq_lengths.append(seq_length)
            self.faa[seq_id] = seq
            self.taxonomy[seq_id] = ('na','na','na','na','na','na','na')
            self.bactNOG[seq_id] = ('na','na','na')
            seq = ""
            seq_id = ""
        inF.close()
        print(" - Total sequences: " + str(self.seq_count))
        
    def load_taxonomy (self, taxonomy_file):
        print("Loading taxonomy annotation file: " + taxonomy_file)
        annotated_count = 0
        inF = io.open(taxonomy_file)                                                            
        for line in inF:
            listLine = line.strip().split("\t")
            if listLine[0] != "ProteinID":
                seq_id = listLine[0]
                dom = listLine[1]
                phy = listLine[2]
                cla = listLine[3]
                odr = listLine[4]
                fam = listLine[5]
                gen = listLine[6]
                spe = listLine[7]
                self.taxonomy[seq_id] = (dom,phy,cla,odr,fam,gen,spe)
                annotated_count += 1
        inF.close()
        print(" - Annotated: "+str(annotated_count))
        print(" - Unannotated: "+str(len(self.taxonomy) - annotated_count))
        
    def load_bactNOG (self, bactNOG_file):
        print("Loading bactNOG annotation file: " + bactNOG_file)
        annotated_count = 0
        inF = io.open(bactNOG_file)                                                              
        for line in inF:
            listLine = line.strip().split("\t")
            if listLine[0] != "ProteinID":
                seq_id = listLine[0]
                nog = listLine[1]
                cat = listLine[2]
                des = listLine[3]
                self.bactNOG[seq_id] = (nog,cat,des)
                annotated_count += 1
        inF.close()
        print(" - Annotated: "+str(annotated_count))
        print(" - Unannotated: "+str(len(self.bactNOG) - annotated_count))
                
#------------------------------------------------------------------#        
#                     Data Access Methods                          #
#------------------------------------------------------------------#

    def get_taxonomy (self, seq_id):
        (dom,phy,cla,odr,fam,gen,spe) = self.taxonomy[seq_id]
        return(dom,phy,cla,odr,fam,gen,spe)

    def get_bactNOG (self, seq_id):
        (nog,cat,des) = self.bactNOG[seq_id]
        return(nog,cat,des)

#------------------------------------------------------------------#        
#                       Report Methods                             #
#------------------------------------------------------------------#

    def nog_category_counts (self, plot="no"):
        print("Counting NOG categories in database\n\tCategory\tCount")
        cat_indexes = {'A':0,'B':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,
                       'J':9,'K':10,'L':11,'M':12,'N':13,'O':14,'P':15,'Q':16,
                       'R':17,'S':18,'T':19,'U':20,'V':21,'W':22,'X':23,'Y':24,'Z':25,'na':26}
        cat_list = ['A','B','C','D','E','F','G','H','I','J','K','L','M',
                    'N','O','P','Q','R','S','T','U','V','W','X','Y','Z','na']
        self.cat_counts = [0.0]*len(cat_list)
        for seq_id in self.bactNOG:
            nog,cat,des = self.bactNOG[seq_id]
            parsed_cat = self.parse_cat(cat)
            for c in parsed_cat:
                self.cat_counts[cat_indexes[c]] += 1/len(parsed_cat)
        for n in range(len(cat_list)):
            category = cat_list[n]
            count = self.cat_counts[n]
            print("\t"+category+"\t"+str(count))
        if plot == "yes":
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ## the data
            N = len(cat_list)                 # number of bactNOG categories
            ## necessary variables
            ind = np.arange(N)                # the x locations for the groups
            width = 0.35                      # the width of the bars
            ## the bars
            rects = ax.bar(ind, self.cat_counts, width,color='black')
            # axes and labels
            ax.set_xlim(-width,len(ind)+width)
            ax.set_ylim(0,int(max(self.cat_counts)+1))
            ax.set_ylabel('Category Counts')
            ax.set_title('bactNOG Categories')
            xTickMarks = cat_list
            ax.set_xticks(ind+width)
            xtickNames = ax.set_xticklabels(xTickMarks)
            plt.setp(xtickNames, rotation=0, fontsize=12)
            plt.show()       
            
    def taxonomy_counts (self, rank, min_count=20, plot="no"):
        print("Counting taxa in database\n - Rank: "+rank)
        rank_list = ['domain','phylum','class','order','family','genus','species']
        rank_index = {'domain':0,'phylum':1,'class':2,'order':3,'family':4,'genus':5,'species':6}
        self.rank_counts = {}
        self.rank_counts[rank] = {}
        taxa_list = []
        taxa_counts = []
        for seq_id in self.taxonomy:
            name = self.taxonomy[seq_id][rank_index[rank]]
            self.rank_counts[rank][name] = self.rank_counts[rank].get(name,0) + 1
        num_taxa = len(self.rank_counts[rank])
        for n in self.rank_counts[rank]:
            count = self.rank_counts[rank][n]
            #print("\t"+n+"\t"+str(count))
            if count >= min_count:
                taxa_list.append(n)
        taxa_list.sort()
        for t in taxa_list:    
            taxa_counts.append(self.rank_counts[rank][t])
        print(" - Total taxa: "+str(num_taxa))
        if plot == "yes":
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ## the data
            N = len(taxa_list)                # number of taxa at selected rank
            ## necessary variables
            ind = np.arange(N)                # the x locations for the groups
            width = 0.35                      # the width of the bars
            ## the bars
            rects = ax.bar(ind, taxa_counts, width,color='black')
            # axes and labels
            ax.set_xlim(-width,len(ind)+width)
            ax.set_ylim(0,int(max(taxa_counts)+1))
            ax.set_ylabel('Sequence Counts')
            ax.set_title('Taxonomic Distribution of CDS')
            xTickMarks = taxa_list
            ax.set_xticks(ind+width)
            xtickNames = ax.set_xticklabels(xTickMarks)
            plt.setp(xtickNames, rotation=270, fontsize=10)
            plt.show()
            
            
            
#------------------------------------------------------------------#    
#              Functions used in multiple Methods                  #
#------------------------------------------------------------------#

    def parse_cat (self,cat):
        parsed_cat = []
        if cat == 'na':
            parsed_cat = ['na']
        else:
            for c in cat:
                parsed_cat.append(c)
        return(parsed_cat)
                
#------------------------------------------------------------------#
#                      IN PROGRESS                                 #
#------------------------------------------------------------------#

"""
#------------------------------------------------------------------#

"""

