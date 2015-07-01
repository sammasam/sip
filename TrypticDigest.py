#!/usr/bin/env python2.7

"""
Sam Bryson
17 June 2015
sammasam@gmail.com

This Class provides code for the creation of a tryptic digest of a protein sequence database

    - this object holds tryptic digest of a protein sequence database in two options
        - set attribute --> unique peps in set
        - dict attribute --> pep:[seq_id1, seq_id2,....seq_idN] dictionary
        
    - methods 
        - load pep_set --> self.pep_set
        - 
"""

import io
import sys

class TrypticDigest():
    
    def __init__(self, name, min_pep = 6, max_pep = 60, cut_after = 'KR', missed_cuts = 2):
        self.name = name
        self.min_pep = min_pep
        self.max_pep = max_pep
        self.cut_after = cut_after
        self.missed_cuts = missed_cuts
        
#------------------------------------------------------------------#        
#                   Initialization Methods                         #
#------------------------------------------------------------------#

    def LoadPepSet (self, inputFile):
        self.pep_list = []
        print("loading peps from: " + inputFile)
        fIN = io.open(inputFile)
        seq = ""
        #seq_id = ""
        for line in fIN:
            line = line.strip()
            if line[0] == ">":
                if seq:
                    self.pep_list += self.GetTrypticPeps(seq)
                    seq = ""
                #seq_id = line[1:]
            else:
                seq += line
        if seq:
            self.pep_list += self.GetTrypticPeps(seq)
            seq = ""
        fIN.close()
        self.pep_set = set(self.pep_list)
        self.pep_list = []

    def GetTrypticPeps (self, seq):
        for aa in self.cut_after:
            cut_seq = seq.replace(aa, aa + "|")
            seq = cut_seq
        trypList = seq.split("|")
        pep = ""
        pepList = []
        for x in range(len(trypList)):
            pep = pep + trypList[x]
            for y in range(x+1,x+1+self.missed_cuts):
                if y < len(trypList):
                    pep = pep + trypList[y]
                    if self.min_pep <= len(pep) <= self.max_pep:
                        pepList.append(pep)
            pep = ""
        return(pepList)
