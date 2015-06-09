#!/usr/bin/env python2.7

"""
Sam Bryson
8 June 2015
sammasam@gmail.com

This Class provides code for the creation of a Sipros proteomics experiment object

    - this object loads multiple SiprosSample objects
        - for data analysis using multiple samples
        - average accros replicates, treatments, time-points, etc
        - perform analysis and create figures
        
    - methods allow for the summary of proteomics data
        - at Protein ID level
        - at hierarchical taxonomic and functional levels
    
    - summary at different level
        - Protein ID counts
        - Balanced spectral counts
        - labeling distribution
        
"""

import io
import sys
from sip.SiprosSample import SiprosSample
import numpy as np
from skbio.diversity.beta import pw_distances as pwd
from skbio.stats.distance import permanova

class SiprosExperiment():
    
    def __init__(self, name):
        self.name = name
        self.samples = []
        self.samples_count = 0
        self.samples_names = []
        self.samples_treatments = []
        self.samples_timepoints = []
        self.samples_locations = []
        self.samples_bsc = []
        self.samples_norm_factors = []
        self.nog_nbsc = {}
        self.protiens = 0
        #self.unique_protCount = 0
        #self.specCount = 0
        #self.unique_specCount = 0
        #self.peptCount = 0
        #self.unique_peptCount = 0

#------------------------------------------------------------------#
#                    Data Input/Calculation Methods                #
#------------------------------------------------------------------#

    def add_sample (self, *argv):
        for sample in argv:
            self.samples.append(sample)
            self.samples_names.append(sample.name)
            self.samples_treatments.append(sample.treatment) 
            self.samples_timepoints.append(sample.timepoint)
            self.samples_locations.append(sample.location)
            self.samples_bsc.append(sample.specCount)
        self.samples_count = len(self.samples)
        self.calculate_normalization_factors()
        print("<>Sample count: "+str(self.samples_count))
        print("Name\tLocation\tTreatment\tTimePoint\tNormFactor")
        for i in range(len(self.samples)):
            n = self.samples_names[i]
            l = self.samples_locations[i]
            t = self.samples_treatments[i]
            m = str(self.samples_timepoints[i])
            f = str(self.samples_norm_factors[i])
            print(n+"\t"+l+"\t"+t+"\t"+m+"\t"+f)
        
        
    def calculate_normalization_factors (self): # calculate normalization factors avg BSC / sample BSC
        #print("<> Calculating normalization factors based on mean spectral counts")
        self.samples_norm_factors = []
        bsc_total = float(sum(self.samples_bsc))
        bsc_avg = bsc_total/self.samples_count
        for bsc in self.samples_bsc:
            self.samples_norm_factors.append(bsc_avg/bsc)

#------------------------------------------------------------------#        
#                       Report Methods                             #
#------------------------------------------------------------------#

    def nog_pbsc_test (self, min_pro_count):
        self.countArray = []
        self.domainList = []
        self.phylumList = []
        self.classList = []
        self.orderList = []
        self.familyList = []
        self.genusList = []
        self.speciesList = []
        self.nogDict = {}
        self.taxaDict = {}
        ranksList = ['domain','phylum','class','order','family','genus','species']
        catsList = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'] #,'na','u']
        for s in range(len(self.samples)):
            sample = self.samples[s]
            name = self.samples_names[s]
            for taxa in sample.summary['species']:
                [D,P,C,O,F,G,S] = taxa.split('|')
                pro_count = 0.0
                spec_count = 0.0
                if S not in ['u', 'na']:
                    nogCounts = [0.0]*len(catsList)
                    for c in range(len(catsList)):
                        nog = catsList[c]
                        if nog in sample.summary['species'][taxa]:
                            [prc, usc, enr] = sample.summary['species'][taxa][nog]
                        else:
                            [prc, usc, enr] = [0,0.0,[0.0]*101]
                        pro_count += prc
                        spec_count += sum(enr)
                        nogCounts[c] = spec_count
                    if pro_count >= min_pro_count:
                        prop_nogCounts = [x/spec_count for x in nogCounts]
                        old_catList = self.nogDict.get(taxa,[[0.0]*len(catsList)]*self.samples_count)
                        old_catList[s] = prop_nogCounts
                        self.nogDict[taxa] = old_catList
                        self.taxaDict[taxa] = [D,P,C,O,F,G,S]  
        for taxa in self.nogDict:
            nog_props = self.nogDict[taxa]
            sum_nog_props = [sum(x) for x in zip(*nog_props)]
            avg_nog_props = [x/self.samples_count for x in sum_nog_props]
            [D,P,C,O,F,G,S] = self.taxaDict[taxa]
            self.countArray.append(avg_nog_props)
            self.domainList.append(D)
            self.phylumList.append(P)
            self.classList.append(C)
            self.orderList.append(O)
            self.familyList.append(F)
            self.genusList.append(G)
            self.speciesList.append(S)
        self.data = np.asarray(self.countArray)
        self.data_eudm = pwd(self.data, self.speciesList, "euclidean")
        print("<> Use this data for performing permanova on NOG category distributions" +
            "\n\t command: permanova(data_eudm,groupList,permutations = 999)" +
            "\n\nAvailable <group>Lists:\n\tgenus, family, order, class, phylum" +
            "\n\n*** from skbio.stats.distance import permanova" +
            "\n*** from skbio.stats.distance import anosim")




#------------------------------------------------------------------#    
#              Functions used in multiple Methods                  #
#------------------------------------------------------------------#

    def update_summary_entry(self):
        self.updated_summary_entry = [0,0.0,[0.0]*101]
        self.updated_summary_entry[0] = self.old_summary_entry[0] + self.new_summary_entry[0]
        self.updated_summary_entry[1] = self.old_summary_entry[1] + self.new_summary_entry[1]
        self.updated_summary_entry[2] = [self.old_summary_entry[2][i] + self.new_summary_entry[2][i] for i in range(0,101)]
        
    def calculate_enrichment_bins (self, enr):
        self.binsList = [0.0]*12
        for x in range(len(enr)):
            b = enr[x]
            if 0<=x<=1:
                a = 0
            elif 2<=x<=10:
                a = 1
            elif 11<=x<=19:
                a = 2
            elif 20<=x<=28:
                a = 3
            elif 29<=x<=37:
                a = 4
            elif 38<=x<=46:
                a = 5
            elif 47<=x<=55:
                a = 6
            elif 56<=x<=64:
                a = 7
            elif 65<=x<=73:
                a = 8
            elif 74<=x<=82:
                a = 9
            elif 83<=x<=91:
                a = 10
            elif 92<=x<=100:
                a = 11
            self.binsList[a] = self.binsList[a] + b
    
    def calculate_enrv (self,enr,lsc):
        enrVals = []
        for x in range(len(enr)):
            y = enr[x]
            enrVals.append(y*x/100)
        enrValSum = sum(enrVals[2:])
        enrSum = sum(enr[2:])
        if enrSum > 0.0:
            self.enrv = enrValSum/enrSum
        else:  
            self.enrv = 0.0




#------------------------------------------------------------------#
#                      IN PROGRESS                                 #
#------------------------------------------------------------------#

"""
#------------------------------------------------------------------#

def taxaXnog_summary (self, rank): # uses self.summary dictionary rank : {taxa: {nog cat: [pro count, usc, [enrList]]}}
        outF = self.name+".pro.taxaXnog_summary."+rank+".txt"
        print("<> Writing summary for taxa and bacNOG categories to file: " + outF)
        print("<> Summary for taxonomic rank: " + rank)
        outFile = io.open(outF, 'a')
        header = unicode("Sample\tRank\tTaxa\tNOG_category\tProIDs\tUSC\tBSC\tLSC\tPBSC\tPLSC\tENRV\tENRF\t0--1\t2--10\t11--19\t" +
                         "20--28\t29--37\t38--46\t47--55\t56--64\t65--73\t74--82\t83--91\t92--100\n")
        outFile.write(header)
        ranksList = ['domain','phylum','class','order','family','genu','species']
        catsList = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','na','u']
        #self.summary[currentRank][taxaName][NOGcategory] = [proCount, unique spec count, [enr 0%,...,100%]]
        for taxa in self.summary[rank]:
            for nog in catsList:
                if nog in self.summary[rank][taxa]:
                    [prc, usc, enr] = self.summary[rank][taxa][nog]
                else:
                    [prc, usc, enr] = [0,0.0,[0.0]*101]
                self.calculate_enrichment_bins(enr) # makes self.binsList
                bsc = sum(enr)
                lsc = sum(enr[2:])
                if bsc > 0.0:
                    pbsc = bsc/self.specCount
                    plsc = lsc/bsc
                    pbins = [(x/bsc) for x in self.binsList]
                    self.calculate_enrv(enr,lsc) # makes self.enrv
                    enrf = self.enrv * plsc
                else:
                    pbsc = 0.0
                    plsc = 0.0
                    pbins = [0.0 for x in self.binsList]
                    self.enrv = 0.0
                    enrf = 0.0
                #print(str(bsc)+"\t"+str(lsc)+"\t"+str(self.enrv)+"\t"+str(enrf))
                descripString = unicode(self.name+"\t"+rank+"\t"+taxa+"\t"+nog+"\t")
                statString = unicode(str(prc)+"\t"+str(usc)+"\t"+str(bsc)+"\t"+str(lsc)+"\t"+str(pbsc)+"\t"+str(plsc)+"\t"+str(self.enrv)+"\t"+str(enrf)+"\t")
                pbinstring = [str(x) for x in pbins]
                binString = unicode("\t".join(pbinstring)+"\n")
                outFile.write(descripString+statString+binString)
        outFile.close()

"""