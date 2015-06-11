#!/usr/bin/env python2.7

"""
Sam Bryson
3 June 2015
sammasam@gmail.com

This Class provides code for the creation of a Sipros proteomics sample object

    - this object loads informaton about identified proteins
        in a sample
        - taxonomic annotations <> requires annotated database
        - fanctional annotations <> requires annotated database
        - labeling information from SIP experiment
        
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
import pandas as pd
import sip

class SiprosSample():
    
    def __init__(self, sample_name, sample_treatment, sample_timepoint, sample_location, taxonomy_file, bactNOG_file, pro_txt, pep_txt):
        self.name = sample_name
        self.treatment = sample_treatment
        self.timepoint = sample_timepoint
        self.location = sample_location
        self.taxonomy_file = taxonomy_file
        self.bactNOG_file = bactNOG_file
        self.pro_txt = pro_txt
        self.pep_txt = pep_txt
        self.proteins = {} # protein name : (sample_proteinID, (nog, cat, des), (dom,phy,cla,odr,fam,gen,spe))
        self.peptides = {}
        self.spectra = {}
        self.enrichments = {}
        self.taxonomy = {}
        self.bactNOG = {}
        self.summary = {}
        self.protCount = 0
        self.unique_protCount = 0
        self.specCount = 0
        self.unique_specCount = 0
        self.peptCount = 0
        self.unique_peptCount = 0
        self.print_metadata()
        self.load_proteins()
        self.load_peptides()
        self.annotate_taxonomy()
        self.annotate_bactNOG()
        self.build_sample_summary()
    
#------------------------------------------------------------------#
#                    Object Initialization Methods                 #
#------------------------------------------------------------------#

    def print_metadata (self):
        print("SiprosSample Object Created:" +
              "\n\tName : "+self.name+"\n\tLocation : "+self.location+
              "\n\tTreatment : "+self.treatment+"\n\tTimepoint : "+self.timepoint)
    
    def load_proteins (self):
        print("Loading pro.txt file: " + self.pro_txt)
        inFile = io.open(self.pro_txt)
        for line in inFile:
            if line[0] != '#':
                listLine = line.strip().split("\t")
                if listLine[8] == 'T':
                    self.proteinString = listLine[0]
                    
                    # count protein IDs #
                    self.protCount += 1
                    if self.proteinString[0] != "{":
                        self.unique_protCount += 1
                    
                    # initialize all dictionaries with protein IDs #
                    self.pro_to_list()
                    for p in self.proteinList:
                        self.proteins[p] = [self.proteinString,('na','na','na'),('na','na','na','na','na','na','na')]
                    self.unique_spectra[self.proteinString] = int(listLine[3])
                    self.peptides[self.proteinString] = [0,0]
                    self.enrichments[self.proteinString] = [0.0]*101 # bins 0% to 100%
                    self.taxonomy[self.proteinString] = ['na','na','na','na','na','na','na'] # dom|phy|cla|ord|fam|gen|spe
                    self.bactNOG[self.proteinString] = ['na','na','na'] # bactNOG|NOGcategory|description
        inFile.close();
        print("\tUnique Protein IDs: " + str(self.unique_protCount))
        print("\tTotal Protein IDs: " + str(self.protCount))
        
    def load_peptides (self):
        print("Loading pep.txt file: " + self.pep_txt)
        inFile = io.open(self.pep_txt)
        for line in inFile:
            if line[0] != '#':
                listLine = line.strip().split("\t")
                if listLine[5] == 'T':
                    self.peptCount += 1
                    if len(proSet) == 1:
                        self.unique_peptCount += 1
                    proIDs = []
                    self.proteinString = listLine[3]
                    self.pro_to_list()
                    for p in self.proteinList:
                        if p in self.proteins:
                            pro = self.proteins[p][0]
                            proIDs.append(pro)
                    proSet = set(proIDs)
                    if len(proSet) > 0:
                        self.enrichmentString = listLine[10]
                        self.parse_enrichments()
                        val = 1.0/len(proSet)
                        
                        # count peptides #
                        for pp in proSet:
                            if len(proSet) == 1:
                                self.peptides[pp][0] += 1
                            else:
                                self.peptides[pp][1] += 1
                                
                        # count spectra #
                        self.specCount += len(self.enrichmentList)
                        if len(proSet) == 1:
                            self.unique_specCount += len(self.enrichmentList)
                        
                        # update self.enrichments dictionary #
                        for enr in self.enrichmentList:
                            for prot in proSet:
                                self.enrichments[prot][enr] += val
        inFile.close();
        print("\tUnique Spectra: " + str(self.unique_specCount))
        print("\Total Spectra: " + str(self.specCount))
        print("\tUnique Peptide IDs: " + str(self.unique_peptCount))
        print("\tTotal Peptide IDs: " + str(self.peptCount))

    def annotate_bactNOG (self):
        print("Loading bactNOG annotation file: " + self.bactNOG_file)
        inFile = io.open(self.bactNOG_file)                                                              
        for line in inFile:
            listLine = line.strip().split("\t")
            if listLine[0] != "ProteinID":
                pro = listLine[0]
                nog = listLine[1]
                cat = listLine[2]
                des = listLine[3]
                if pro in self.proteins:
                    self.proteins[pro][1] = (nog,cat,des)
        inFile.close()
        print("Annotating protein IDs with bactNOGs")
        for p in self.bactNOG:
            self.proteinString = p
            self.pro_to_list()
            self.get_LCF()
            self.bactNOG[p] = self.lcf
            
    def annotate_taxonomy (self):
        print("Loading taxonomy annotation file: " + self.taxonomy_file)
        inFile = io.open(self.taxonomy_file)                                                              
        for line in inFile:
            listLine = line.strip().split("\t")
            if listLine[0] != "ProteinID":
                pro = listLine[0]
                dom = listLine[1]
                phy = listLine[2]
                cla = listLine[3]
                odr = listLine[4]
                fam = listLine[5]
                gen = listLine[6]
                spe = listLine[7]
                if pro in self.proteins:
                    self.proteins[pro][2] = (dom,phy,cla,odr,fam,gen,spe)
        inFile.close()
        print("Annotating protein IDs with taxonomy")
        for p in self.taxonomy:
            self.proteinString = p
            self.pro_to_list()
            self.get_LCA()
            self.taxonomy[p] = self.lca
            
    def build_sample_summary (self):
        self.summary = {}
        print("Calculating summary of spectra and protein IDs by taxa and NOG category")
        for p in self.spectra:
            self.new_summary_entry = []
            usc = self.spectra[p] # integer of unique spectral counts
            enr = self.enrichments[p] # list of spectral counts by index/bin 0% to 100%
            pep = self.peptides[p]
            [nog,cat,des] = self.bactNOG[p]
            [dom,phy,cla,odr,fam,gen,spe] = self.taxonomy[p]
            
            # maybe add in counts for each functional and taxonomic category #
            
            prc = 1.0 # for counting total protein IDs
            self.cats = cat
            self.parse_categories()
            numCats = len(self.catList)
            enr = [e/numCats for e in enr]
            prc = 1.0/numCats
            usc = usc/numCats
            self.new_summary_entry = [prc,usc,enr]
            D = dom
            P = dom+"|"+phy
            C = dom+"|"+phy+"|"+cla
            O = dom+"|"+phy+"|"+cla+"|"+odr
            F = dom+"|"+phy+"|"+cla+"|"+odr+"|"+fam
            G = dom+"|"+phy+"|"+cla+"|"+odr+"|"+fam+"|"+gen
            S = dom+"|"+phy+"|"+cla+"|"+odr+"|"+fam+"|"+gen+"|"+spe
            phylogeny = [D,P,C,O,F,G,S]
            ranksList = ['domain','phylum','class','order','family','genus','species']
            for i in range(len(ranksList)):
                currentRank = ranksList[i]
                taxaName = phylogeny[i]
                for NOGcategory in self.catList:
                    #print(currentRank+"\t"+taxaName+"\t"+NOGcategory)
                    if currentRank not in self.summary:
                        self.summary[currentRank] = {}
                    if taxaName not in self.summary[currentRank]:
                        self.summary[currentRank][taxaName] = {}
                    if NOGcategory not in self.summary[currentRank][taxaName]:
                        self.summary[currentRank][taxaName][NOGcategory] = self.new_summary_entry
                    else:
                        self.old_summary_entry = self.summary[currentRank][taxaName][NOGcategory]
                        self.update_summary_entry()
                        self.summary[currentRank][taxaName][NOGcategory] = self.updated_summary_entry
                    

#------------------------------------------------------------------#        
#                       Report Methods                             #
#------------------------------------------------------------------#

    def protein_summary (self):
        outF = self.name+".pro.protein_summary.txt"
        print("<> writing protein annotations and enrichments file: " + outF)
        outFile = io.open(outF, 'a')
        header = unicode("ProteinID\tDomain\tPhyla\tClass\tOrder\tFamily\tGenus\tSpecies\t" +
                         "bactNOG\tNOG_category\tNOG_description\t" +
                         "UniqueSpectra(USC)\tBalancedSpectra(BSC)\tBins[0-100]\n")
        outFile.write(header)
        for p in self.spectra:
            usc = self.spectra[p]
            enr = self.enrichments[p]
            bsc = sum(enr)
            enrString ="|".join([str(e) for e in enr])
            [nog,cat,des] = self.bactNOG[p]
            [dom,phy,cla,odr,fam,gen,spe] = self.taxonomy[p]
            outFile.write(unicode(p+"\t"+dom+"\t"+phy+"\t"+cla+"\t"+odr+"\t"+fam+"\t"+gen+"\t"+spe+"\t"+
                                  nog+"\t"+cat+"\t"+des+"\t"+str(usc)+"\t"+str(bsc)+"\t"+enrString+"\n"))
        outFile.close()

    def taxaXnog_summary (self, rank): # uses self.summary dictionary rank : {taxa: {nog cat: [pro count, usc, [enrList]]}}
        outF = self.name+".pro.taxaXnog_summary."+rank+".txt"
        print("Writing summary for taxa and bacNOG categories to file: " + outF)
        print("Summary for taxonomic rank: " + rank)
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
    
    def taxa_summary (self, rank): # uses self.summary dictionary rank : {taxa: {nog cat: [pro count, usc, [enrList]]}}
        outF = self.name+".pro.taxa_summary."+rank+".txt"
        print("Writing summary for taxa to file: " + outF)
        print("Summary for taxonomic rank: " + rank)
        outFile = io.open(outF, 'a')
        header = unicode("Sample\tRank\tTaxa\tProIDs\tUSC\tBSC\tLSC\tPBSC\tPLSC\tENRV\tENRF\t0--1\t2--10\t11--19\t" +
                         "20--28\t29--37\t38--46\t47--55\t56--64\t65--73\t74--82\t83--91\t92--100\n")
        outFile.write(header)
        #ranksList = ['domain','phylum','class','order','family','genu','species']
        #self.summary[currentRank][taxaName][NOGcategory] = [proCount, unique spec count, [enr 0%,...,100%]]
        for taxa in self.summary[rank]:
            prc = 0.0
            usc = 0.0
            enr = [0.0]*101
            for nog in self.summary[rank][taxa]:
                [in_prc, in_usc, in_enr] = self.summary[rank][taxa][nog]
                prc += in_prc
                usc = in_usc
                enr = [enr[i]+in_enr[i] for i in range(len(enr))]
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
            descripString = unicode(self.name+"\t"+rank+"\t"+taxa+"\t")
            statString = unicode(str(prc)+"\t"+str(usc)+"\t"+str(bsc)+"\t"+str(lsc)+"\t"+str(pbsc)+"\t"+str(plsc)+"\t"+str(self.enrv)+"\t"+str(enrf)+"\t")
            pbinstring = [str(x) for x in pbins]
            binString = unicode("\t".join(pbinstring)+"\n")
            outFile.write(descripString+statString+binString)
        outFile.close()
    

#------------------------------------------------------------------#    
#              Functions used in multiple Methods                  #
#------------------------------------------------------------------#
        
    def pro_to_list (self): #input pro IDs {idA,idB,...,idN} -> return [idA,...idN]
        self.proteinList = []
        tempString = self.proteinString
        if tempString[0] == "{":    
            tempString = tempString[1:-1]
        self.proteinList = tempString.split(',')
        self.proteinList = [p for p in self.proteinList if 'Rev' not in p]

    def parse_enrichments (self): #input string "{C13_1Pct,C13_1Pct}" -> output int list [1,1]
        self.enrichmentList = []
        if "Null" in self.enrichmentString:
            self.enrichmentList = [1]*self.enrichmentString.count("Null")
        else:
            self.enrichmentString = self.enrichmentString[1:-1]
            self.enrichmentList = self.enrichmentString.replace('C13_','').replace('Pct','').split(',')
            self.enrichmentList = [int(e) for e in self.enrichmentList]

    def get_LCF (self): # prots <= {pro_1,pro_2,...pro_n}
        self.lcf = ['u','u','u']
        nogS = []
        catS = []
        desS = []
        for p in self.proteinList:
            (nog,cat,des) = self.proteins[p][1]
            nogS.append(nog)
            catS.append(cat)
            desS.append(des)
        func = [nogS,catS,desS]
        for x in range(len(func)):
            y = func[x];
            z = set(y);
            if len(z) == 1:
                self.lcf[x] = y[0];
            else:
                self.lcf[x] = 'u'

    def get_LCA (self): # prots <= {pro_1,pro_2,...pro_n}
        self.lca = ['u','u','u','u','u','u','u']
        domS = []
        phyS = []
        claS = []
        odrS = []
        famS = []
        genS = []
        speS = []
        for p in self.proteinList:
            (dom,phy,cla,odr,fam,gen,spe) = self.proteins[p][2]
            domS.append(dom)
            phyS.append(phy)
            claS.append(cla)
            odrS.append(odr)
            famS.append(fam)
            genS.append(gen)
            speS.append(spe);
        taxa = [domS,phyS,claS,odrS,famS,genS,speS]
        for x in range(len(taxa)):
            y = taxa[x];
            z = set(y);
            if len(z) == 1:
                self.lca[x] = y[0]
            else:
                self.lca[x] = 'u'
        
    def parse_categories(self): # could be "A", "AB", or "na"
        if 'na' in self.cats:
            self.catList = ['na']
        else: 
            self.catList = [c for c in self.cats]
            
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

"""
