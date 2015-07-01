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
import matplotlib.pyplot as plt
import numpy as np

class SiprosSample():
    
    def __init__(self, sample_name, sample_treatment, sample_timepoint, sample_location, sample_database, pro_txt, pep_txt):
        self.name = sample_name
        self.treatment = sample_treatment
        self.timepoint = sample_timepoint
        self.location = sample_location
        self.sample_database = sample_database
        self.pro_txt = pro_txt
        self.pep_txt = pep_txt
        self.proteins = {}  # protein_ID => taxonomy => (dom,phy,cla,odr,fam,gen,spe)
                            #               bactNOG_annotation => (nog, cat, des)
                            #               peptides => ['ASDFGHJK',......]
                            #               enrichments => [0.0]*101
                            #               spectra => [unique_count,shared_count]
        self.protein_id_lookup_dict = {}
        self.protein_count = 0
        self.unique_protein_count = 0
        self.peptide_count = 0
        self.spectra_count = 0
        self.unique_spectra_count = 0
        self.normalization_factor = 1
        self.print_metadata()
        self.load_proteins()
        self.load_peptides()
    
#------------------------------------------------------------------#
#                    Object Initialization Methods                 #
#------------------------------------------------------------------#

    def print_metadata (self):
        print("SiprosSample Object Created:" +
              "\n\tName : "+self.name+"\n\tLocation : "+self.location+
              "\n\tTreatment : "+self.treatment+"\n\tTimepoint : "+str(self.timepoint))
    
    def load_proteins (self):
        print("Loading pro.txt file: " + self.pro_txt)
        inFile = io.open(self.pro_txt)
        for line in inFile:
            if line[0] != '#':
                line_list = line.strip().split("\t")
                if line_list[8] == 'T':
                    self.protein_count += 1
                    protein_string = line_list[0]
                    protein_list = self.pro_to_list(protein_string)
                    taxonomy_annotation = self.annotate_taxonomy(protein_list)
                    bactNOG_annotation = self.annotate_bactNOG(protein_list)
                    if len(protein_list) == 1:
                        self.unique_protein_count += 1
                    self.proteins[protein_string] = {}
                    self.proteins[protein_string]['taxonomy'] = taxonomy_annotation  # ['D','P','C','O','F','G','S'] 
                    self.proteins[protein_string]['bactNOG'] = bactNOG_annotation   # [bactNOG_number, bactNOG_category, description]
                    self.proteins[protein_string]['peptides'] = []
                    self.proteins[protein_string]['enrichments'] = [0.0]*101
                    self.proteins[protein_string]['peptide_counts'] = [int(line_list[1]),int(line_list[2])] # [unique, total]
                    self.proteins[protein_string]['spectra_counts'] = [int(line_list[3]),int(line_list[4])] # [unique, total]
                    # self.proteins[protein_string]['BSC']  # Balanced Spectra Counts
                    # self.proteins[protein_string]['PBSC'] # Proportion of Balanced Spectra Counts
                    # self.proteins[protein_string]['PESC'] # Proportion of Spectra Counts with Enrichment
                    # self.proteins[protein_string]['AESC'] # Average %Enrichment of Labeled Spectral
                    for prot in protein_list:
                        self.protein_id_lookup_dict[prot] = protein_string
        inFile.close()
        print("\tUnique Protein IDs: " + str(self.unique_protein_count))
        print("\tTotal Protein IDs: " + str(self.protein_count))
        
    def load_peptides (self):
        print("Loading pep.txt file: " + self.pep_txt)
        inFile = io.open(self.pep_txt)
        for line in inFile:
            if line[0] != '#':
                line_list = line.strip().split("\t")
                if line_list[5] == 'T':
                    peptide = line_list[2][1:-1]
                    self.peptide_count += 1
                    protein_string = line_list[3]
                    protein_list = self.pro_to_list(protein_string)
                    pro_list = []
                    for pro in protein_list:    # not every protein with peptide ID is a protein ID passing 1 unique + 1 other peptide
                        if pro in self.protein_id_lookup_dict:
                            pro_list.append(self.protein_id_lookup_dict[pro])
                    pro_set = set(pro_list)
                    pro_count = len(pro_set)
                    enrichment_string = line_list[10]
                    enrichment_list = self.parse_enrichments(enrichment_string)
                    if len(pro_set) == 1:
                            self.unique_spectra_count += len(enrichment_list)
                    if pro_count > 0:
                        val = 1.0/pro_count
                        self.spectra_count += len(enrichment_list)
                        for p in pro_set:
                            for enr in enrichment_list:
                                self.proteins[p]['enrichments'][enr] += val
                            self.proteins[p]['peptides'].append(peptide)
        inFile.close()
        print("\tUnique Spectra: " + str(self.unique_spectra_count))
        print("\tTotal Spectra: " + str(self.spectra_count))
        print("\tTotal Peptide IDs: " + str(self.peptide_count))
        
    def annotate_taxonomy (self,protein_list):
        taxonomy_annotation = ['u','u','u','u','u','u','u']
        dom_list = []
        phy_list = []
        cla_list = []
        odr_list = []
        fam_list = []
        gen_list = []
        spe_list = []
        for protein_id in protein_list:
            (dom,phy,cla,odr,fam,gen,spe) = self.sample_database.get_taxonomy(protein_id)
            dom_list.append(dom)
            phy_list.append(phy)
            cla_list.append(cla)
            odr_list.append(odr)
            fam_list.append(fam)
            gen_list.append(gen)
            spe_list.append(spe);
        taxa = [dom_list,phy_list,cla_list,odr_list,fam_list,gen_list,spe_list]
        for x in range(len(taxa)):
            y = taxa[x];
            z = set(y);
            if len(z) == 1:
                taxonomy_annotation[x] = y[0]
            else:
                taxonomy_annotation[x] = 'u'
        return(taxonomy_annotation)

    def annotate_bactNOG (self, protein_list):
        bactNOG_annotation = ['u','u','u']
        bactNOG_numbers = []
        bactNOG_categories = []
        bactNOG_descriptions = []
        for protein_id in protein_list:
            (nog,cat,des) = self.sample_database.get_bactNOG(protein_id)
            bactNOG_numbers.append(nog)
            bactNOG_categories.append(cat)
            bactNOG_descriptions.append(des)
        bactNOG_annotation_lists = [bactNOG_numbers,bactNOG_categories,bactNOG_descriptions]
        for x in range(len(bactNOG_annotation_lists)):
            y = bactNOG_annotation_lists[x]
            z = set(y)
            if len(z) == 1:
                bactNOG_annotation[x] = y[0]
            else:
                bactNOG_annotation[x] = 'u'
        return(bactNOG_annotation)

#------------------------------------------------------------------#        
#                       Report Methods                             #
#------------------------------------------------------------------#

    def taxonomy_counts (self, rank, min_count=20, plot="no", data="spectra_counts"):
        print("Counting taxa in proteomics sample\n - Rank: "+rank)
        rank_list = ['domain','phylum','class','order','family','genus','species']
        rank_index = {'domain':0,'phylum':1,'class':2,'order':3,'family':4,'genus':5,'species':6}
        self.rank_counts = {}
        self.rank_counts[rank] = {}
        taxa_list = []
        taxa_counts = []
        for protein_string in self.proteins:
            protein_taxonomy = self.proteins[protein_string]['taxonomy']    # [dom,phy,cla,odr,fam,gen,spe]
            name = protein_taxonomy[rank_index[rank]]
            if data == 'spectra_counts':
                self.rank_counts[rank][name] = self.rank_counts[rank].get(name,0.0) + sum(self.proteins[protein_string]['enrichments'])
            if data == 'protein_counts':
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
            if data == 'spectra_counts':
                ax.set_ylabel('Balanced Spectra Counts')
            if data == 'protein_counts':
                ax.set_ylabel('Protein Id Counts')
            ax.set_title('Taxonomic Distribution')
            xTickMarks = taxa_list
            ax.set_xticks(ind+width)
            xtickNames = ax.set_xticklabels(xTickMarks)
            plt.setp(xtickNames, rotation=270, fontsize=10)
            plt.show()

    def nog_category_counts (self, plot="no", data="spectra_counts"):
        print("bactNOG category counts in proteomics sample\n\tCategory\tCount")
        cat_indexes = {'A':0,'B':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,
                       'J':9,'K':10,'L':11,'M':12,'N':13,'O':14,'P':15,'Q':16,
                       'R':17,'S':18,'T':19,'U':20,'V':21,'W':22,'X':23,'Y':24,'Z':25,'na':26,'u':27}
        cat_list = ['A','B','C','D','E','F','G','H','I','J','K','L','M',
                    'N','O','P','Q','R','S','T','U','V','W','X','Y','Z','na','u']
        self.protein_cat_counts = [0.0]*len(cat_list)
        self.spectra_cat_counts = [0.0]*len(cat_list)
        for protein_string in self.proteins:
            [nog,cat,des] = self.proteins[protein_string]['bactNOG']
            parsed_cat = self.parse_cat(cat)
            for c in parsed_cat:
                self.spectra_cat_counts[cat_indexes[c]] += sum(self.proteins[protein_string]['enrichments'])/len(parsed_cat) # add fraction of total spectra
                self.protein_cat_counts[cat_indexes[c]] += 1/len(parsed_cat)
        for n in range(len(cat_list)):
            category = cat_list[n]
            if data == 'spectra_counts':
                count = self.spectra_cat_counts[n]
            if data == 'protein_counts':
                count = self.protein_cat_counts[n]
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
            if data == 'spectra_counts':
                rects = ax.bar(ind, self.spectra_cat_counts, width,color='black')
            if data == 'protein_counts':
                rects = ax.bar(ind, self.protein_cat_counts, width,color='black')
            # axes and labels
            ax.set_xlim(-width,len(ind)+width)
            if data == 'spectra_counts':
                ax.set_ylim(0,int(max(self.spectra_cat_counts)+1))
            if data == 'protein_counts':
                ax.set_ylim(0,int(max(self.protein_cat_counts)+1))
            ax.set_ylabel('Category Counts')
            ax.set_title('bactNOG Categories')
            xTickMarks = cat_list
            ax.set_xticks(ind+width)
            xtickNames = ax.set_xticklabels(xTickMarks)
            plt.setp(xtickNames, rotation=0, fontsize=12)
            plt.show()       

    def enrichment_counts (self, plot="no"):
        enr_list = ['0--1','2--10','11--19','20--28','29--37','38--46','47--55','56--64','65--73','74--82','83--91','92--100']
        self.spectra_enr_counts = [0.0]*len(enr_list)
        for protein_string in self.proteins:
            enrichment_bins = self.calculate_enrichment_bins(self.proteins[protein_string]['enrichments'])
            self.spectra_enr_counts = [ x + y for x,y in zip(self.spectra_enr_counts,enrichment_bins)]
        if plot == "yes":
            print("plotting enrichment histogram")
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ## the data
            N = len(enr_list) - 1             # number of enrichment bins
            ## necessary variables
            ind = np.arange(N)                # the x locations for the groups
            width = 0.35                      # the width of the bars
            ## the bars
            rects = ax.bar(ind, self.spectra_enr_counts[1:], width,color='black')
            # axes and labels
            ax.set_xlim(-width,len(ind)+width)
            ax.set_ylim(0,int(max(self.spectra_enr_counts[1:])+1))
            ax.set_ylabel('Spectra Counts')
            ax.set_title('Enrichment Bins')
            xTickMarks = enr_list[1:]
            ax.set_xticks(ind+width)
            xtickNames = ax.set_xticklabels(xTickMarks)
            plt.setp(xtickNames, rotation=225, fontsize=12)
            plt.show()       
        

#------------------------------------------------------------------#    
#              Functions used in multiple Methods                  #
#------------------------------------------------------------------#
        
    def pro_to_list (self, protein_string): #input pro IDs {idA,idB,...,idN} -> return [idA,...idN]
        protein_list = []
        if protein_string[0] == "{":    
            protein_string = protein_string[1:-1]
        protein_list = protein_string.split(',')
        protein_list = [p for p in protein_list if 'Rev' not in p]
        return(protein_list)

    def parse_enrichments (self, enrichment_string): #input string "{C13_1Pct,C13_1Pct}" -> output int list [1,1]
        enrichment_list = []
        if "Null" in enrichment_string:
            enrichment_list = [1]*enrichment_string.count("Null")
        else:
            enrichment_string = enrichment_string[1:-1]
            enrichment_list = enrichment_string.replace('C13_','').replace('Pct','').split(',')
            enrichment_list = [int(e) for e in enrichment_list]
        return(enrichment_list)
    
    def parse_cat (self,cat):
        parsed_cat = []
        if 'na' in cat:
            parsed_cat = ['na']
        else:
            parsed_cat = [c for c in cat]
        return(parsed_cat)
            
    def update_summary_entry(self):
        self.updated_summary_entry = [0,0.0,[0.0]*101]
        self.updated_summary_entry[0] = self.old_summary_entry[0] + self.new_summary_entry[0]
        self.updated_summary_entry[1] = self.old_summary_entry[1] + self.new_summary_entry[1]
        self.updated_summary_entry[2] = [self.old_summary_entry[2][i] + self.new_summary_entry[2][i] for i in range(0,101)]
        
    def calculate_enrichment_bins (self, enr):
        bins_list = [0.0]*12
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
            bins_list[a] = bins_list[a] + b
        return(bins_list)
    
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
                        
#------------------------------------------------------------------#

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

#------------------------------------------------------------------#

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



#------------------------------------------------------------------#
                        
                        
"""
