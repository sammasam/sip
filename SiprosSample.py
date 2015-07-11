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
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import gridspec
import numpy as np
import scipy.stats as stats #t_stat, p_val = stats.ttest_ind(sample1, sample2, equal_var=False)
                            #chisq,  p_val = scipy.stats.chisquare(f_obs, f_exp)


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

    def taxonomy_counts (self, rank, min_count=0, plot="no", data="spectra_counts", print_txt = "", save_fig = ""):
        #print("Counting taxa in proteomics sample\n - Rank: "+rank)
        rank_list = ['domain','phylum','class','order','family','genus','species']
        rank_index = {'domain':0,'phylum':1,'class':2,'order':3,'family':4,'genus':5,'species':6}
        # get counts for each nog category
        self.taxa_count_dict = {}
        self.taxa_count_dict[self.name] = {}
        for protein_string in self.proteins:
            protein_taxonomy = self.proteins[protein_string]['taxonomy']    # [dom,phy,cla,odr,fam,gen,spe]
            name = protein_taxonomy[rank_index[rank]]
            if data == 'spectra_counts':
                self.taxa_count_dict[self.name][name] = self.taxa_count_dict[self.name].get(name, 0.0) + sum(self.proteins[protein_string]['enrichments'])
            if data == 'protein_counts':
                self.taxa_count_dict[self.name][name] = self.taxa_count_dict[self.name].get(name, 0.0) + 1
        # filter taxa list by min count for plotting
        taxa_list = []
        for taxa in self.taxa_count_dict[self.name]:
            count = self.taxa_count_dict[self.name][taxa]
            if count >= min_count:
                taxa_list.append(taxa)
            taxa_set = set(taxa_list)
            taxa_list = list(taxa_set)
            taxa_list.sort()
        # put counts in an array for plotting
        self.taxa_count_list = []   
        for taxa in taxa_list:
            data_point = [self.name, taxa, self.taxa_count_dict[self.name][taxa]]
            self.taxa_count_list.append(data_point)
        if plot == "yes":
            self.make_bar_plot (self.taxa_count_list, taxa_list, start_index = 0, error_list = "", fig_name = save_fig, y_max = "", x_label = "Taxa", y_label = data)

#------------------------------------------------------------------#

    def nog_category_counts (self, plot="no", data="spectra_counts",  print_txt = "", save_fig = ""):
        #print("bactNOG category counts in proteomics sample\n\tCategory\tCount")
        cat_indexes = {'A':0,'B':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,
                       'J':9,'K':10,'L':11,'M':12,'N':13,'O':14,'P':15,'Q':16,
                       'R':17,'S':18,'T':19,'U':20,'V':21,'W':22,'X':23,'Y':24,'Z':25,'na':26,'u':27}
        cat_list = ['A','B','C','D','E','F','G','H','I','J','K','L','M',
                    'N','O','P','Q','R','S','T','U','V','W','X','Y','Z','na','u']
        # get counts for each nog category
        self.nog_category_count_dict = {}
        self.nog_category_count_dict[self.name] = {}
        for cat in cat_list:
            self.nog_category_count_dict[self.name][cat] = 0.0
        for protein_string in self.proteins:
            [nog,cat,des] = self.proteins[protein_string]['bactNOG']
            parsed_cat = self.parse_cat(cat)
            for c in parsed_cat:
                if data == "spectra_counts":
                    self.nog_category_count_dict[self.name][c] += sum(self.proteins[protein_string]['enrichments'])/len(parsed_cat)
                if data == "protein_counts":
                    self.nog_category_count_dict[self.name][c] += 1/len(parsed_cat)
        # put counts in an array for plotting
        self.nog_category_count_list = []   
        for cat in cat_list:
            data_point = [self.name, cat, self.nog_category_count_dict[self.name][cat]]
            self.nog_category_count_list.append(data_point)
        if plot == "yes":
            self.make_bar_plot (self.nog_category_count_list, cat_list, start_index = 0, error_list = "", fig_name = save_fig, y_max = "", x_label = "bactNOG Categories", y_label = data)

#------------------------------------------------------------------#

    def enrichment_counts (self, plot = "no", max_y = "", print_txt = "", save_fig = ""):
        enr_list = ['0--1','2--10','11--19','20--28','29--37','38--46','47--55','56--64','65--73','74--82','83--91','92--100']
        # get spectral counts for each enrichment bin
        self.enrichment_count_dict = {}
        self.enrichment_count_dict[self.name] = {}
        for enr in enr_list:
            self.enrichment_count_dict[self.name][enr] = 0.0
        for protein_string in self.proteins:
            enrichment_bins = self.calculate_enrichment_bins(self.proteins[protein_string]['enrichments'])
            for i in range(len(enrichment_bins)):
                self.enrichment_count_dict[self.name][enr_list[i]] += enrichment_bins[i]
        # put counts in an array for plotting
        self.enrichment_count_list = []   
        for enr in enr_list[1:]:
            data_point = [self.name, enr, self.enrichment_count_dict[self.name][enr]/self.spectra_count]
            self.enrichment_count_list.append(data_point)
        if plot == "yes":
            self.make_bar_plot (self.enrichment_count_list, enr_list, start_index = 1, error_list = "", fig_name = save_fig, y_max = max_y , x_label = "%Enrichment Bins", y_label = "Spectral Counts")
        
#------------------------------------------------------------------#

    def taxa_x_nog_category_counts (self, rank, min_count=0, plot="no", data="spectra_counts", max_y = "", print_txt = "", save_fig = ""):
        
        # add in options for taxa and cat lists ie. plot_taxa = "", and plot_cats = ""
        # do this for all the report methods
        
        rank_list = ['domain','phylum','class','order','family','genus','species']
        rank_index = {'domain':0,'phylum':1,'class':2,'order':3,'family':4,'genus':5,'species':6}
        cat_indexes = {'A':0,'B':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,
                       'J':9,'K':10,'L':11,'M':12,'N':13,'O':14,'P':15,'Q':16,
                       'R':17,'S':18,'T':19,'U':20,'V':21,'W':22,'X':23,'Y':24,'Z':25,'na':26,'u':27}
        cat_list = ['A','B','C','D','E','F','G','H','I','J','K','L','M',
                    'N','O','P','Q','R','S','T','U','V','W','X','Y','Z','na','u']
        # get counts for each nog category for each taxa
        self.taxa_x_nog_category_count_dict = {}
        self.taxa_x_nog_category_count_dict[self.name] = {}
        for protein_string in self.proteins:
            protein_taxonomy = self.proteins[protein_string]['taxonomy']    # [dom,phy,cla,odr,fam,gen,spe]
            taxa_name = protein_taxonomy[rank_index[rank]]
            if taxa_name not in self.taxa_x_nog_category_count_dict[self.name]:
                self.taxa_x_nog_category_count_dict[self.name][taxa_name] = [0.0]*len(cat_list)
            [nog,cat,des] = self.proteins[protein_string]['bactNOG']
            parsed_cat = self.parse_cat(cat)
            for c in parsed_cat:                    ####IF cat in cat list
                if data == 'spectra_counts':
                    self.taxa_x_nog_category_count_dict[self.name][taxa_name][cat_indexes[c]] += sum(self.proteins[protein_string]['enrichments'])/len(parsed_cat)
                if data == 'protein_counts':
                    self.taxa_x_nog_category_count_dict[self.name][taxa_name][cat_indexes[c]] += 1/len(parsed_cat)
        # filter taxa list by min count for plotting
        self.taxa_x_nog_category_count_plot_dict = {}
        taxa_list = []
        for taxa_name in self.taxa_x_nog_category_count_dict[self.name]:
            nog_counts_list = self.taxa_x_nog_category_count_dict[self.name][taxa_name]
            nog_count_sum = sum(nog_counts_list)
            if nog_count_sum >= min_count:
                taxa_list.append(taxa_name)
                self.taxa_x_nog_category_count_plot_dict[taxa_name] = [x/nog_count_sum for x in nog_counts_list]
            taxa_set = set(taxa_list)
            taxa_list = list(taxa_set)
            taxa_list.sort()                  ####IF plot_taxa: taxa_list = <>
        self.make_multi_bar_plot (self.taxa_x_nog_category_count_plot_dict, taxa_list, cat_list, start_index = 0,
                                  error_list = "", fig_name = "NOG counts by Taxa", y_max = "")
           
            
                
                
        
            


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
        #parsed_cat = []
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
            
    def make_bar_plot (self, data_array, labels_list, start_index = 0, error_list = "", fig_name = "", y_max = "", x_label = "", y_label = ""):
        dpoints = np.asarray(data_array)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        group_names = np.unique(dpoints[:,0])
        division_names = np.unique(dpoints[:,1])
        # the space between each set of bars
        space = 0.5
        n = len(group_names)
        width = (1 - space) / (len(group_names))
        # Create a set of bars at each position
        for i,cond in enumerate(group_names):
            print "group:", cond
            indeces = range(1, len(division_names)+1)
            vals = dpoints[dpoints[:,0] == cond][:,2].astype(np.float)
            #errs = dpoints[dpoints[:,0] == cond][:,3].astype(np.float) use 4 th spot in data list for erreor
            pos = [j - (1 - space) / 2. + i * width for j in indeces]
            ax.bar(pos, vals, width=width, label=cond, color=cm.Accent(float(i) / n)) # yerr = errs
        # Set the max value for the y axis
        y_vals = np.unique(dpoints[:,2])
        y_vals = [float(x) for x in y_vals]
        if y_max:
            ax.set_ylim(0, y_max)
        else:
            ax.set_ylim(0,max(y_vals)+0.01)        
        # Set the x-axis tick labels to be equal to the categories
        ax.set_xticks(indeces)
        ax.set_xticklabels(labels_list[start_index:])
        plt.setp(plt.xticks()[1], rotation=90)
        # Add the axis labels
        ax.set_ylabel(y_label)
        ax.set_title(x_label)
        # Add a legend
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='upper right')
        if fig_name:
            plt.savefig(fig_name)
        plt.show()
            
    def make_multi_bar_plot (self, data_dict, title_list, labels_list, start_index = 0, error_list = "", fig_name = "", y_max = ""):
        labels_list = labels_list[start_index:]
        ind = np.arange(len(labels_list))
        width = 0.9
        f, ax = plt.subplots(len(title_list), 1, sharey=True, sharex=True)
        f.subplots_adjust(bottom=0) #make room for the legend
        f.subplots_adjust(hspace=0)
        if y_max:
            max_y = y_max
        else:
            max_y = 1
        #plt.yticks(np.arange(0,max_y ,2))
        plt.xticks(ind,labels_list[start_index:])
        plt.suptitle(fig_name)
        for i in range(len(title_list)):
            vals = data_dict[title_list[i]]
            vals = vals[start_index:]
            ax[i].bar(ind, vals, width=0.9, color=cm.Accent(float(i)/len(vals)), align = 'center')
            ax[i].set_ylim(0, 1)
            ax[i].set_ylabel(title_list[i], verticalalignment='center', rotation='horizontal', horizontalalignment='center',size=6, x=1.0,y=0.5)
        plt.figure(figsize=(8, 10), dpi=300)
        if fig_name:
            plt.savefig(fig_name)
        plt.show()



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
