#!/usr/bin/env python2.7

"""
Sam Bryson
1 July 2015
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
import sip
import numpy as np
from skbio.diversity.beta import pw_distances as pwd
from skbio.stats.distance import permanova
from matplotlib.mlab import PCA
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d


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
        
        self.groups = []    # used in self.group_samples(self, sample_list, sample_name)
        self.group_names = []
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
            self.samples_bsc.append(sample.spectra_count)
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
            
    def group_samples (self, sample_list, sample_name):
        self.groups.append(sample_list)
        self.group_names.append(sample_name)
            
            
#------------------------------------------------------------------#        
#                       Report Methods                             #
#------------------------------------------------------------------#

    def make_bar_plot (self, data_array, labels_list, start_index = 0, error_bars = "no", fig_name = "", y_max = "", x_label = "", y_label = ""):
        dpoints = np.asarray(data_array)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        group_names = np.unique(dpoints[:,0])
        division_names = np.unique(dpoints[:,1])
        # the space between each set of bars
        space = 0.1
        n = len(group_names)
        width = (1 - space) / (len(group_names))
        # Create a set of bars at each position
        for i,cond in enumerate(group_names):
            print "group:", cond
            indeces = range(1, len(division_names)+1)
            vals = dpoints[dpoints[:,0] == cond][:,2].astype(np.float)
            pos = [j - (1 - space) / 2. + i * width for j in indeces]
            if error_bars == "yes":
                errs = dpoints[dpoints[:,0] == cond][:,3].astype(np.float) #use 4 th spot in data list for error value
                ax.bar(pos, vals, width=width, label=cond, color=cm.Accent(float(i) / n) ,
                       yerr = errs, error_kw = dict(ecolor='black', lw=1, capsize=1, capthick=1))
            else:
                ax.bar(pos, vals, width=width, label=cond, color=cm.Accent(float(i) / n))
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
        
        
#------------------------------------------------------------------#
#------------------------------------------------------------------#

    def taxonomy_counts (self, rank, minimum_count=0, plot="no", data_type="spectra_counts", max_y = "", print_txt = "", save_fig = ""):
        #print("Counting taxa in proteomics samples\n - Rank: "+rank)
        rank_list = ['domain','phylum','class','order','family','genus','species']
        rank_index = {'domain':0,'phylum':1,'class':2,'order':3,'family':4,'genus':5,'species':6}
        # get counts for each taxa from each sample
        self.taxa_counts_dict = {}
        taxa_filter_dict = {}
        for sample in self.samples:
            sample.taxonomy_counts(rank, data=data_type)
            self.taxa_counts_dict[sample.name] = {}
            for taxa in sample.taxa_count_dict[sample.name]:
                count = sample.taxa_count_dict[sample.name][taxa]
                self.taxa_counts_dict[sample.name][taxa] = count
                if taxa in taxa_filter_dict:
                    taxa_filter_dict[taxa].append(count)
                else:
                    taxa_filter_dict[taxa] = [count]
        # filter taxa by minimum count
        taxa_list = []
        for taxa in taxa_filter_dict:
            count_list = taxa_filter_dict[taxa]
            if max(count_list) >= minimum_count:
                taxa_list.append(taxa)
        taxa_set = set(taxa_list)
        taxa_list = list(taxa_set)
        taxa_list.sort()
        # put counts in an array for plotting
        self.taxa_counts_list = []
        for sample in self.samples:
            for taxa in taxa_list:
                data_point = [sample.name, taxa, self.taxa_counts_dict[sample.name].get(taxa,0.0)]
                self.taxa_counts_list.append(data_point)
        if plot == "yes":
            self.make_bar_plot (self.taxa_counts_list, taxa_list, start_index = 0, error_bars = "",
                                fig_name = save_fig, y_max = max_y, x_label = "Taxa", y_label = data_type)

#------------------------------------------------------------------#

    def group_average_taxonomy_counts (self, rank, minimum_count=0, plot="no", data_type="spectra_counts", max_y = "", print_txt = "", save_fig = ""):
        #print("Counting taxa in proteomics samples\n - Rank: "+rank)
        rank_list = ['domain','phylum','class','order','family','genus','species']
        rank_index = {'domain':0,'phylum':1,'class':2,'order':3,'family':4,'genus':5,'species':6}
        # get counts for each taxa from each sample
        self.group_taxa_counts_dict = {}
        taxa_filter_dict = {}
        for g in range(len(self.group_names)):
            group_name = self.group_names[g]
            group_list = self.groups[g]
            self.group_taxa_counts_dict[group_name] = {}
            for s in range(len(group_list)):
                sample = group_list[s]
                sample.taxonomy_counts(rank, data=data_type)
                for taxa in sample.taxa_count_dict[sample.name]:
                    if data_type == "spectra_counts":
                        count = sample.taxa_count_dict[sample.name][taxa]/sample.spectra_count
                    else:
                        count = sample.taxa_count_dict[sample.name][taxa]
                    if taxa in self.group_taxa_counts_dict[group_name]:
                        self.group_taxa_counts_dict[group_name][taxa].append(count)
                    else:
                        self.group_taxa_counts_dict[group_name][taxa] = [count]
                    if taxa in taxa_filter_dict:
                        taxa_filter_dict[taxa].append(count)
                    else:
                        taxa_filter_dict[taxa] = [count]
        # filter taxa by minimum count
        taxa_list = []
        for taxa in taxa_filter_dict:
            count_list = taxa_filter_dict[taxa]
            if max(count_list) >= minimum_count:
                taxa_list.append(taxa)
        taxa_set = set(taxa_list)
        taxa_list = list(taxa_set)
        taxa_list.sort()
        # put counts in an array for plotting
        self.group_taxa_counts_list = []
        for group_name, group_list in zip(self.group_names, self.groups):
            for taxa in taxa_list:
                count_list = np.asarray(self.group_taxa_counts_dict[group_name].get(taxa, [0.0]*len(group_list)))
                average_count = np.mean(count_list)
                if average_count > 0.0:
                    std_dev_count = np.std(count_list)
                else:
                    std_dev_count = 0.0
                data_point = [group_name, taxa, average_count, std_dev_count]
                self.group_taxa_counts_list.append(data_point)
        if plot == "yes":
            self.make_bar_plot (self.group_taxa_counts_list, taxa_list, start_index = 0, error_bars = "yes",
                                fig_name = save_fig, y_max = max_y, x_label = "Taxa", y_label = data_type)

#------------------------------------------------------------------#

    def nog_category_counts (self, plot="no", data_type="spectra_counts", max_y = "", print_txt = "", save_fig = ""):
        #print("bactNOG category counts in proteomics sample\n\tCategory\tCount")
        cat_indexes = {'A':0,'B':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,
                       'J':9,'K':10,'L':11,'M':12,'N':13,'O':14,'P':15,'Q':16,
                       'R':17,'S':18,'T':19,'U':20,'V':21,'W':22,'X':23,'Y':24,'Z':25,'na':26,'u':27}
        cat_list = ['A','B','C','D','E','F','G','H','I','J','K','L','M',
                    'N','O','P','Q','R','S','T','U','V','W','X','Y','Z','na','u']
        # get counts for each nog category
        self.nog_category_counts_dict = {}
        for sample in self.samples:
            sample.nog_category_counts(data=data_type)
            self.nog_category_counts_dict[sample.name] = {}
            for cat in cat_list:
                count = sample.nog_category_count_dict[sample.name][cat]
                self.nog_category_counts_dict[sample.name][cat] = count
        # put counts in an array for plotting
        self.nog_category_counts_list = []
        for sample in self.samples:
            for cat in cat_list:
                data_point = [sample.name, cat, self.nog_category_counts_dict[sample.name][cat]]
                self.nog_category_counts_list.append(data_point)
        if plot == "yes":
            self.make_bar_plot (self.nog_category_counts_list, cat_list, start_index = 0, error_bars = "",
                                fig_name = save_fig, y_max = max_y, x_label = "bactNOG Categories", y_label = data_type)

#------------------------------------------------------------------#
    
    def group_average_nog_category_counts (self, plot="no", data_type="spectra_counts", max_y = "", print_txt = "", save_fig = ""):
        #print("bactNOG category counts in proteomics sample\n\tCategory\tCount")
        cat_indexes = {'A':0,'B':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,
                       'J':9,'K':10,'L':11,'M':12,'N':13,'O':14,'P':15,'Q':16,
                       'R':17,'S':18,'T':19,'U':20,'V':21,'W':22,'X':23,'Y':24,'Z':25,'na':26,'u':27}
        cat_list = ['A','B','C','D','E','F','G','H','I','J','K','L','M',
                    'N','O','P','Q','R','S','T','U','V','W','X','Y','Z','na','u']
        # get counts for each nog category
        self.group_average_nog_category_counts_dict = {}
        for g in range(len(self.group_names)):
            group_name = self.group_names[g]
            group_list = self.groups[g]
            self.group_average_nog_category_counts_dict[group_name] = {}
            for s in range(len(group_list)):
                sample = group_list[s]
                sample.nog_category_counts(data=data_type)
                for cat in cat_list:
                    count = sample.nog_category_count_dict[sample.name][cat]
                    if data_type == "spectra_counts":
                        count = sample.nog_category_count_dict[sample.name][cat]/sample.spectra_count
                    else:
                        count = sample.nog_category_count_dict[sample.name][cat]
                    if cat in self.group_average_nog_category_counts_dict[group_name]:
                        self.group_average_nog_category_counts_dict[group_name][cat].append(count)
                    else:
                        self.group_average_nog_category_counts_dict[group_name][cat] = [count]
        # put counts in an array for plotting
        self.group_average_nog_category_counts_list = []
        for group_name, group_list in zip(self.group_names, self.groups):
            for cat in cat_list:
                count_list = np.asarray(self.group_average_nog_category_counts_dict[group_name].get(cat, [0.0]*len(group_list)))
                average_count = np.mean(count_list)
                if average_count > 0.0:
                    std_dev_count = np.std(count_list)
                else:
                    std_dev_count = 0.0
                data_point = [group_name, cat, average_count, std_dev_count]
                self.group_average_nog_category_counts_list.append(data_point)
        if plot == "yes":
            self.make_bar_plot (self.group_average_nog_category_counts_list, cat_list, start_index = 0, error_bars = "yes",
                                fig_name = save_fig, y_max = max_y, x_label = "bactNOG Categories", y_label = data_type)
    
#------------------------------------------------------------------#
    
    def enrichment_counts (self, plot = "no", max_y = "", print_txt = "", save_fig = ""):
        enr_list = ['0--1','2--10','11--19','20--28','29--37','38--46','47--55','56--64','65--73','74--82','83--91','92--100']
        # get spectral counts for each enrichment bin
        self.enrichment_counts_dict = {}
        for sample in self.samples:
            sample.enrichment_counts(max_y = max_y)
            self.enrichment_counts_dict[sample.name] = {}
            for enr in enr_list:
                count = sample.enrichment_count_dict[sample.name][enr]
                self.enrichment_counts_dict[sample.name][enr] = count
        # put counts in an array for plotting
        self.enrichment_counts_list = []
        for sample in self.samples:
            for enr in enr_list[1:]:
                data_point = [sample.name, enr, self.enrichment_counts_dict[sample.name][enr]/sample.spectra_count]
                self.enrichment_counts_list.append(data_point)
        if plot == "yes":
            self.make_bar_plot (self.enrichment_counts_list, enr_list, start_index = 1, error_bars = "",
                                fig_name = save_fig, y_max = max_y , x_label = "%Enrichment Bins",
                                y_label = "Proportion of Balanced Spectral Counts")
    
#------------------------------------------------------------------#

    def group_average_enrichment_counts (self, plot = "no", max_y = "", print_txt = "", save_fig = ""):
        enr_list = ['0--1','2--10','11--19','20--28','29--37','38--46','47--55','56--64','65--73','74--82','83--91','92--100']
        # get spectral counts for each enrichment bin
        self.group_average_enrichment_counts_dict = {}
        for g in range(len(self.group_names)):
            group_name = self.group_names[g]
            group_list = self.groups[g]
            self.group_average_enrichment_counts_dict[group_name] = {}
            for s in range(len(group_list)):
                sample = group_list[s]
                sample.enrichment_counts(max_y = max_y)
                for enr in enr_list:
                    count = sample.enrichment_count_dict[sample.name][enr]/sample.spectra_count
                    if enr in self.group_average_enrichment_counts_dict[group_name]:
                        self.group_average_enrichment_counts_dict[group_name][enr].append(count)
                    else:
                        self.group_average_enrichment_counts_dict[group_name][enr] = [count]
        # put counts in an array for plotting
        self.group_average_enrichment_counts_list = []
        for group_name, group_list in zip(self.group_names, self.groups):
            for enr in enr_list[1:]:
                count_list = np.asarray(self.group_average_enrichment_counts_dict[group_name].get(enr, [0.0]*len(group_list)))
                average_count = np.mean(count_list)
                if average_count > 0.0:
                    std_dev_count = np.std(count_list)
                else:
                    std_dev_count = 0.0
                data_point = [group_name, enr, average_count, std_dev_count]
                self.group_average_enrichment_counts_list.append(data_point)
        if plot == "yes":
            self.make_bar_plot (self.group_average_enrichment_counts_list, enr_list, start_index = 1, error_bars = "yes",
                                fig_name = save_fig, y_max = max_y , x_label = "%Enrichment Bins",
                                y_label = "Proportion of Balanced Spectral Counts")

#------------------------------------------------------------------#
    """ 
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
        print("Use this data for:" +
              "\tPerforming PCoa" +
              "\t <> from skbio.stats.ordination import PCoA" +
              "\n\tPerforming permanova on NOG category distributions" +
              "\n\t command: permanova(data_eudm,groupList,permutations = 999)" +
              "\n\nAvailable <group>Lists:\n\tgenus, family, order, class, phylum" +
              "\n\n*** from skbio.stats.distance import permanova" +
              "\n*** from skbio.stats.distance import anosim")
        
            
    """

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