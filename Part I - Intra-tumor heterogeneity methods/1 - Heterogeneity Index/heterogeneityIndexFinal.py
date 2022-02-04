#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 23:34:53 2021

@author: david
"""


#####  Libraries  #####


# To open system files
import os

# To access all files of a specific type without using ifs
import glob

# Convert from a list of lists to a list
import itertools as it

# Numpy for data management
import numpy as np

# Pandas for data management
import pandas as pd

# Seaborn for plotting and styling
import seaborn as sns
sns.set_style("darkgrid")

# Matplotlib for additional customization
from matplotlib import pyplot as plt

# For additional customization
import matplotlib as mat

# Reverse the matplotlib settings back to default (inclunding colors of graphs)  
mat.rc_file_defaults()

# To have one color for each mice samples plot
import colorcet as cc

# To test run times
# import time

# start_time = time.time()
# print("--- %s seconds ---" % (time.time() - start_time))







#####  Constants  #####


# Change the values of the constants accordingly


# 0 to not include CNAs and 1 to consider the effect of CNAs in SNVs
USE_CNAS = 1

# if USE_CNAS = 1, when CNAS_ERRORS = 0 the CNAs with major copy number = 2 and minor copy number = 0 are used
# when CNAS_ERRORS = 1, the CNAs with major copy number = 2 and minor copy number = 0 are considered sequencing errors and are not used
# when USE_CNAS = 0, the CNAS_ERRORS constant has no effect and can be either 0 or 1
CNAS_ERRORS = 1




#--------------------------------------------------------------------------------------------------




#####  Methods  #####


def loadSNVs(filename): 
    """
    Load each mouse sample SNVs to a dataframe.
  
    Parameters  
            filename (str): the name of each csv file with a mouse sample SNVs to be load into a dataframe.

    Returns
            genomeInfo (dataframe): each mouse sample SNVs information.
    """

    genomeInfo = pd.read_csv(filename, delimiter=';',  dtype = {'CHROM': 'category', 'POS': 'uint32', 'FILTER': 'category', 'TUMOR': str},
                             names = ['CHROM', 'POS', 'FILTER', 'TUMOR'], header = 0, usecols = [0, 1, 6, 9])

    # SNVs in chr X will not be included because of CNAs caller methods errors (with some CNAs with major and minor copy number = 0 in chr X) 
    genomeInfo = genomeInfo[(genomeInfo['FILTER'] == 'PASS') & (genomeInfo['CHROM'] != 'chrX')]
  
    # remove the filter column
    genomeInfo.drop(columns = ['FILTER'], axis=1, inplace=True)
  
    # reset indexes to start on 0 and increment by 1
    genomeInfo.index = range(len(genomeInfo.index))
  
    
    return genomeInfo
  
    
  

def loadCNAs(filename):
    """
    Load each mouse sample CNAs to a dataframe.

    Parameters
            filename (str): the name of each csv file with a mouse sample CNAs to be load into a dataframe.

    Returns
            cnasAttributes (dataframe): each mouse sample CNAs information.
    """

    # skipfooter = 1, remove last line of each dataframe, remove chrY since the data is about breast cancer cell line (women only have the X chromosome)
    # engine = "python" to be able to use skipfooter
    cnasAttributes = pd.read_csv(filename, delimiter=';',  
                                 dtype = {'CHROM': 'category', 'START_POS': 'uint32', 'END_POS': 'uint32', 'MAJOR_COPY' : pd.UInt8Dtype(), 'MINOR_COPY' : pd.UInt8Dtype()},
                                 names = ['CHROM', 'START_POS', 'END_POS', 'MAJOR_COPY', 'MINOR_COPY'], header = 0,
                                 usecols = [0, 1, 2, 8, 9], skipfooter = 1, engine = "python")


    # drop rows with nan values
    cnasAttributes = cnasAttributes.dropna() 
        
    # change data type to reduce memory usage
    cnasAttributes['MAJOR_COPY'] = cnasAttributes['MAJOR_COPY'].astype('uint8')
    cnasAttributes['MINOR_COPY'] = cnasAttributes['MINOR_COPY'].astype('uint8')
    
    
    # eliminate CNAs sequencing errors
    cnasAttributes = cnasAttributes.drop(cnasAttributes[(cnasAttributes['MAJOR_COPY'] == 0) & (cnasAttributes['MINOR_COPY'] == 0)].index)
    
    # CNAs with major copy number 2 and minor copy number 0 will not be used (considered as possible sequencing errors)
    if CNAS_ERRORS == 1:
        cnasAttributes = cnasAttributes.drop(cnasAttributes[((cnasAttributes['MAJOR_COPY'] == 2) & (cnasAttributes['MINOR_COPY'] == 0))].index)
    
    
    # reset indexes to start on 0 and increment by 1
    cnasAttributes.index = range(len(cnasAttributes.index))
  

    return cnasAttributes




def removeNonDiploidSNVs(sampleSNVsInfo, sampleCNAsInfo):
    """
    Remove the SNVs that are in regions with CNAs with non-normal copy numbers in each mouse sample.

    Parameters
            sampleSNVsInfo (dataframe): the SNVs of a specific sample.
            sampleCNAsInfo (dataframe): the CNAs of a specific sample.
                
    Returns
            diploidSNVsOnly (dataframe): only SNVs that are in diploid regions, since here we cannot assess the effect of CNAs on these, so we just remove them.
    """

    # Intersect the SNVs and CNAs of a same mouse sample (how = 'left' in the case a SNV is not intersected by any CNA, all the SNVs are kept)
    snvsIntersectedByCNAs = pd.merge(sampleSNVsInfo, sampleCNAsInfo, on = ['CHROM'], how = 'left')
    
    
    # Get the chromosome and position of the SNVs to remove
    # Sometimes, a same SNV in a same chromosome is not intersected by a CNA, but it is intersected by other, so the SNV should be remove (it is a non-diploid region)
    # Only SNVs intersected by CNAs with major and minor copy number that are not one are removed
    snvsToRemove = snvsIntersectedByCNAs[((snvsIntersectedByCNAs['START_POS']) <= (snvsIntersectedByCNAs['POS'])) & \
                                         ((snvsIntersectedByCNAs['POS']) <= (snvsIntersectedByCNAs['END_POS'])) &
                                         ((snvsIntersectedByCNAs['MAJOR_COPY'] != 1) | (snvsIntersectedByCNAs['MINOR_COPY'] != 1))] [['CHROM', 'POS']] 
                                         

    snvsToRemove = (snvsToRemove['CHROM'].astype(str) + snvsToRemove['POS'].astype(str)).unique().tolist()

    # find which SNVs have the same chromosome and position as the ones to remove
    # for example, an entry of a SNV to remove with a CNA in the same chromosome that did not intersect it
    snvsToRemove = np.array(snvsIntersectedByCNAs.loc[ \
            ((snvsIntersectedByCNAs['CHROM'].astype(str) + snvsIntersectedByCNAs['POS'].astype(str)).isin(snvsToRemove))].index).astype('uint32')
    
    
    # drop the SNVs in non-diploid regions
    diploidSNVsOnly = snvsIntersectedByCNAs.drop(snvsToRemove).drop_duplicates(subset = ['CHROM', 'POS'])

    diploidSNVsOnly.drop(columns = ['START_POS', 'END_POS', 'MAJOR_COPY', 'MINOR_COPY'], axis=1, inplace=True)

      
    # reset indexes to start on 0 and increment by 1
    diploidSNVsOnly.index = range(len(diploidSNVsOnly.index))


    return diploidSNVsOnly




def getDFsToPlot():
    """
    Create the heterogeneity indexes, SNVs sequencing total coverage and SNVs VAFs dataframes to be plot.
    
    Returns
            miceHI (dataframe): the dataframe with the heterogeneity indexes values per sample.
            totalCoverage (dataframe): the dataframe with total coverage of each SNV and its sample.
            samplesSNVsVAF (dataframe): the dataframe with the VAF of each SNV and its sample.
    """
    
    groups = ['Ctrl', 'Sunit']

    miceID = ['Mouse 49', 'Mouse 55', 'Mouse 61', 'Mouse 62']

    regionsID = ['Region 1', 'Region 2', 'Region 3', 'Region 4']


    # To identify each mouse sample in the total coverage and VAF plots, with the format M'xx'R'y', where xx = 49, 55, 61 or 62, and y = 1, 2, 3 or 4
        
    sampleIDs = list(it.chain.from_iterable([[mouse[0] + mouse[-2:] + region[0] + region[-1:] for region in regionsID ] for mouse in miceID]))

    sampleIDs = list(it.chain.from_iterable([[sample] * len(sampleSNVs) for sample, sampleSNVs in zip (sampleIDs, snvsInfo)]))


    # Get the reference and alternative alleles read counts of each SNV.   
    allelesReads = [list(map(lambda alleleInfo: \
                     (int(alleleInfo[alleleInfo.find(':') + 1 : alleleInfo.find(',')]), 
                        int(alleleInfo[alleleInfo.find(',') + 1 : alleleInfo.find(':', alleleInfo.find(','))])), \
                           list(sampleSNVsReads['TUMOR']))) for sampleSNVsReads in snvsInfo]


    # SNVs reads total coverage dataframe
    
    samplesTotalCoverage = [[snvReads[0] + snvReads[1] for snvReads in snvSampleReads] for snvSampleReads in allelesReads]
                                 
    totalCoverage = pd.DataFrame(list(zip(list(it.chain.from_iterable(samplesTotalCoverage)), sampleIDs)), columns = ['Total coverage', 'Sample'])
                                
    
    # Get the variant allele frequencies of each SNV based on the alleles read counts
    vafs = [[(readsInfo[1] / (readsInfo[0] + readsInfo[1])) * 100 for readsInfo in snvSampleReads] for snvSampleReads in allelesReads]


    # SNVs variant allele frequencies dataframe
    
    samplesSNVsVAF = pd.DataFrame(list(zip(list(it.chain.from_iterable(vafs)), sampleIDs)), columns = ['Variant Allele Frequency', 'Sample'])
    
    
    
    # Estimating the heterongeity index using the Shannon's index
        

    # take only the counts (number of SNVs in a bin) and not bin_edges (using 10 bins)
    # numpy.histogram - compute the histogram of a dataset. 
    # it returns two values: first an array with the values of the histogram and then the bin edges based on the min and max VAF values divided in 10 intervals
    miceSampleCounts = [list(np.histogram(snvVAF, bins = 10)[0]) for snvVAF in vafs]
     
    
    # probability of a SNV being in each of the bins of a mouse sample.
    # len(tumorInfo[miceSampleCounts.index(mouseCounts)), the number of SNVs in a sample, is the same as the sum of the SNVs in each bin
    snvBinsProbability = [[snvBinCount/sum(mouseCounts) for snvBinCount in mouseCounts] for mouseCounts in miceSampleCounts]
    
    
    # Calculated based on the population diversity Shannon's Index
    # The if condition is to account for np.log(0) which is infinite, since there are some bins which are empty, do not have any SNVs
    heterogeneityIndices = [-sum(list(map(lambda snvBinsProbability: 0 if snvBinsProbability == 0.0 else snvBinsProbability * np.log(snvBinsProbability), sampleSNVBinsProbability))) \
                                for sampleSNVBinsProbability in snvBinsProbability]
    
    
    # The dataframe with the heterogeneity indices information about each mouse sample 
    # With the corresponding groups: Mouse 49 and 55 - Ctrl, Mouse 61 and 62 - Sunit.
    miceHI = pd.DataFrame(list(zip(heterogeneityIndices, \
                            list(it.chain.from_iterable([[miceID[i]] * 4 for i in range(0, len(miceID))])), regionsID * 4, [groups[0]] * 8 + [groups[1]] * 8)), \
                              columns = ['TH index', 'Mouse', 'Region', 'Group'])


    return miceHI, totalCoverage, samplesSNVsVAF




def plotTotCoverageAllSamples():
    """
    Prints the KDE plot with the indication of the total coverage which corresponds to the maximum density using all samples.
    """
    
    # KDE plot with the normal distribution of the SNVs coverage in all samples
    g = sns.displot(data = snvsCoverage, x = "Total coverage", kind = "kde", cut = 0)
    
    
    # Get all the data points of the plot
    data_points = g.fig.get_children()[1].get_lines()[0].get_xydata()
    
    
    # Get the x and y points corresponding to the highest y value, in this case, the kernel density 
    max_xy = data_points[np.argmax(data_points[:, 1])]
    
    # Round the maximum coverage x-axis value to 2 decimal places
    maxDensity = "Total coverage with  \nmax density is {max_d}".format(max_d = round(max_xy[0], 2))
    
    
    # Annotate the figure with the indication of which is the total coverage value corresponding to the maximum kernel density
    g.fig.get_children()[1].annotate(maxDensity, fontsize = 14, xy=(max_xy[0], max_xy[1]), xycoords='data',
                                     xytext=(0.92, 0.952), textcoords='axes fraction',
                                     arrowprops=dict(facecolor='black', shrink=0.07),
                                     horizontalalignment='right', verticalalignment='center')
    
    figure_name = "Genome sequencing data reads total coverage distribution using all samples"
    
    
    # Define axes labels size
    for ax in g.axes.flat:
        
        
        xlabel = ax.get_xlabel()
        
        ylabel = ax.get_ylabel()
        
                    
        ax.set_xlabel(xlabel, fontsize = 14, labelpad = 8)
        ax.set_ylabel(ylabel, fontsize = 14, labelpad = 8)
        
        # default of axis is both
        ax.tick_params(labelsize = 13)
    
    
    plt.suptitle(figure_name, fontsize = 15, y = 1.09).set_weight("bold")
    
    
    # save figure
    
    figure_name = figure_name + '.png'
      
    plt.savefig(figure_name, bbox_inches = 'tight')    
    
    plt.show()




def plotVafAllSamples():
    """
    Prints the KDE plot with the normal distribution and the histogram of the variant allele frequency of the SNVs using all samples
    """

    # KDE plot with the normal distribution of the SNVs VAF in all samples
    g = sns.displot(data = snvsVAFs, x = 'Variant Allele Frequency', kind = "kde", cut = 0)


    # Define axes labels size    
    for ax in g.axes.flat:
        
        
        xlabel = ax.get_xlabel()
        
        ylabel = ax.get_ylabel()
        
                    
        ax.set_xlabel(xlabel, fontsize = 14, labelpad = 8)
        ax.set_ylabel(ylabel, fontsize = 14, labelpad = 8)
        
        # default of axis is both
        ax.tick_params(labelsize = 13)
    
    
    figure_name = "Variant Allele Frequency distribution using all samples"
    
    plt.suptitle(figure_name, fontsize = 15, y = 1.08).set_weight("bold")
    
    
    # save figure
    
    figure_name = figure_name + '.png'
      
    plt.savefig(figure_name, bbox_inches = 'tight')    
    
    plt.show()
    
    
    

def plotSamplesCoverageDifPlots():
    """
    Prints the KDE plot with the indication of the total coverage which corresponds to the maximum density for each sample in different plots.
    """
    
    # KDE plot with the normal distribution of the SNVs coverage in different samples output in different plots
    g = sns.displot(snvsCoverage, x = 'Total coverage', kind = "kde", cut = 0, col = "Sample", col_wrap = 4, aspect = 1.2, \
                    height = 5.6, legend = False, facet_kws={'sharey' : False, 'sharex' : False})
        
        
    # Figure and subplots settings    
    g.fig.set_size_inches(21, 12)
        
    g.fig.subplots_adjust(hspace = 1.0)
      
    g.fig.subplots_adjust(wspace = 0.45) 
    
    
    # Set the tittles for each of the subplots
    g.set_xlabels("Total coverage").set_titles("Sample {col_name}")  
    
    
    # For each of the subplots
    for ax in g.fig.get_children()[1:]:
        
        # Get all the data points of that subplot
        data_points = ax.get_lines()[0].get_xydata()
        
        
        # Get the x and y points corresponding to the highest y value, in this case, the kernel density 
        max_xy = data_points[np.argmax(data_points[:, 1])]
        
        # Round the maximum coverage x-axis value to 2 decimal places
        maxDensity = "Total coverage with  \nmax density is {max_d}".format(max_d = round(max_xy[0], 2))
        
    
        # Annotate each subplot with the indication of which is the total coverage value corresponding to the maximum kernel density
        ax.annotate(maxDensity, fontsize = 16, xy=(max_xy[0], max_xy[1]),  xycoords='data',
                    xytext=(1.00, 0.952), textcoords='axes fraction',
                    arrowprops=dict(facecolor='black', shrink=0.07),
                    horizontalalignment='right', verticalalignment='center')
    
    
    figure_name = "Genome sequencing data reads total coverage distribution by sample"
    
    
    # Define axes labels size
    for ax in g.axes.flat:
            
        xlabel = ax.get_xlabel()
        
        ylabel = ax.get_ylabel()
        
        
        if xlabel == '':
            xlabel = 'Total coverage'
    
    
        if ylabel == '':
            ylabel = 'Density'
    
    
        ax.set_xlabel(xlabel, fontsize = 15, labelpad = 8)
        ax.set_ylabel(ylabel, fontsize = 15, labelpad = 8)
    
        ax.set_title(ax.get_title(), fontsize = 16, y = 1.20).set_weight("bold")
        
        
        # default of axis is both
        ax.tick_params(labelsize = 14)
    
    
    plt.suptitle(figure_name, fontsize = 18, y = 1.10).set_weight("bold")
    
    
    # save figure
    
    figure_name = figure_name + '.png'
      
    plt.savefig(figure_name, bbox_inches = 'tight')    
    
    plt.show()



    
def plotSamplesCoverageOnePlot():
    """
    Prints the KDE plot with the indication of the total coverage which corresponds to the maximum density for each sample output in the same plot.
    """
  
    # KDE plot with the normal distribution of the SNVs coverage in different samples output in the same plot
    g = sns.displot(data = snvsCoverage, x = "Total coverage", kind = "kde", cut = 0, legend = False, 
                    hue = "Sample", hue_order = list(snvsCoverage['Sample'].unique())[::-1], 
                    palette = sns.color_palette(cc.glasbey, n_colors=16)[::-1], height = 8, aspect = 1.5)
    
    
    figure_name = "Genome sequencing SNVs reads coverage distribution of different samples in one plot"
    
    
    # Define axes labels size
    for ax in g.axes.flat:
        
        
        xlabel = ax.get_xlabel()
        
        ylabel = ax.get_ylabel()
        
                    
        ax.set_xlabel(xlabel, fontsize = 14.5, labelpad = 10)
        ax.set_ylabel(ylabel, fontsize = 14.5, labelpad = 10)
        
        # default of axis is both
        ax.tick_params(labelsize = 14)
    
    
    # Figure settings
    g.fig.set_size_inches(6, 6.8)
    
    
    sampleLegend = plt.gca().legend(bbox_to_anchor=(1.04, 0.47), loc='center left', borderaxespad=1.5, fontsize = 14, \
                                    labels = list(snvsCoverage['Sample'].unique()), title = r'$\bf{Sample}$', title_fontsize = 14, \
                                    shadow = True, facecolor = 'white', borderpad = 0.9, labelspacing = 1.0, handletextpad = 0.7,
                                    handlelength = 2.5)    
    
        
    # Increase legend title pad from legend labels
        
    legendTitle = sampleLegend.get_title()
    
    legendTitle.set_y(3)
    
    
    # Increase legend handles width
    for legHandle in sampleLegend.legendHandles:
        legHandle.set_linewidth(3.5)
    
            
    plt.suptitle(figure_name, fontsize = 15, y = 1.08).set_weight("bold")
    
         
    # save figure
    
    figure_name = figure_name + '.png'
      
    plt.savefig(figure_name, bbox_inches = 'tight')    
    
    plt.show()
    



def plotSamplesVAFDifPlots():
    """
    Prints the KDE plot with the normal distribution of the variant allele frequency of the SNVs for each sample
    """
    
    # KDE plot with the normal distribution of the SNVs VAF in each sample output in the same plot
    g = sns.displot(data = snvsVAFs, x = 'Variant Allele Frequency', kind = "kde", cut = 0, col = "Sample", col_wrap = 4, \
                    aspect = 1.4, height = 5.6, legend = False, facet_kws={'sharey' : False, 'sharex' : False})
        
        
    # Figure and subplots settings       
    g.fig.set_size_inches(21, 12)
        
    g.fig.subplots_adjust(hspace = 0.8)
      
    g.fig.subplots_adjust(wspace = 0.25) 

    
    # Set the titles for each of the subplots
    g.set_xlabels("Variant Allele Frequency").set_titles("Sample {col_name}")  
    
    
    figure_name = "Variant Allele Frequency distribution by sample"
    
    
    # Define axes labels size
    for ax in g.axes.flat:
        
        xlabel = ax.get_xlabel()
        
        ylabel = ax.get_ylabel()
        
        
        if xlabel == '':
            xlabel = 'Variant Allele Frequency'
    
    
        if ylabel == '':
            ylabel = 'Count'
    
    
        ax.set_xlabel(xlabel, fontsize = 15, labelpad = 8)
        ax.set_ylabel(ylabel, fontsize = 15, labelpad = 8)
    
        ax.set_title(ax.get_title(), fontsize = 16, y = 1.04).set_weight("bold")
        
        
        # default of axis is both
        ax.tick_params(labelsize = 14)
    
    
    plt.suptitle(figure_name, fontsize = 18, y = 1.08).set_weight("bold")
    
    
    # save figure
    
    figure_name = figure_name + '.png'
      
    plt.savefig(figure_name, bbox_inches = 'tight')    
    
    plt.show()




def plotSamplesVAFOnePlot():   
    """
    Prints the KDE plot with the normal distribution of the variant allele frequency of the SNVs from each sample output in one plot
    """
    
    # KDE plot with the normal distribution of the SNVs VAF in each sample with each output in one plot simultaneously
    g = sns.displot(data = snvsVAFs, x = 'Variant Allele Frequency', kind = "kde", cut = 0, hue = "Sample", legend = False, \
                    hue_order = list(snvsVAFs['Sample'].unique())[::-1], palette = sns.color_palette(cc.glasbey, n_colors=16)[::-1],)
    
    
    # Define axes labels size
    for ax in g.axes.flat:
        
        
        xlabel = ax.get_xlabel()
        
        ylabel = ax.get_ylabel()
        
                    
        ax.set_xlabel(xlabel, fontsize = 14.5, labelpad = 10)
        ax.set_ylabel(ylabel, fontsize = 14.5, labelpad = 10)
        
        # default of axis is both
        ax.tick_params(labelsize = 14)
    
    
    figure_name = "Variant Allele Frequency distribution of the different samples SNVs in one plot"
    
    
    # Figure settings
    g.fig.set_size_inches(6, 6.8)
    
    
    sampleLegend = plt.gca().legend(bbox_to_anchor=(1.04, 0.47), loc='center left', borderaxespad=1.5, fontsize = 14, \
                                    labels = list(snvsVAFs['Sample'].unique()), title = r'$\bf{Sample}$', title_fontsize = 14, \
                                    shadow = True, facecolor = 'white', borderpad = 0.9, labelspacing = 1.0, handletextpad = 0.7,
                                    handlelength = 2.5)    
    
        
    # Increase legend title pad from legend labels
        
    legendTitle = sampleLegend.get_title()
    
    legendTitle.set_y(3)
    
    
    # Increase legend handles width
    for legHandle in sampleLegend.legendHandles: 
        legHandle.set_linewidth(3.5) 
    
            
    plt.suptitle(figure_name, fontsize = 15, y = 1.09).set_weight("bold")
    
         
    # save figure
    
    figure_name = figure_name + '.png'
      
    plt.savefig(figure_name, bbox_inches = 'tight')    
    
    plt.show()
    
    
        
        
def plotHeterogeneityIndex():
    """
    Plot the heterogeneity index values for each of the mice samples
    """
    
    # if only SNVs are used or if CNAs are also used
    if USE_CNAS == 0:
        figure_name = "TH index comparison between experiment groups (only SNVs)"
         
    else:    
        # If CNAs with major copy 2 and minor copy 0 are considered sequencing errors
        if CNAS_ERRORS == 0:
            figure_name = "TH index comparison between experiment groups (SNVs and CNAs)"    
        else:
            figure_name = "TH index comparison between experiment groups (SNVs and CNAs, excluding possible sequencing errors)"
            
        
    # Plot the mice data points colored by group type
    g = sns.catplot(x = "Mouse", y = "TH index", hue = "Group", palette=sns.color_palette("tab10", 2), \
        data = miceSamplesHI, legend = False, sharex = False, sharey = False, aspect = 1.2) 
    
    
    # set the axes labels
    g.set_axis_labels("Mouse", "Tumor Heterogeneity index")
        
        
    # Define axes labels size
    for ax in g.axes.flat:
        
        xlabel = ax.get_xlabel()
        
        ylabel = ax.get_ylabel()
        
        ax.set_xlabel(xlabel, fontsize = 15, labelpad = 13)
        ax.set_ylabel(ylabel, fontsize = 15, labelpad = 13)
    
    
        # default of axis is both
        ax.tick_params(labelsize = 13)
        
        
        # Itearate through the mice samples data points
        for mp in ax.get_children()[0:4]:
    
            # increase data points size
            mp.set_sizes([42,42,42,42])
        
    
    groupLegend = plt.gca().legend(bbox_to_anchor=(1.03, 0.5), loc='center left', borderaxespad=1.0, fontsize = 13, labels = ['Control', 'Sunitinib'], title = r'$\bf{Group}$',
                                   title_fontsize = 14, shadow = True, facecolor = 'white', borderpad = 0.7, labelspacing = 0.6, handletextpad = 0.2)    
    
    
    # same color as data points
    groupLegend.legendHandles[0].set_color(plt.gca().get_children()[1].get_facecolor()[0:3])
    groupLegend.legendHandles[1].set_color(plt.gca().get_children()[2].get_facecolor()[0:3])
    
    
    # set markers' sizes
    groupLegend.legendHandles[0].set_sizes([10])
    groupLegend.legendHandles[1].set_sizes([10])
    
    
    # line width of legend handles    
    for legHan in groupLegend.legendHandles:
        legHan.set_linewidth(4.5)
        
        
    # Increase legend title pad from legend labels
        
    legendTitle = groupLegend.get_title()
    
    legendTitle.set_y(3)
       

    plt.suptitle(figure_name, fontsize = 16, y = 1.10).set_weight("bold")
    
    
    # save figure
    
    figure_name = figure_name + '.png'
      
    plt.savefig(figure_name, bbox_inches = 'tight')    
    
    plt.show()
        


    
    
    
    
#--------------------------------------------------------------------------------------------------




#####  Load SNVs and CNAs input data  #####


# Current directory: Dados/variantes_somaticas/Heterogeneity Index


os.chdir("../excel_format")


# load the SNVs
snvsInfo = list(map(loadSNVs, sorted(glob.iglob("*.csv", recursive=False))))




# if true, we do not use CNAs and only SNVs
if USE_CNAS == 1:
    
    
    os.chdir("../../CNVS/excel_format_germline")
    

    # load the CNAs
    cnasInfo = list(map(loadCNAs, sorted(glob.iglob("*.csv", recursive=False))))


    # Intersect SNVs and CNAs, removing only the SNVs in non-diploid regions, given that we can not assess the effect of CNAs in SNVs with an algorithm using the heterogeneity index
    # Zones with copy number = 2 and major and minor copy numbers = 1 are normal diploid, and so in these cases SNVs are not removed
    snvsInfo = list(map(removeNonDiploidSNVs, snvsInfo, cnasInfo))
    
    
    
    
    
    
    
#--------------------------------------------------------------------------------------------------

    


#####  Dataframes of the heterogeneity indexes, total coverage and the VAFs of SNVs to be used for plotting  #####


miceSamplesHI, snvsCoverage, snvsVAFs = getDFsToPlot()




#--------------------------------------------------------------------------------------------------




# Change to the images directory where the output plots will be

if USE_CNAS == 0: 
    # if not using CNAs
    os.chdir("../Heterogeneity Index/Images")

else: 
    os.chdir("../../variantes_somaticas/Heterogeneity Index/Images")


# The analysis of the SNVs total coverage and VAFs distributions using both all samples or each sample differently are not affected by using only SNVs or also CNAs
# The difference in the number of SNVs used is not significative from case to case




#####  Data analysis using all samples  #####


# Plotting the genome sequencing data reads total coverage distribution using all samples for the distribution

plotTotCoverageAllSamples()


# Plotting the variant allele frequency of SNVs using all samples for the distribution 

plotVafAllSamples()




#####  Data analysis by sample  #####


# Plotting the genome sequencing data reads total coverage distribution by sample in different plots

plotSamplesCoverageDifPlots()


# The different sample coverage distributions are plot in the same axes       

plotSamplesCoverageOnePlot()


# Plotting the variant allele frequency of SNVs by sample in different plots

plotSamplesVAFDifPlots()


# The different VAFs distibutions are plot in the same axes       

plotSamplesVAFOnePlot()







#####  Plotting the mice samples heterogeneity indexes  #####


# Plots only using SNVs, using SNVs and all CNAs, and using SNVs and CNAs, but considering the CNAs with major copy 2 and minor copy 0 as segment errors
# The different plots are based on the constants defined values

plotHeterogeneityIndex()






