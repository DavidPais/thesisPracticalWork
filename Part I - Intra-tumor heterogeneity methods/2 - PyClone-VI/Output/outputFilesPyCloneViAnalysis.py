#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 19:31:18 2021

@author: david
"""


#####  Libraries  #####


# To open system files
import os

# To access all files of a specific type without using ifs
import glob

# To access regular expressions methods
import regex as re

# Pandas for data management
import pandas as pd

# Mathematical functions 
import numpy as np 

# Seaborn for plotting and styling
import seaborn as sns
sns.set_style("ticks")

# Matplotlib for additional customization
from matplotlib import pyplot as plt

# For additional customization
import matplotlib as mat

# Reverse the matplotlib settings back to default (inclunding colors of graphs)  
mat.rc_file_defaults()

# For the QQ plot
import statsmodels.api as sm

# For the Wilcoxon rank-sum test (or Mann-Whitney U test)
from scipy.stats import mannwhitneyu

# To keep the precision in operations between floats
from decimal import Decimal

# To test run times
# import time

# start_time = time.time()
# print("--- %s seconds ---" % (time.time() - start_time))




#--------------------------------------------------------------------------------------------------




#####  Methods  #####


def processOutputFiles(filename, withoutCNAsErrors):
    """
    Load and process the information of the output files from PyClone-VI.

    Parameters
            filename (str): the name of the PyClone-VI execution file to be processed.
            withoutCNAsErrors (bool): whether or not to consider CNAs errors.

    Returns
            pyCloneViExecution (dataframe): the processed PyClone-VI execution file information.
    """
    
    
    pyCloneViExecution = pd.read_csv(filename, delimiter = '\t', dtype = {'mutation_id' : 'category', 'sample_id' : 'category', 'cluster_id' : 'category',
                                                                          'cellular_prevalence' : 'float16', 'cellular_prevalence_std' : 'float16', 
                                                                          'cluster_assignment_prob' : 'float16'})


    # Define the execution_id, mouse, group and sample columns
    pyCloneViExecution['execution_id'], pyCloneViExecution['mouse'], pyCloneViExecution['group'], pyCloneViExecution['sample'] = \
        [filename[:-4], pyCloneViExecution['sample_id'].str[:3].astype('category'), \
                     pyCloneViExecution['sample_id'].str[3 : -2].astype('category'), pyCloneViExecution['sample_id'].str[-2:].astype('category')]

    # originally, filename in the format: ".\x\y", where, e.g., x = ".\allMice\" and y = "allMice_snvsInfo5c10r", we only want to keep y
    # get everything after last slash:
    pyCloneViExecution['execution_id'] = pyCloneViExecution['execution_id'].str.extract(r'([^\\]+$)', expand = False)
    
    

    # Different processing of num_clusters and num_restarts columns based on the file name including or not CNAs errors
    if withoutCNAsErrors == False:
    
        pyCloneViExecution[['num_clusters', 'num_restarts']] = pyCloneViExecution['execution_id'].str.extract(r'o(.*?)r', expand = False).str.split("c", expand = True)

    else:
        
    # extractall results in a multi index dataframe, with the level match having the different matches for each row
    # column name defined in the extractall with ?P<parameters> where the pattern is found
    # xs selects a particular level of the multi index, and the entry we want first (the 1)
    # select the column and apply the split
    
        pyCloneViExecution[['num_clusters', 'num_restarts']] = pyCloneViExecution['execution_id'].str.extractall(r's(?P<parameters>.*?)r').xs(1, level = 'match') \
                ['parameters'].str.split("c",expand = True)


    
    # select only those mutations which have more than 60% certainty of belonging to a cluster
    pyCloneViExecution = pyCloneViExecution[pyCloneViExecution['cluster_assignment_prob'] > 0.6]


    # change data types to reduce memory usage
    pyCloneViExecution = pyCloneViExecution.astype({'execution_id': 'category', 'num_clusters': 'uint8', 'num_restarts' : 'uint16'}) 
   
    # order columns 
    pyCloneViExecution = pyCloneViExecution[['execution_id', 'mutation_id', 'sample_id', 'mouse', 'sample', 'group', 'num_clusters', 'num_restarts', \
                                             'cluster_id', 'cellular_prevalence', 'cellular_prevalence_std', 'cluster_assignment_prob']]
    

    # sort the dataframe according the sample and the number of cluster 
    pyCloneViExecution.sort_values(['sample_id','cluster_id'], inplace = True)
        
    
    # reset indexes to start on 0 and increment by 1
    pyCloneViExecution.index = range(len(pyCloneViExecution.index))

    
    return pyCloneViExecution




def miceStatistics(executionSettingsSample, separatedMouseFiles):
    """
    Process statistics about the PyClone-VI execution data files including the cluster mean cellular prevalence and the samples heterogeneity index based on these.

    Parameters
            executionSettingsSample (dataframe): the execution dataframe being processed.
            separatedMouseFiles (bool): if all mice together in a file or different mouse files dataframe lists are being processed.

    If separatedMouseFiles == False:
        Returns
                miceClusterCCFs (dataframe): the mean cellular prevalence of each cluster in each sample of an execution file.
                miceSamplesHI (dataframe): the heterogeneity index per sample based on the clusters mean cellular prevalence values.
            
    Else:
        Returns
                miceNumClustersObtained (dataframe): the number of clusters that were obtained with each different PyClone-VI execution settings.
                miceSamplesHI (dataframe): the heterogeneity index per sample based on the clusters mean cellular prevalence values.
    """
    
    
    # observed = True to not include all categories in all rows, that is, all combinations, just the ones that really exist in each sample
    miceClusterCCFs = executionSettingsSample.groupby(['execution_id', 'sample_id', 'mouse', 'sample', 'group', 'num_clusters', 'num_restarts', 'cluster_id'], \
                                                        as_index = False, observed = True).agg({'cellular_prevalence' : 'mean'})

    miceClusterCCFs.rename(columns = {'cellular_prevalence' : 'mean_cellular_prevalence'}, inplace = True)

    miceClusterCCFs = miceClusterCCFs.astype({'num_clusters' : 'uint8', 'num_restarts' : 'uint16'})



    # the heterogeneity index was calculated based on the CCF in absolute value, and then the CCF was used in percentage for representation
    miceSamplesHI = miceClusterCCFs.groupby(['execution_id', 'sample_id', 'mouse', 'sample', 'group', 'num_clusters', 'num_restarts'], as_index = False, observed = True) \
        ['mean_cellular_prevalence'].apply( \
       
            # for each sample clusters group of groupby, get the values of each cluster mean cellular prevalence to calculate the heterogeneity index
            lambda sampleClusters: -sum(list(map( \
                lambda clusterMeanPrevallence: 0 if clusterMeanPrevallence == 0.0 \
                    else clusterMeanPrevallence * np.log(clusterMeanPrevallence), list(sampleClusters)))))
          
    miceSamplesHI.rename(columns = {'mean_cellular_prevalence' : 'heterogeneity_index'}, inplace = True)
  
    miceSamplesHI = miceSamplesHI.astype({'num_clusters' : 'uint8', 'num_restarts' : 'uint16', 'heterogeneity_index' : 'float32'})
                                                                                                                            
                                                                                                                                      
                                                                                             
    if separatedMouseFiles == True:
        
        miceNumClustersObtained = miceClusterCCFs.groupby(['execution_id', 'sample_id', 'mouse', 'sample', 'group', 'num_clusters', 'num_restarts'], \
                                   as_index = False, observed = True).agg({'cluster_id' : 'count'})                                                                                                                             
                                                                                                                                          
        miceNumClustersObtained.rename(columns = {'cluster_id' : 'output_clusters_obtained'}, inplace = True)
         
        # the num_clusters is the input setting used to run the PyClone-VI
        # the output_clusters_obtained is the number of clusters obtained for each setting and mouse sample after PyClone-VI execution
        miceNumClustersObtained = miceNumClustersObtained.astype({'num_clusters' : 'uint8', 'num_restarts' : 'uint16', 'output_clusters_obtained' : 'uint8'})
                                                                                                                                 
    
        # We can not compare the different clusters CCF when using different files for each mouse
        return miceNumClustersObtained, miceSamplesHI

 
    return miceClusterCCFs, miceSamplesHI




def allMiceClusterCCFsDistributionPlot(allMiceTypeDF, outputDir):
    """
    Plot the cellular_prevalence (or CCFs since purity = 1) distribution of the mutations belonging to each of the mice clusters, using all mice together.

    Parameters
            allMiceTypeDF (dataframe): all the different mice mutations CCFs dataframe type (with or without considering CNAs errors).
            outputDir (str): the directory to where the plots will be output.
    """
    
    
    os.chdir(outputDir)
    
    # all together in one dataframe
    allMiceTypeDF = pd.concat(allMiceTypeDF, ignore_index = True)
    
    
    # Converting the CCF to percentage
    allMiceTypeDF['cellular_prevalence'] = allMiceTypeDF['cellular_prevalence'] * 100
    
    
    # Plot for each execution settings
    for execution in allMiceTypeDF['execution_id'].unique():
            
        miceExecSettings = allMiceTypeDF[allMiceTypeDF['execution_id'] == execution]
        
        
        # KDE plot with the normal distribution of the SNVs CCFs of each experimental group in all samples
        g = sns.displot(miceExecSettings, x = "cellular_prevalence", hue = "group", col = "cluster_id", kind = "kde", cut = 0, \
                        palette=sns.color_palette("tab10", 2), legend = False, facet_kws={'sharey' : False, 'sharex' : False} ) 
      
        
        g.set_xlabels("Cluster Mutations Cellular Fraction").set_titles("Cluster {col_name}")  
            
        g.fig.tight_layout()   
    
        g.fig.subplots_adjust(wspace = 0.55)
    
        
        # Define axes labels size
        for ax in g.axes.flat:
            
            xlabel = ax.get_xlabel()
            
            ylabel = ax.get_ylabel()
            
            
            if(xlabel == ''):
                xlabel = 'Cluster Mutations Cellular Fraction'
    
    
            if(ylabel == ''):
                ylabel = 'Density'
    
    
            ax.set_xlabel(xlabel, fontsize = 17, labelpad = 15)
            ax.set_ylabel(ylabel, fontsize = 17, labelpad = 15)
        
            ax.set_title(ax.get_title(), fontsize = 18, y = 1.045)
            
            
            # default of axis is both
            ax.tick_params(labelsize = 15)
            
            
            # avoid scientific notation in y-axis
            ax.ticklabel_format(style = 'plain', axis = 'y')
        
        
        figure_name = "CCF mutations distribution per cluster with {n_c} clusters used while fitting and {n_r} random restarts of variational inference"  \
                         .format(n_c = miceExecSettings['num_clusters'].values[0], n_r = miceExecSettings['num_restarts'].values[0])
    
        plt.suptitle(figure_name, fontsize = 18, y = 1.13).set_weight("bold")
    
        
        group_legend = plt.gca().legend(handles = plt.gca().get_lines()[::-1], labels = ['Control', 'Sunitinib'], bbox_to_anchor=(1.08, 0.5), loc='center left', \
                                        borderaxespad=1.0, fontsize = 16, title = r'$\bf{Group}$', title_fontsize = 17, shadow = True, facecolor = 'white', borderpad = 0.6, \
                                        labelspacing = 0.6, handletextpad = 0.5)    
    
            
        for legHan in group_legend.legendHandles:
            legHan.set_linewidth(4)
    
    
        # Increase legend title pad from legend labels
            
        legendTitle = group_legend.get_title()
        
        legendTitle.set_y(5)
    
    
    
        # save figure
        figure_name = figure_name + '.png'
          
        plt.savefig(figure_name, bbox_inches = 'tight')    
        
        plt.show()
    



def allMiceMeanClusterCCFsPlot(allMiceTypeDF, outputDir):
    """
    Plot the average cellular_prevalence (or CCFs since purity = 1) distribution of each of the mice clusters, using all mice together.

    Parameters
            allMiceTypeDF (dataframe): all the different mice clusters average CCFs dataframe type (with or without considering CNAs errors).
            outputDir (str): the directory to where the plots will be output.
    """
    
    
    os.chdir(outputDir)
    
    # Converting the CCF values to percentage
    allMiceTypeDF['mean_cellular_prevalence'] = allMiceTypeDF['mean_cellular_prevalence'] * 100
    
    
    # Plot for each execution settings
    for execution in allMiceTypeDF['execution_id'].unique():
        
        miceExecSettings = allMiceTypeDF[allMiceTypeDF['execution_id'] == execution]
        
        
        # Plot the mice averge clusters CCFs data points colored by group type
        g = sns.catplot(x = "mouse", y = "mean_cellular_prevalence", hue = "group", palette=sns.color_palette("tab10"), \
                        col = "cluster_id", data = miceExecSettings, legend = False, sharex = False, sharey = False)
        
            
        g.set_axis_labels("Mouse ID", "Cancer Cell Fraction").set_titles ("Cluster {col_name}")
        
        
        # Define axes labels size
        for ax in g.axes.flat:
            
            xlabel = ax.get_xlabel()
            
            ylabel = ax.get_ylabel()
            
            
            if(xlabel == ''):
                xlabel = 'Mouse ID'
        
            if(ylabel == ''):
                ylabel = 'Cancer Cell Fraction'
        
            
            ax.set_xlabel(xlabel, fontsize = 17, labelpad = 15)
            ax.set_ylabel(ylabel, fontsize = 17, labelpad = 15)
        
            ax.set_title(ax.get_title(), fontsize = 18, y = 1.04)
            
            
            # default of axis is both
            ax.tick_params(labelsize = 15)
            
            
            # define ymin and ymax defined on the maximum point value 
            
            # get the maximum of all maximums of the 4 mice
            new_max = np.max([np.max(ax.collections[miceSamplePoint].get_offsets(), axis = 0)[1] for miceSamplePoint in range(0, 4)])
            
            new_min = np.min([np.min(ax.collections[miceSamplePoint].get_offsets(), axis = 0)[1] for miceSamplePoint in range(0, 4)])
        
        
            ax.set_ylim(new_min * 0.95, new_max * 1.05)
         
            
            
            for mp in ax.get_children()[0:4]:
        
                # increase mouse data points size
                mp.set_sizes([50,50,50,50])
            
        
            
        g.fig.tight_layout()   
        
        g.fig.subplots_adjust(wspace = 0.40)
            
        
        figure_name = "Average cluster CCF with {n_c} clusters used while fitting and {n_r} random restarts of variational inference" \
                         .format(n_c = miceExecSettings['num_clusters'].values[0], n_r = miceExecSettings['num_restarts'].values[0])
        
        plt.suptitle(figure_name, fontsize = 18, y = 1.12).set_weight("bold")
        
        
        group_legend = plt.gca().legend(bbox_to_anchor=(1.10, 0.5), loc='center left', borderaxespad=1.0, fontsize = 15, labels = ['Control', 'Sunitinib'], title = r'$\bf{Group}$', \
                                        title_fontsize = 16, shadow = True, facecolor = 'white', borderpad = 0.6, labelspacing = 0.6, handletextpad = 0.2)    
        
        
        # same color as data points
        group_legend.legendHandles[0].set_color(plt.gca().get_children()[1].get_facecolor()[0:3])
        group_legend.legendHandles[1].set_color(plt.gca().get_children()[2].get_facecolor()[0:3])
        
        
        # set markers' sizes
        group_legend.legendHandles[0].set_sizes([20])
        group_legend.legendHandles[1].set_sizes([20])
        
        
        for legHan in group_legend.legendHandles:
            legHan.set_linewidth(5.5)
        
            
        # Increase legend title pad from legend labels
            
        legendTitle = group_legend.get_title()
        
        legendTitle.set_y(4.3)
    
        
    
        # save figure
        figure_name = figure_name + '.png'
          
        plt.savefig(figure_name, bbox_inches = 'tight')    
        
        plt.show()
        
        
        

def allMiceClusterSampleHIPlot(allMiceTypeDF, outputDir):
    """
    Plot the clusters-based heterogeneity index for each of the mice, using all mice together.

    Parameters
            allMiceTypeDF (dataframe): all the different mice clusters heterogeneity index dataframe type (with or without considering CNAs errors).
            outputDir (str): the directory to where the plots will be output.
    """
    
    
    os.chdir(outputDir)
    
    
    # Checking if the data is normally distributed so we can use other statistical measures like t-test
    
    
    #####  Sample tumor heterogeneity index distribution by execution  #####
    
    
    g = sns.FacetGrid(allMiceTypeDF[['execution_id','heterogeneity_index']], col = 'execution_id', col_wrap = 3, aspect = 1.2, height = 5.6, sharex = False, sharey  = False)
    
    # Histogram plot to access the normality of data
    g.map(sns.histplot, "heterogeneity_index", kde = True, legend = False)
    
    
    
    g.fig.set_size_inches(24, 16)
        
    g.fig.subplots_adjust(hspace = 0.67)
      
    g.fig.subplots_adjust(wspace = 0.28) 
    
    
    
    g.set_xlabels('Sample TH index')
    
    titles = ["Execution with " + str(allMiceTypeDF[allMiceTypeDF['execution_id'] == execution]['num_clusters'].values[0]) + " clusters for fitting and " + \
              str(allMiceTypeDF[allMiceTypeDF['execution_id'] == execution]['num_restarts'].values[0]) + " random restarts" for execution in allMiceTypeDF['execution_id'].unique()]    
        
    
    # Define axes labels size
    for ax, title in zip(g.axes.flat, titles):
        
        ax.set_title(title, fontsize = 17.4, y = 1.10)
    
    
        xlabel = ax.get_xlabel()
        ylabel = ax.get_ylabel()
        
        
        if(xlabel == ''):
            xlabel = 'Sample TH index'
    
        if(ylabel == ''): 
            ylabel = 'Density'
    
         
        ax.set_xlabel(xlabel, fontsize = 16.8, labelpad = 11)
        ax.set_ylabel(ylabel, fontsize = 16.8, labelpad = 11)
    
        
        # default of axis is both
        ax.tick_params(labelsize = 16, pad = 13)
    
    
    
    figure_name = "Sample tumor heterogeneity index distribution by execution"
    
    plt.suptitle(figure_name, fontsize = 19, y = 1.062).set_weight("bold")
    
    
    
    # save figure
    figure_name = figure_name + '.png'
      
    plt.savefig(figure_name, bbox_inches = 'tight')    
    
    plt.show()
    
    
    
    
    
    
    
    #####  Quantile-Quantile plot of the sample tumor heterogeneity index and normal distributions  #####
    
    
    g = sns.FacetGrid(allMiceTypeDF[['execution_id','heterogeneity_index']], col = 'execution_id', col_wrap = 3, aspect = 1.2, height = 5.6, sharex = False, sharey  = False)
    
    
    titles = ["Execution with " + str(allMiceTypeDF[allMiceTypeDF['execution_id'] == execution]['num_clusters'].values[0]) + " clusters for fitting and " + \
              str(allMiceTypeDF[allMiceTypeDF['execution_id'] == execution]['num_restarts'].values[0]) + " random restarts" for execution in allMiceTypeDF['execution_id'].unique()]  
    
    
    for ax, execution in zip(g.axes.flat, allMiceTypeDF['execution_id'].unique()):
        
        g.map(sm.qqplot, data = allMiceTypeDF[ allMiceTypeDF['execution_id'] == execution] [['heterogeneity_index']].to_numpy(), line = 's', ax = ax, 
              fmt = 'none') # to avoid conflict between fmt and kwargs user warning: https://matplotlib.org/stable/_modules/matplotlib/axes/_base.html
             
    
    
    # Define axes labels size
    for ax, title in zip(g.axes.flat, titles):
            
        ax.set_title(title, fontsize = 17.4, y = 1.10)
    
    
    
        xlabel = ax.get_xlabel()
        ylabel = ax.get_ylabel()
        
        
        if(xlabel == ''):
            xlabel = 'Theoretical Quantiles'
    
        if(ylabel == ''): 
            ylabel = 'Sample Quantiles'
    
         
        ax.set_xlabel(xlabel, fontsize = 16.9, labelpad = 11)
        ax.set_ylabel(ylabel, fontsize = 16.9, labelpad = 11)
    
        
        # default of axis is both
        ax.tick_params(labelsize = 16, pad = 13)
    
    
    
    g.fig.set_size_inches(25, 16)
        
    g.fig.subplots_adjust(hspace = 0.675)
      
    g.fig.subplots_adjust(wspace = 0.285) 
    
    
    
    figure_name = "Quantile-Quantile plot of the sample tumor heterogeneity index and normal distributions"
    
    plt.suptitle(figure_name, fontsize = 19, y = 1.062).set_weight("bold")
    
    
    # save figure
    figure_name = figure_name + '.png'
      
    plt.savefig(figure_name, bbox_inches = 'tight')    
    
    plt.show()    
        
        
    
    
    
    

    #####   Wilcoxon rank-sum test (or Mann-Whitney U test)  #####
    
    
    # The Mann-Whitney U test is a nonparametric statistical significance test for determining whether two independent samples were drawn from a population with the same distribution.
    
    p_values = [0] * len(allMiceTypeDF['execution_id'].unique())
    
    
    # The different number of pyClone-VI executions is 12
    for execution, execID in zip(allMiceTypeDF['execution_id'].unique(), range(0, len(allMiceTypeDF['execution_id'].unique()))):
    
        sampleX = allMiceTypeDF[ (allMiceTypeDF['group'] == 'Ctrl') & (allMiceTypeDF['execution_id'] == execution)] [['execution_id', 'heterogeneity_index']].reset_index(drop = True) \
                    [['heterogeneity_index']].to_numpy()
        sampleY = allMiceTypeDF[ (allMiceTypeDF['group'] == 'Sunit') & (allMiceTypeDF['execution_id'] == execution)] [['execution_id', 'heterogeneity_index']].reset_index(drop = True) \
                    [['heterogeneity_index']].to_numpy()
    
    
        print("Execution with {num_clusters} clusters for fitting and {num_restarts} random restarts:".format(
                num_clusters = str(allMiceTypeDF[allMiceTypeDF['execution_id'] == execution]['num_clusters'].unique()[0]),
                num_restarts = str(allMiceTypeDF[allMiceTypeDF['execution_id'] == execution]['num_restarts'].unique()[0])))

    
        # compare samples
        stat, p = mannwhitneyu(sampleX, sampleY) # alternative = None (default)
        
        p_values[execID] = float(round(Decimal(p), 3))
        
        
        # The Mann-Whitney U statistic corresponding with sample x
        print('Statistics = %.3f, p = %.3f' % (stat, p))
        
        # interpretation
        alpha = 0.05
        if p > alpha:
            print('Same distribution (fail to reject H0)\n')
        else:
            print('Different distribution (reject H0)\n')         
        
        
        
        
    #####  Samples' heterogeneity index value based on the clusters CCF  #####
        
        
    # point is the mean value of heterogeneity index of all the samples of that mice for that execution
    
    g = sns.catplot(x = "mouse", y = "heterogeneity_index", hue = "group", palette=sns.color_palette("tab10"), \
                    col = "execution_id", data = allMiceTypeDF, legend = False, sharex = False, sharey = False, \
                    col_wrap= 3, kind = "point", join = False, height = 5.7, aspect = 1.8)
        
    
    g.fig.set_size_inches(38, 22.5)
        
    g.fig.subplots_adjust(hspace = 0.6)
    
    g.fig.subplots_adjust(wspace = 0.43)
        
    
    
    g.set_axis_labels("Mouse ID", "Sample TH Index")
    
    
    titles = ["Execution with " + str(allMiceTypeDF[allMiceTypeDF['execution_id'] == execution]['num_clusters'].values[0]) + " clusters for fitting and " + \
              str(allMiceTypeDF[allMiceTypeDF['execution_id'] == execution]['num_restarts'].values[0]) + " random restarts " for execution in allMiceTypeDF['execution_id'].unique()]    
        
        
        
    for ax, title, execID in zip(g.axes.flat, titles, range(0, len(allMiceTypeDF['execution_id'].unique()))):
        
        p_value = "(" + r'$\bf{{p = {p_value}}}$'.format(p_value = str(p_values[execID])) + ")" 
        
        
        ax.set_title(title + p_value, fontsize = 25,  y = 1.10)
    
       
        xlabel = ax.get_xlabel()
        
        ylabel = ax.get_ylabel()
        
        
        if(xlabel == ''):
            xlabel = 'Mouse ID'
    
        if(ylabel == ''): 
            ylabel = 'Sample TH Index'
    
                    
        ax.set_xlabel(xlabel, fontsize = 24, labelpad = 16.5)
        ax.set_ylabel(ylabel, fontsize = 24, labelpad = 16.5)
    
    
        # default of axis is both
        ax.tick_params(labelsize = 24)
    
        # control group mean points
        ax.get_children()[0].set_sizes([115])
    
        # sunitinib mean points
        ax.get_children()[1].set_sizes([115])
    
    
        # increase lines width
        for mp in ax.get_children()[2:10]:
            mp.set_linewidth(4.5)
    
    
    
    figure_name = "Comparison of the control and treated groups samples' heterogeneity index variation in different executions" 
    
    
    plt.suptitle(figure_name, fontsize = 28, y = 1.065).set_weight("bold")
    
    
    
    group_legend = plt.gca().legend(bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=1.0, fontsize = 24.5, \
                                    labels = ['Control', 'Sunitinib', '(' + r'$\bf{{\alpha = 0.05}}$' + ')'], title = r'$\bf{Group}$', title_fontsize = 26, shadow = True, \
                                    facecolor = 'white', borderpad = 0.8, labelspacing = 0.6, handletextpad = 0.7)    
    
    
    # same color as data points
    group_legend.legendHandles[0].set_color(plt.gca().get_children()[0].get_facecolor()[0:3])
    group_legend.legendHandles[1].set_color(plt.gca().get_children()[1].get_facecolor()[0:3])
    group_legend.legendHandles[2].set_visible(False)
    
    
    # set markers' sizes
    group_legend.legendHandles[0].set_markersize(23)
    group_legend.legendHandles[1].set_markersize(23)
    
    
    # the sunitinib bar was more thick
    for legHan in group_legend.legendHandles[ : -1]:
        legHan.set_linewidth(6)
    
        
    # Increase legend title pad from legend labels
        
    legendTitle = group_legend.get_title()
    
    legendTitle.set_y(6.5)
 
    
    # save figure
    figure_name = figure_name + '.png'
      
    plt.savefig(figure_name, bbox_inches = 'tight')    
    
    plt.show()
    



def difMiceNumberClusters(difMiceTypeDF, outputDir):
    """
    Plot the number of clusters obtained in each execution for each mouse separately.

    Parameters
            difMiceTypeDF (dataframe): the different mouse number of clusters obtained dataframe type (with or without considering CNAs errors).
            outputDir (str): the directory to where the plots will be output.
    """
    
    
    os.chdir(outputDir)
    
        
    # All mice with the same execution setting identifier to allow the different mouse files number of clusters obtained to be compared in a single subplot 
    
    difMiceTypeDF['execution_id'] = difMiceTypeDF['execution_id'].str[4:]
    
    
    # Bar plot comparing the number of clusters obtained after each PyClone-VI execution setting
    g = sns.catplot(x = "mouse", y = "output_clusters_obtained", hue = "group",  palette=sns.color_palette("tab10"), col = "execution_id", \
                    data = difMiceTypeDF, legend = False, kind = "bar", sharex = False, sharey = False, col_wrap= 3, height = 5.7, aspect = 1.8)
    
        
    g.fig.subplots_adjust(hspace = 0.65)
    
    g.fig.subplots_adjust(wspace = 0.4)
    
    
    g.set_axis_labels("Mouse ID", "Number of Clones")
    
    
    titles = ["Execution with " + str(difMiceTypeDF[difMiceTypeDF['execution_id'] == execution]['num_clusters'].values[0]) + " clusters for fitting and " + \
              str(difMiceTypeDF[difMiceTypeDF['execution_id'] == execution]['num_restarts'].values[0]) + " random restarts" for execution in difMiceTypeDF['execution_id'].unique()]
       
        
    # Define axes labels size
    for ax, title in zip(g.axes.flat, titles):
      
        ax.set_title(title, fontsize = 24.3, y = 1.12)    
    
    
        xlabel = ax.get_xlabel()
    
        ylabel = ax.get_ylabel()
    
    
        if(xlabel == ''):
            xlabel = 'Mouse ID'
    
        if(ylabel == ''): 
            ylabel = 'Number of Clones'
    
                
        ax.set_xlabel(xlabel, fontsize = 24, labelpad = 22)
        ax.set_ylabel(ylabel, fontsize = 24, labelpad = 22)
    
    
        # default of axis is both
        ax.tick_params(labelsize = 24)
         
       

        # Center data bars on ax ticks 
       
        n_bar = 0
       
    
        for bar_cord in range(0,8): 
        
            if (bar_cord == 4):
                n_bar = 0
        
        
            # gets left coordinate of each bar
            l_coord = str(ax.patches[bar_cord].get_x())
           
            # gets bar width to find right coordinate of each bar
            bar_width = str(ax.patches[bar_cord].get_width())
           
            
            coords = [l_coord, bar_width] 
           
        
            r_coord = sum(Decimal(i) for i in coords)
           
        
            # middle point of the data bar
            pos_shift = abs((r_coord - Decimal(l_coord)) / Decimal(2))
           
            coords[1] = pos_shift
           
            
            # get the databar half tick
            xlocs = ax.get_xticks()
        
            tick_coord = xlocs[n_bar]
         
            
            n_bar = n_bar + 1
        
                
            # if the bar is located to the left of the tick (the bar needs to move right)
            if ((Decimal(int(tick_coord)) + Decimal(0.001)).compare(abs(Decimal(r_coord))) == 1):
                # left coordinate + position shift
                ax.patches[bar_cord].set_x(float(sum(Decimal(i) for i in coords)))
            
            # if the bar is located to the right of the tick (the bar needs to move left)
            else:
                ax.patches[bar_cord].set_x(float(Decimal(l_coord) - pos_shift))
        
        
            # min and max number of clusters possible
            ax.set_yticks([0, 1, 2, 3, 4, 5])
    
        
      
    figure_name = "Comparison of the number of clones obtained for each mouse in each execution"
        
        
    plt.suptitle(figure_name, fontsize = 27, y = 1.065).set_weight("bold")
    
        
    group_legend = plt.gca().legend(bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=1.0, fontsize = 23.4, labels = ['Control', 'Sunitinib'], title =  r'$\bf{Group}$', \
                                    title_fontsize = 24.4, shadow = True, facecolor = 'white', borderpad = 0.8, labelspacing = 0.6, handletextpad = 0.7)    
    
    
    # same color as data points
    group_legend.legendHandles[0].set_color(plt.gca().get_children()[0].get_facecolor()[0:3])
    group_legend.legendHandles[1].set_color(plt.gca().get_children()[5].get_facecolor()[0:3])
    
    
    # set markers' sizes
    group_legend.legendHandles[0].set_markersize(20)
    group_legend.legendHandles[1].set_markersize(20)
    
    
        
    for legHan in group_legend.legendHandles:
        legHan.set_linewidth(6.5)
    
        
    # Increase legend title pad from legend labels
        
    legendTitle = group_legend.get_title()
    
    legendTitle.set_y(5)
    

    # save figure
    figure_name = figure_name + '.png'
      
    plt.savefig(figure_name, bbox_inches = 'tight')    
    
    plt.show()




def difMiceClusterSampleHIPlot(difMiceTypeDF, outputDir):
    """
    Plot the clusters-based heterogeneity index for each mouse separately.

    Parameters
            difMiceTypeDF (dataframe): the different mouse clusters heterogeneity index dataframe type (with or without considering CNAs errors).
            outputDir (str): the directory to where the plots will be output.
    """
    
    
    os.chdir(outputDir)
    
    
    # here we compare the samples' heterogeneity index using different mice files and not all mices in one file
    # even though the 4 mice are all in the x axis
    
    difMiceTypeDF['execution_id'] = difMiceTypeDF['execution_id'].str[4:]
    
    
    # point is the mean value of heterogeneity index of the samples of that mice for that execution
    
    g = sns.catplot(x = "mouse", y = "heterogeneity_index", hue = "group", palette=sns.color_palette("tab10"), \
                    col = "execution_id", data = difMiceTypeDF, legend = False, \
                    sharex = False, sharey = False, col_wrap = 3,  kind = "point", join = False, height = 5.7, aspect = 1.8)
        
    
    g.fig.subplots_adjust(hspace = 0.65)
    
    g.fig.subplots_adjust(wspace = 0.4)
        
    
    g.set_axis_labels("Mouse ID", "Sample TH Index")
    
    
    titles = ["Execution with " + str(difMiceTypeDF[difMiceTypeDF['execution_id'] == execution]['num_clusters'].values[0]) + " clusters for fitting and " + \
                        str(difMiceTypeDF[difMiceTypeDF['execution_id'] == execution]['num_restarts'].values[0]) + " random restarts" \
                                                  for execution in difMiceTypeDF['execution_id'].unique()]    
        
        
    # Define axes labels size
    for ax, title in zip(g.axes.flat, titles):
        
        ax.set_title(title, fontsize = 24.3,  y = 1.12)
    
       
        xlabel = ax.get_xlabel()
        
        ylabel = ax.get_ylabel()
        
        
        if(xlabel == ''):
            xlabel = 'Mouse ID'
    
        if(ylabel == ''):
            ylabel = 'Sample TH Index'
    
                    
        ax.set_xlabel(xlabel, fontsize = 24, labelpad = 22)
        ax.set_ylabel(ylabel, fontsize = 24, labelpad = 22)
    
    
        # default of axis is both
        ax.tick_params(labelsize = 24)
    
        # control group mean points
        ax.get_children()[0].set_sizes([115])
    
        # sunitinib mean points
        ax.get_children()[1].set_sizes([115])
    
    
        # increase lines width
        for mp in ax.get_children()[2:10]:
            mp.set_linewidth(4.5)
    
    
    
    figure_name = "Comparison of the control and treated groups samples' heterogeneity index variation in different executions" 
    
    
    plt.suptitle(figure_name, fontsize = 27, y = 1.065).set_weight("bold")
    
    
    group_legend = plt.gca().legend(bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=1.0, fontsize = 23.4, labels = ['Control', 'Sunitinib'], title =  r'$\bf{Group}$', \
                                    title_fontsize = 24.4, shadow = True, facecolor = 'white', borderpad = 0.8, labelspacing = 0.6, handletextpad = 0.7)    
    
    
    # same color as data points
    group_legend.legendHandles[0].set_color(plt.gca().get_children()[0].get_facecolor()[0:3])
    group_legend.legendHandles[1].set_color(plt.gca().get_children()[1].get_facecolor()[0:3])
    
    # set markers' sizes
    group_legend.legendHandles[0].set_markersize(20)
    group_legend.legendHandles[1].set_markersize(20)
    
        
    for legHan in group_legend.legendHandles:
        legHan.set_linewidth(5.5)
    
        
    # Increase legend title pad from legend labels
        
    legendTitle = group_legend.get_title()
    
    legendTitle.set_y(5)
    
    
    # save figure
    figure_name = figure_name + '.png'
      
    plt.savefig(figure_name, bbox_inches = 'tight')    
    
    plt.show()
    
    
    
    



#--------------------------------------------------------------------------------------------------




#####  Load PyClone-VI files results with all mice together files and each separate mouse files  #####


# Current directory: Dados/CNVS


os.chdir("./outputFilesPyCloneVI")


# Load all mice together and each separate mouse files including or not CNAs errors
miceSNVsInfo, miceSNVsInfoWithoutCNAsErrors = \
    [ 
         # considering or not the CNAs with major copy number 2 and minor copy number 0 as sequencing errors
         [ list(map(lambda filename: processOutputFiles(filename, withoutCNAsErrors), \
                sorted(glob.iglob(os.path.join(".", foldername, dirType, "") + "/*.tsv", recursive=False), key = orderFunc))) \
           
           for foldername in sorted(os.listdir(), key = str.lower)] \
     
                        for withoutCNAsErrors, dirType, orderFunc in zip([False, True], ["", "withoutCNAsErrors"], \
                                
                                    # first, for each filename get between o and r, with group 1 removing the included o
                                    # then split on c to get a list with the number of clusters and restarts ordering, converting them into ints
                                    [lambda x: list(map(int, re.search(r'o(.*?)r', x).group(1).split("c"))) , \
                                     
                                     # first get all occurences between s and r, then take the last one ([-1]) 
                                     # then split on c to get a list with the number of clusters and restarts ordering, converting them into ints
                                     lambda x: list(map(int, re.findall(r's(.*?)r', x, overlapped = True)[-1].split("c")))]) ]    
    
    
# Get the different lists of dataframes with the files proccessed
allMiceSNVsInfo, allMiceSNVsInfoWithoutCNAsErrors, differentMiceSNVsInfo, differentMiceSNVsInfoWithoutCNAsErrors = \
    [miceSNVsInfo[0], miceSNVsInfoWithoutCNAsErrors[0], miceSNVsInfo[1:], miceSNVsInfoWithoutCNAsErrors[1:]]


# remove variables to save space
del miceSNVsInfo, miceSNVsInfoWithoutCNAsErrors







#--------------------------------------------------------------------------------------------------




#####  Calculate data statistics for all mice together files and each separate mouse files, with and without considering CNAs errors  #####


# Obtain the cluster mean cellular prevalence and the sample heterogeneity index values dataframes for all mice together dataframes, with or without considering CNAs errors
allMiceStatistics, allMiceStatisticsWithoutCNAsErrors = \
    zip(*[list(map(lambda executionSettingsSample: miceStatistics(executionSettingsSample, False), allMiceDFs)) \
                          for allMiceDFs in zip(allMiceSNVsInfo, allMiceSNVsInfoWithoutCNAsErrors)])


# Unzip the returned pair of lists from the miceStatistics method, defining four variables, each of the two data statistics with or without considering CNAs 
(allMiceSNVsClusterCCFs, allMiceSNVsSampleHI), (allMiceSNVsClusterCCFsWithoutCNAsErrors, allMiceSNVsSampleHIWithoutCNAsErrors) \
    = [ map(lambda executionSettingsSample: pd.concat(executionSettingsSample, ignore_index = True), map(list, zip(*allMiceStatistics))), 
        map(lambda executionSettingsSample: pd.concat(executionSettingsSample, ignore_index = True), map(list, zip(*allMiceStatisticsWithoutCNAsErrors))) ]


# remove variables to save space
del allMiceStatistics, allMiceStatisticsWithoutCNAsErrors




# The number of clusters obtained in each execution and the sample heterogeneity index values dataframes for separated mouse files, with or without considering CNAs errors
sepMiceSNVs = [[list(map(lambda mouseExecSetting: miceStatistics(mouseExecSetting, True), mouse)) \
                    for mouse in differentMiceDFs] for differentMiceDFs in zip(differentMiceSNVsInfo, differentMiceSNVsInfoWithoutCNAsErrors)]


# Join all executions for each separate mouse statistics in 2 lists, with or without considering CNAs errors
sepMiceSNVs = [[list(map(list, zip(*differentMiceDFs))) for  differentMiceDFs in mouse] for mouse in sepMiceSNVs]


# Join all executions dataframes in one for each mice list position with both statistics calculated, with or without considering CNAs errors
sepMiceSNVs = [[list(map(lambda mouseStatisticsType: pd.concat(mouseStatisticsType, ignore_index = True), differentMiceDFs))
                    for differentMiceDFs in mouse] for mouse in sepMiceSNVs]


difMiceSNVsClusters, difMiceSNVsSampleHI, difMiceSNVsClustersWithoutCNAsErrors, difMiceSNVsSampleHIWithoutCNAsErrors = \
     [pd.concat([mouseStatistics[0][0] for mouseStatistics in sepMiceSNVs], ignore_index = True), 
      pd.concat([mouseStatistics[0][1] for mouseStatistics in sepMiceSNVs], ignore_index = True),
      pd.concat([mouseStatistics[1][0] for mouseStatistics in sepMiceSNVs], ignore_index = True),
      pd.concat([mouseStatistics[1][1] for mouseStatistics in sepMiceSNVs], ignore_index = True)] 


# remove variables to save space
del differentMiceSNVsInfo, differentMiceSNVsInfoWithoutCNAsErrors, sepMiceSNVs







#--------------------------------------------------------------------------------------------------




#####  Plot the PyClone-VI results statistics using all mice together executions (considering or not CNAs errors)  #####


# 1 - Using the mutations CCF values distribution of the mice clusters 

for allMiceDFsType, changeToDir in zip([allMiceSNVsInfo, allMiceSNVsInfoWithoutCNAsErrors], 
                                     ["../outputPyCloneViAnalysisImages/allMice/1 - clustersCCFsDistribution" , "./withoutCNAsErrors"]):

    allMiceClusterCCFsDistributionPlot(allMiceDFsType, changeToDir)
    
    
    

# 2 - Using the average CCF values of the mice clusters                          

for allMiceDFsType, changeToDir in zip([allMiceSNVsClusterCCFs, allMiceSNVsClusterCCFsWithoutCNAsErrors], 
                                     ["../../2 - clustersCCFsMean",  "./withoutCNAsErrors"]):

    allMiceMeanClusterCCFsPlot(allMiceDFsType, changeToDir)




# 3 - Using the sample heterogeneity index values based on the CCF of the mice clusters

for allMiceDFsType, changeToDir in zip([allMiceSNVsSampleHI, allMiceSNVsSampleHIWithoutCNAsErrors], 
                                     ["../../3 - clustersSampleHI",  "./withoutCNAsErrors"]):

    allMiceClusterSampleHIPlot(allMiceDFsType, changeToDir)







#####  Plot the PyClone-VI results statistics using different mouse executions (considering or not CNAs errors)  #####


# We can not use the clusters average CCF nor the mutations' CCF distribution of a cluster plots, 
# because the mice are in different files and we can not compare the clusters
# since a cluster in a mouse has no relation with a cluster in other mouse, 
# using the mice in different PyClone-VI executions
# this is different when using all mice together in a same file, when it is possible
    


# 1 - Comparing the number of clusters between each of the mice
for difMiceDFsType, changeToDir in zip([difMiceSNVsClusters, difMiceSNVsClustersWithoutCNAsErrors],
                                           ["../../../differentMice/1 - numClusters", "./withoutCNAsErrors"]):
    
    difMiceNumberClusters(difMiceDFsType, changeToDir)




# 2 - Using the sample heterogeneity index values based on the CCF of the clusters of each of the mice

for difMiceDFsType, changeToDir in zip([difMiceSNVsSampleHI, difMiceSNVsSampleHIWithoutCNAsErrors],
                                           ["../../2 - clustersSampleHI", "./withoutCNAsErrors"]):
    
    difMiceClusterSampleHIPlot(difMiceDFsType, changeToDir)






