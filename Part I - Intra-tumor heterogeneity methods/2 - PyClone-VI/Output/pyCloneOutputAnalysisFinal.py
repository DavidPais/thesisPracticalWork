#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 19:31:18 2021

@author: david
"""


#####  Libraries  #####


# Numpy for data management
import numpy as np

# Pandas also for data management
import pandas as pd

# Seaborn for plotting and styling
import seaborn as sns
sns.set_style("ticks")

# Matplotlib for additional customization
from matplotlib import pyplot as plt

# To open system files
import os

# For additional customization
import matplotlib as mat

# Reverse the matplotlib settings back to default (inclunding colors of graphs)  
mat.rc_file_defaults()

# To keep the precision in operations between floats
from decimal import Decimal

# For the QQ plot
import statsmodels.api as sm

# For the Wilcoxon rank-sum test (or Mann-Whitney U test)
from scipy.stats import mannwhitneyu

#--------------------------------------------------------------------------------------------------




#####  Methods  #####



# Calculate clusters' Shannon's Index

"""clusterCCF""" # cluster cell fraction (proportion of malignant cells in the cluster in the sample)

# returns: part of the shannon's index formula result
def shannonIndex(clusterCCF):
    return clusterCCF * np.log(clusterCCF)




#--------------------------------------------------------------------------------------------------





#####  Load pyClone files results with all mice and using all CNAs  #####


# Current directory: Dados/CNVS

os.chdir("./pyCloneOutputFiles/allMice/With all CNAs")


allMiceOutput = {}


# key to order files by number of clusters and then, if two files have the same number of clusters, by number of restarts
for filename in sorted(os.listdir(os.getcwd()), key = lambda x: list(map(int, x[ x.find( x[x.find("o")], x.find("o") + 1) + 1 : x.find("r")].split("c")))[0::]):
    
    # Not .h5 files, only the tsv files that resulted from the write_results pyClone algorithm method
    if(filename.endswith('.tsv')):
        
   
        with open(os.path.join(os.getcwd(), filename), 'r'):
        
             data = pd.read_csv(filename, sep='\t', header = 0)
             
           
             data['sample_id'] = [s[s.find('\'') + 1 : s.find( s[s.find('\'')] , s.find('\'') + 1)] for s in data['sample_id']] 
           
             data['mutation_id'] = [m[m.find('\'') + 1 : m.find( m[m.find('\'')] , m.find('\'') + 1)] for m in data['mutation_id']] 
           
            
             # select only those mutations which have more than 60% certainty of belonging to a cluster
             data = data [data['cluster_assignment_prob'] > 0.6]
            
             data['execution_id'] = filename[ : -4]
            
            
             data['mouse'] = [s[ : 3] for s in data['sample_id']]
            
             data['group'] = [s[3 : - 2] for s in data['sample_id']]   
            
             data['sample'] = [s[-2 : ] for s in data['sample_id']]
             
             
             
             data['num_clusters'] = [list(map(int, x[ x.find( x[x.find("o")], x.find("o") + 1) + 1 : x.find("r")].split("c")))[0] for x in data['execution_id']]
             
             data['num_restarts'] = [list(map(int, x[ x.find( x[x.find("o")], x.find("o") + 1) + 1 : x.find("r")].split("c")))[1] for x in data['execution_id']]
            
            
            
            
             data = data[['execution_id', 'mutation_id', 'sample_id', 'mouse', 'sample', 'group', 'num_clusters', 'num_restarts', 'cluster_id', 'cellular_prevalence', 'cellular_prevalence_std', \
                           'cluster_assignment_prob']]
            
            
                
             # drop rows with nan values
             data = data.dropna() 
       
             # reset indexes to start on 0 and increment by 1
             data = data.reset_index(drop = True)
        
            
            
             allMiceOutput[filename] = data
           



# (dictionary name).keys() to get the keys of the dictionary (in this case, the execution filenames)
allMiceDictKeys = sorted(list(allMiceOutput.keys()) , key = lambda x: list(map(int, x[ x.find( x[x.find("o")], x.find("o") + 1) + 1 : x.find("r")].split("c")))[0::])

# the sample_id set ma,es is the same for all execution run files of pyClone
samples = list(allMiceOutput[allMiceDictKeys[0]]['sample_id'].unique())


    



#--------------------------------------------------------------------------------------------------






#####  Load pyClone files results with all mice and not using all CNAs  #####



os.chdir("../Not all CNAs")


allMiceOutputNotAll = {}


# key to order files by number of clusters and then, if two files have the same number of clusters, by number of restarts
for filename in sorted(os.listdir(os.getcwd()), key = lambda x: list(map(int, x[ x.find("s") + 1 : x.find("r")].split("c")))[0::]):
    
    # Not .h5 files, only the tsv files that resulted from the write_results pyClone algorithm method
    if(filename.endswith('.tsv')):
        

        with open(os.path.join(os.getcwd(), filename), 'r'):
        
             data = pd.read_csv(filename, sep='\t', header = 0)
             
           
             data['sample_id'] = [s[s.find('\'') + 1 : s.find( s[s.find('\'')] , s.find('\'') + 1)] for s in data['sample_id']] 
           
             data['mutation_id'] = [m[m.find('\'') + 1 : m.find( m[m.find('\'')] , m.find('\'') + 1)] for m in data['mutation_id']] 
           
            
             # select only those mutations which have more than 60% certainty of belonging to a cluster
             data = data [data['cluster_assignment_prob'] > 0.6]
            
             data['execution_id'] = filename[ : -4]
            
            
             data['mouse'] = [ s[ : 3] for s in data['sample_id']]
            
             data['group'] = [ s[3 : - 2] for s in data['sample_id']]   
            
             data['sample'] = [ s[-2 : ] for s in data['sample_id']]
             

            
             data['num_clusters'] = [list(map(int, x[ x.find("s") + 1 : x.find("r")].split("c")))[0] for x in data['execution_id']]
             
             data['num_restarts'] = [list(map(int, x[ x.find("s") + 1 : x.find("r")].split("c")))[1] for x in data['execution_id']]
            
            
            

             data = data[['execution_id', 'mutation_id', 'sample_id', 'mouse', 'sample', 'group', 'num_clusters', 'num_restarts', 'cluster_id', 'cellular_prevalence', 'cellular_prevalence_std', \
                           'cluster_assignment_prob']]
            
            
            
             # drop rows with nan values
             data = data.dropna() 
       
             # reset indexes to start on 0 and increment by 1
             data = data.reset_index(drop = True)
        
            
            
             allMiceOutputNotAll[filename] = data
           



# (dictionary name).keys() to get the keys of the dictionary (in this case, the execution filenames)
allMiceNotAllDictKeys = sorted(list(allMiceOutputNotAll.keys()) , key = lambda x: list(map(int, x[ x.find("s") + 1 : x.find("r")].split("c")))[0::])






#--------------------------------------------------------------------------------------------------
    





#####  Load pyClone files results with different mice and using all CNAs  #####



os.chdir("../..")


differentMiceOutput = {}


# iterating each of the mice folders with the execution files
for foldername in sorted(os.listdir(), key = str.lower)[1:]:
    
    # last path component is empty, so a directory separator ('/') will be put at the end along with the concatenated value
    path = os.path.join(".", foldername, "With all CNAs", "")
    
    
    differentMiceOutput[foldername] = {}


    # key to order files by number of clusters and then, if two files have the same number of clusters, by number of restarts
    for filename in sorted(os.listdir(path), key = lambda x: list(map(int , x[x.find( x[ x.find("o")], x.find("o") + 1) + 1 : x.find("r")].split("c")[0::]))): 
        
        # Not .h5 files, only the tsv files that resulted from the write_results pyClone algorithm method
        if(filename.endswith('.tsv')):
        
            with open(os.path.join(path, filename), 'r') as file:
        
      
                data = pd.read_csv(file, sep='\t', header = 0)
  
  
                data['sample_id'] = [s[s.find('\'') + 1 : s.find( s[s.find('\'')] , s.find('\'') + 1)] for s in data['sample_id']] 
           
                data['mutation_id'] = [m[m.find('\'') + 1 : m.find( m[m.find('\'')] , m.find('\'') + 1)] for m in data['mutation_id']] 
        
            
                # select only those mutations which have more than 60% certainty of belonging to a cluster
                data = data [data['cluster_assignment_prob'] > 0.6]
            
                data['execution_id'] = filename[ : -4]
            
            
                data['mouse'] = [ s[ : 3] for s in data['sample_id']]
            
                data['group'] = [ s[3 : - 2] for s in data['sample_id']]   
            
                data['sample'] = [ s[-2 : ] for s in data['sample_id']]
             

                data['num_clusters'] = [list(map(int, x[ x.find( x[x.find("o")], x.find("o") + 1) + 1 : x.find("r")].split("c")))[0] for x in data['execution_id']]
             
                data['num_restarts'] = [list(map(int, x[ x.find( x[x.find("o")], x.find("o") + 1) + 1 : x.find("r")].split("c")))[1] for x in data['execution_id']]
            
            
            

                data = data[['execution_id', 'mutation_id', 'sample_id', 'mouse', 'sample', 'group', 'num_clusters', 'num_restarts', 'cluster_id', 'cellular_prevalence', 'cellular_prevalence_std', \
                           'cluster_assignment_prob']]
            
            
    
                # drop rows with nan values
                data = data.dropna() 
           
                # reset indexes to start on 0 and increment by 1
                data = data.reset_index(drop = True)
        
    
  
                differentMiceOutput[foldername][filename] = data
                
                
                        

        
        
#--------------------------------------------------------------------------------------------------






#####  Load pyClone files results with different mice and not using all CNAs  #####



differentMiceOutputNotAll = {}



# iterating each of the mice folders with the execution files
for foldername in sorted(os.listdir(), key = str.lower)[1:]:
    
    # last path component is empty, so a directory separator ('/') will be put at the end along with the concatenated value
    path = os.path.join(".", foldername, "Not all CNAs", "")
    
    
    differentMiceOutputNotAll[foldername] = {}
        

    # key to order files by number of clusters and then, if two files have the same number of 
    for filename in sorted(os.listdir(path), key = lambda x: list(map(int , x[ x.find("s") + 1 : x.find("r")].split("c")))): 
        
        # Not .h5 files, only the tsv files that resulted from the write_results pyClone algorithm method
        if(filename.endswith('.tsv')):
        
            with open(os.path.join(path, filename), 'r') as file:
        
               
                data = pd.read_csv(file, sep='\t', header = 0)
  
  
                data['sample_id'] = [s[s.find('\'') + 1 : s.find( s[s.find('\'')] , s.find('\'') + 1)] for s in data['sample_id']] 
           
                data['mutation_id'] = [m[m.find('\'') + 1 : m.find( m[m.find('\'')] , m.find('\'') + 1)] for m in data['mutation_id']] 
           
            
                # select only those mutations which have more than 60% certainty of belonging to a cluster
                data = data [data['cluster_assignment_prob'] > 0.6]
            
                data['execution_id'] = filename[ : -4]
            
            
                data['mouse'] = [ s[ : 3] for s in data['sample_id']]
            
                data['group'] = [ s[3 : - 2] for s in data['sample_id']]   
            
                data['sample'] = [ s[-2 : ] for s in data['sample_id']]
                         
             
                data['num_clusters'] = [list(map(int, x[ x.find("s") + 1 : x.find("r")].split("c")))[0] for x in data['execution_id']]
             
                data['num_restarts'] = [list(map(int, x[ x.find("s") + 1 : x.find("r")].split("c")))[1] for x in data['execution_id']]
            
            
            
            
                data = data[['execution_id', 'mutation_id', 'sample_id', 'mouse', 'sample', 'group', 'num_clusters', 'num_restarts', 'cluster_id', 'cellular_prevalence', 'cellular_prevalence_std', \
                           'cluster_assignment_prob']]
            
            
            
                # drop rows with nan values
                data = data.dropna() 
           
                # reset indexes to start on 0 and increment by 1
                data = data.reset_index(drop = True)
        
    
  
                differentMiceOutputNotAll[foldername][filename] = data
                




      
#--------------------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------------------






#####  Process files content with all mice and using all CNAs  #####



# Cellular_prevalence - proportion of malignant cells with the mutation in the sample. 
# This is also called cancer cell fraction (CCF) in the literature.               

# So this is the CCF of the cluster per sample


cluster_ccf = {}

sample_HI = {}



# iterating each execution file
for exec_name in allMiceDictKeys:

    h = 0

    # Here the cluster mean ccf is equal to the ccfs of all mutations since pyClone has the same cellular_prevalence for all the elements of a specific cluster
    cluster_ccf[exec_name] = allMiceOutput[exec_name].groupby(['execution_id', 'sample_id', 'group', 'mouse', 'sample', 'num_clusters', 'num_restarts', 'cluster_id'],\
                                                              as_index = False).agg(({'cellular_prevalence' : 'mean'}))


    l = {}

    exec_id = cluster_ccf[exec_name]['execution_id'].values[0]


    # for each of the samples
    for s in samples:
 
         sample_ccfs = list(cluster_ccf[exec_name][cluster_ccf[exec_name]['sample_id'] == s]['cellular_prevalence'])


         group = cluster_ccf[exec_name][cluster_ccf[exec_name]['sample_id'] == s]['group'].values[0]
         
         mouse = cluster_ccf[exec_name][cluster_ccf[exec_name]['sample_id'] == s]['mouse'].values[0]
         
         sample = cluster_ccf[exec_name][cluster_ccf[exec_name]['sample_id'] == s]['sample'].values[0]

         
         n_clusters = cluster_ccf[exec_name][cluster_ccf[exec_name]['sample_id'] == s]['num_clusters'].values[0]
         
         n_restarts = cluster_ccf[exec_name][cluster_ccf[exec_name]['sample_id'] == s]['num_restarts'].values[0]
         

         # heterogeneity index estimation for each sample's clusters
         heterogeneityIndex = -sum(list(map(shannonIndex, sample_ccfs)))
    
    
         l[h] = pd.Series(data = {'execution_id' : exec_id, 'sample_id' : s, 'mouse': mouse, 'sample': sample, 'num_clusters': n_clusters, 'num_restarts': n_restarts, \
                                  'group' : group, 'heterogeneity_index' : heterogeneityIndex })   
         
            
         h += 1
         
         
    
    sample_HI[exec_name] = pd.DataFrame.from_dict(l, "index")
    



# Join all dataframes in one for plotting

# comment the next lines if doing for execution in different dataframes (dictionary with execution keys and each entry is a dataframe of that execution)

cluster_ccf = pd.concat(list(cluster_ccf.values()), ignore_index = True)

sample_HI = pd.concat(list(sample_HI.values()), ignore_index = True)


# =============================================================================
# # For the distribution graph, comment the next line if separation by execution dicts is needed
allMiceOutput = pd.concat(list(allMiceOutput.values()), ignore_index = True)
# =============================================================================





#--------------------------------------------------------------------------------------------------





#####  Process files content with all mice and not using all CNAs  #####


clusterNotAll_ccf = {}

sampleNotAll_HI = {}


# iterating each execution file
for exec_name in allMiceNotAllDictKeys:

    h = 0

    # Here the cluster mean ccf is equal to the ccfs of all mutations since pyClone has the same cellular_prevalence for all the elements of a specific cluster
    clusterNotAll_ccf[exec_name] = allMiceOutputNotAll[exec_name].groupby(['execution_id', 'sample_id', 'group', 'mouse', 'sample', 'num_clusters', 'num_restarts', 'cluster_id'],\
                                                              as_index = False).agg(({'cellular_prevalence' : 'mean'}))


    l = {}

    exec_id = clusterNotAll_ccf[exec_name]['execution_id'].values[0]


    # for each of the samples
    for s in samples:
 
         sample_ccfs = list(clusterNotAll_ccf[exec_name][clusterNotAll_ccf[exec_name]['sample_id'] == s]['cellular_prevalence'])


         group = clusterNotAll_ccf[exec_name][clusterNotAll_ccf[exec_name]['sample_id'] == s]['group'].values[0]
         
         mouse = clusterNotAll_ccf[exec_name][clusterNotAll_ccf[exec_name]['sample_id'] == s]['mouse'].values[0]
         
         sample = clusterNotAll_ccf[exec_name][clusterNotAll_ccf[exec_name]['sample_id'] == s]['sample'].values[0]

         
         n_clusters = clusterNotAll_ccf[exec_name][clusterNotAll_ccf[exec_name]['sample_id'] == s]['num_clusters'].values[0]
         
         n_restarts = clusterNotAll_ccf[exec_name][clusterNotAll_ccf[exec_name]['sample_id'] == s]['num_restarts'].values[0]
         

         # heterogeneity index estimation for each sample's clusters
         heterogeneityIndex = -sum(list(map(shannonIndex, sample_ccfs)))
    
    
         l[h] = pd.Series(data = {'execution_id' : exec_id, 'sample_id' : s, 'mouse': mouse, 'sample': sample, 'num_clusters': n_clusters, 'num_restarts': n_restarts, \
                                  'group' : group, 'heterogeneity_index' : heterogeneityIndex })   
         
            
         h += 1
         
         
    
    sampleNotAll_HI[exec_name] = pd.DataFrame.from_dict(l, "index")
    



# Join all dataframes in one for plotting

# comment the next lines if doing for execution in different dataframes (dictionary with execution keys and each entry is a dataframe of that execution)

clusterNotAll_ccf = pd.concat(list(clusterNotAll_ccf.values()), ignore_index = True)

sampleNotAll_HI = pd.concat(list(sampleNotAll_HI.values()), ignore_index = True)


# =============================================================================
# # For the distribution graph, comment the next line if separation by execution dicts is needed
allMiceOutputNotAll = pd.concat(list(allMiceOutputNotAll.values()), ignore_index = True)
# =============================================================================





#--------------------------------------------------------------------------------------------------





#####  Process files content with different mice and using all CNAs  #####


mice_ccf = {}             

mice_sample_HI = {}


                
# iterating each mice
for m in list(differentMiceOutput.keys()):
    
    mice_ccf[m] = {}
    
    mice_sample_HI[m] = {}
    
    
    # iterating each execution file
    for exec_key in list(differentMiceOutput[m].keys()):
        
        mice_samples = list(differentMiceOutput[m][exec_key]['sample_id'].unique())

        
        l = {}
        
        h = 0
                

        mice_ccf[m][exec_key] = differentMiceOutput[m][exec_key].groupby(['execution_id', 'sample_id', 'group', 'mouse', 'sample', \
                                             'num_clusters', 'num_restarts', 'cluster_id'], as_index = False).agg(({'cellular_prevalence' : 'mean'}))
   
    
    
        exec_id = mice_ccf[m][exec_key]['execution_id'].values[0]
        
     
        # for each of the samples of the mice
        for ms in mice_samples:

            mice_sample_ccfs = list(mice_ccf[m][exec_key][ mice_ccf[m][exec_key]['sample_id'] == ms ] ['cellular_prevalence'])


            group = mice_ccf[m][exec_key][ mice_ccf[m][exec_key]['sample_id'] == ms ]['group'].values[0]
      
            mouse = mice_ccf[m][exec_key][ mice_ccf[m][exec_key]['sample_id'] == ms ]['mouse'].values[0]
         
            sample = mice_ccf[m][exec_key][ mice_ccf[m][exec_key]['sample_id'] == ms ]['sample'].values[0]
         
            
            n_clusters = mice_ccf[m][exec_key][ mice_ccf[m][exec_key]['sample_id'] == ms ]['num_clusters'].values[0]
         
            n_restarts = mice_ccf[m][exec_key][ mice_ccf[m][exec_key]['sample_id'] == ms ]['num_restarts'].values[0]  
      
        
        
            # heterogeneity index estimation for each sample's clusters
            heterogeneityIndex = -sum(list(map(shannonIndex, mice_sample_ccfs)))
            
    
            l[h] = pd.Series(data = {'execution_id' : exec_id, 'sample_id' : ms, 'mouse': mouse, 'sample': sample, 'num_clusters': n_clusters, 'num_restarts': n_restarts, \
                                     'group' : group, 'heterogeneity_index' : heterogeneityIndex })   
         
            
            h += 1
         
        

        # for a given mouse of the mice dictionary entries, every execution dictionary series rows are joined in a dataframe
        mice_sample_HI[m][exec_key] = pd.DataFrame.from_dict(l, "index")
    


    # uncomment both below if the dataframes must be separated by execution, an execution per dataframe instead of all executions together in one dataframe for plotting
    
    mice_ccf[m] = pd.concat(list(mice_ccf[m].values()), ignore_index = True)

    mice_sample_HI[m] = pd.concat(list(mice_sample_HI[m].values()), ignore_index = True)



# uncomment if the dataframes must be separated by mice_id 
# comment if using the first 2 graphs of comparing different mice (uncomment if not using)

mice_sample_HI = pd.concat(list(mice_sample_HI.values()), ignore_index = True)



# used for the process of obtaining the len of the clusters of each mouse in each execution (comparing the ctrl and sunit groups)

execs = pd.concat(list(mice_ccf.values()), ignore_index = True)

execs['execution_id'] = [x.split("_")[1] for x in execs['execution_id']]



# obtain number of clusters of different mice by sample in different executions
# for a same execution compared the number of clusters obtained by different mice

mice_clusters_info = {}


# execs is mice_ccf joining all the different mice values in one dataframe
for e in list(execs['execution_id'].unique()):
        
    mice_execution = execs[ execs['execution_id'] == e]


    l = {}
    
    h = 0


    # generating mice_clusters_info to get the number of clusters of each mouse in each execution
    for s in list(mice_execution['sample_id'].unique()):

        mice_nclusters = mice_execution[ mice_execution['sample_id'] == s]

        group = mice_execution[ mice_execution['sample_id'] == s]['group'].values[0]
         
        mouse = mice_execution[ mice_execution['sample_id'] == s]['mouse'].values[0]
         
        sample = mice_execution[ mice_execution['sample_id'] == s]['sample'].values[0]
         
        n_clusters = mice_execution[ mice_execution['sample_id'] == s]['num_clusters'].values[0]
         
        n_restarts = mice_execution[ mice_execution['sample_id'] == s]['num_restarts'].values[0]
     
        
        len_clusters = len(mice_nclusters['cluster_id'].values)

    

        l[h] = pd.Series(data = { 'execution_id' : e, 'sample_id' : s, 'mouse': mouse, 'sample': sample, 'num_clusters': n_clusters, 'num_restarts': n_restarts, \
                                          'group' : group, 'len_clusters' : len_clusters})


        h += 1


        mice_clusters_info[e] = pd.DataFrame.from_dict(l, "index")
        


mice_clusters_info = pd.concat(list(mice_clusters_info.values()), ignore_index = True)







#--------------------------------------------------------------------------------------------------






#####  Process files content with different mice and not using all CNAs  #####


miceNotAll_ccf = {}             

miceNotAll_sample_HI = {}


                
# iterating each mice
for m in list(differentMiceOutputNotAll.keys()):
    
    miceNotAll_ccf[m] = {}
    
    miceNotAll_sample_HI[m] = {}
    
    
    # iterating each execution file
    for exec_key in list(differentMiceOutputNotAll[m].keys()):
        
        mice_samples = list(differentMiceOutputNotAll[m][exec_key]['sample_id'].unique())

        
        l = {}
        
        h = 0
                

        miceNotAll_ccf[m][exec_key] = differentMiceOutputNotAll[m][exec_key].groupby(['execution_id', 'sample_id', 'group', 'mouse', 'sample', \
                                             'num_clusters', 'num_restarts', 'cluster_id'], as_index = False).agg(({'cellular_prevalence' : 'mean'}))
   
    
    
        exec_id = miceNotAll_ccf[m][exec_key]['execution_id'].values[0]
        
     
        # for each of the samples of the mice
        for ms in mice_samples:

            mice_sample_ccfs = list(miceNotAll_ccf[m][exec_key][ miceNotAll_ccf[m][exec_key]['sample_id'] == ms ] ['cellular_prevalence'])


            group = miceNotAll_ccf[m][exec_key][ miceNotAll_ccf[m][exec_key]['sample_id'] == ms ]['group'].values[0]
      
            mouse = miceNotAll_ccf[m][exec_key][ miceNotAll_ccf[m][exec_key]['sample_id'] == ms ]['mouse'].values[0]
         
            sample = miceNotAll_ccf[m][exec_key][ miceNotAll_ccf[m][exec_key]['sample_id'] == ms ]['sample'].values[0]
         
            
            n_clusters = miceNotAll_ccf[m][exec_key][ miceNotAll_ccf[m][exec_key]['sample_id'] == ms ]['num_clusters'].values[0]
         
            n_restarts = miceNotAll_ccf[m][exec_key][ miceNotAll_ccf[m][exec_key]['sample_id'] == ms ]['num_restarts'].values[0]  
      
        
        
            # heterogeneity index estimation for each sample's clusters
            heterogeneityIndex = -sum(list(map(shannonIndex, mice_sample_ccfs)))
            
    
            l[h] = pd.Series(data = {'execution_id' : exec_id, 'sample_id' : ms, 'mouse': mouse, 'sample': sample, 'num_clusters': n_clusters, 'num_restarts': n_restarts, \
                                     'group' : group, 'heterogeneity_index' : heterogeneityIndex })   
         
            
            h += 1
         
        

        # for a given mouse of the mice dictionary entries, every execution dictionary series rows are joined in a dataframe
        miceNotAll_sample_HI[m][exec_key] = pd.DataFrame.from_dict(l, "index")
    


    # uncomment both below if the dataframes must be separated by execution, an execution per dataframe instead of all executions together in one dataframe for plotting
    
    miceNotAll_ccf[m] = pd.concat(list(miceNotAll_ccf[m].values()), ignore_index = True)

    miceNotAll_sample_HI[m] = pd.concat(list(miceNotAll_sample_HI[m].values()), ignore_index = True)



# uncomment if the dataframes must be separated by mice_id 
# comment if using the first 2 graphs of comparing different mice (uncomment if not using)

miceNotAll_sample_HI = pd.concat(list(miceNotAll_sample_HI.values()), ignore_index = True)



# used for the process of obtaining the len of the clusters of each mouse in each execution (comparing the ctrl and sunit groups)

execsNotAll = pd.concat(list(miceNotAll_ccf.values()), ignore_index = True)

execsNotAll['execution_id'] = [x.split("_")[1] for x in execsNotAll['execution_id']]



# obtain number of clusters of different mice by sample in different executions
# for a same execution compared the number of clusters obtained by different mice

mice_clusters_infoNotAll = {}


# execs is mice_ccf joining all the different mice values in one dataframe
for e in list(execsNotAll['execution_id'].unique()):
        
    mice_execution = execsNotAll[ execsNotAll['execution_id'] == e]


    l = {}
    
    h = 0


    # generating mice_clusters_info to get the number of clusters of each mouse in each execution
    for s in list(mice_execution['sample_id'].unique()):

        mice_nclusters = mice_execution[ mice_execution['sample_id'] == s]

        group = mice_execution[ mice_execution['sample_id'] == s]['group'].values[0]
         
        mouse = mice_execution[ mice_execution['sample_id'] == s]['mouse'].values[0]
         
        sample = mice_execution[ mice_execution['sample_id'] == s]['sample'].values[0]
         
        n_clusters = mice_execution[ mice_execution['sample_id'] == s]['num_clusters'].values[0]
         
        n_restarts = mice_execution[ mice_execution['sample_id'] == s]['num_restarts'].values[0]
     
        
        len_clusters = len(mice_nclusters['cluster_id'].values)

    

        l[h] = pd.Series(data = { 'execution_id' : e, 'sample_id' : s, 'mouse': mouse, 'sample': sample, 'num_clusters': n_clusters, 'num_restarts': n_restarts, \
                                          'group' : group, 'len_clusters' : len_clusters})


        h += 1


        mice_clusters_infoNotAll[e] = pd.DataFrame.from_dict(l, "index")
        


mice_clusters_infoNotAll = pd.concat(list(mice_clusters_infoNotAll.values()), ignore_index = True)



#--------------------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------------------







#####  Visualize data with all mice and using all CNAs  #####


os.chdir("../pyCloneAnalysisImages/allMiceFile/With all CNAs")




# 1 - Using the clusters average ccf                            


os.chdir("./1 - Clusters average CCF")


# Converting the CCF to percentage
cluster_ccf['cellular_prevalence'] = cluster_ccf['cellular_prevalence'] * 100




for e in list(cluster_ccf['execution_id'].unique()):
    
    
    execution = cluster_ccf[ cluster_ccf['execution_id'] == e]
    

        
    g = sns.catplot(x = "mouse", y = "cellular_prevalence", hue = "group", palette=sns.color_palette("tab10"), \
      col = "cluster_id", data = execution, legend = False, \
          sharex = False, sharey = False)


    # refer on the thesis that cancer cell fraction is in percentage
    g.set_axis_labels("Mouse ID", "Cancer Cell Fraction").\
        set_titles ("Cluster {col_name}")
    
    
    
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
        
        
        
        new_max = np.max([np.max(ax.collections[i].get_offsets(), axis = 0)[1] for i in range(0, 4)])
        
        new_min = np.min([np.min(ax.collections[i].get_offsets(), axis = 0)[1] for i in range(0, 4)])
        
        

        ax.set_ylim(new_min * 0.95, new_max * 1.05)
     
        
        for mp in ax.get_children()[0:4]:

             # increase data points size
             mp.set_sizes([50,50,50,50])
        
  
    
    g.fig.tight_layout()   
    
    
    g.fig.subplots_adjust(hspace = 0.57)
    
    g.fig.subplots_adjust(wspace = 0.47)
        



    figure_name = "Average cluster CCF with {n_c} clusters used while fitting and {n_r} random restarts of variational inference" \
                     .format(n_c = execution['num_clusters'].values[0], n_r = execution['num_restarts'].values[0])
    
    plt.suptitle(figure_name, fontsize = 18, y = 1.12).set_weight("bold")

    

    group_legend = plt.gca().legend(bbox_to_anchor=(1.10, 0.5), loc='center left', borderaxespad=1.0, fontsize = 15, labels = ['Control', 'Sunitinib'], title = "Group", title_fontsize = 16, shadow = True, \
            facecolor = 'white', borderpad = 0.6, labelspacing = 0.6, handletextpad = 0.2)    
    
    
  
    # same color as data points
    group_legend.legendHandles[0].set_color(plt.gca().get_children()[1].get_facecolor()[0:3])
    group_legend.legendHandles[1].set_color(plt.gca().get_children()[2].get_facecolor()[0:3])
    
    
    # set markers' sizes
    group_legend.legendHandles[0].set_sizes([20])
    group_legend.legendHandles[1].set_sizes([20])


        
    for legHan in group_legend.legendHandles:
        legHan.set_linewidth(5.5)
    
        
    group_legend._legend_box.sep = 13
    


    # save figure
    figure_name = figure_name + '.png'
      
    plt.savefig(figure_name, bbox_inches = 'tight')    

    plt.show()






#--------------------------------------------------------------------------------------------------






# 2 - Using the samples' heterogeneity index value based on the clusters CCF

os.chdir("../2 - Heterogeneity index value")



# Checking if the data is normally distributed so we can use other statistical measures like t-test


#####  Sample tumor heterogeneity index distribution by execution  #####


g = sns.FacetGrid(sample_HI[['execution_id','heterogeneity_index']], col = 'execution_id', col_wrap = 3, aspect = 1.2, height = 5.6, 
                      sharex = False, sharey  = False)

g.map(sns.histplot, "heterogeneity_index", kde = True, legend = False)



g.fig.set_size_inches(24, 16)
    
g.fig.subplots_adjust(hspace = 0.67)
  
g.fig.subplots_adjust(wspace = 0.28) 




g.set_xlabels('Sample TH index')

titles = ["Execution with " + str(sample_HI[sample_HI['execution_id'] == e]['num_clusters'].values[0]) + " clusters for fitting and " + \
                    str(sample_HI[sample_HI['execution_id'] == e]['num_restarts'].values[0]) + " random restarts" \
                                              for e in list(sample_HI['execution_id'].unique())]    
    



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


g = sns.FacetGrid(sample_HI[['execution_id','heterogeneity_index']], col = 'execution_id', col_wrap = 3, aspect = 1.2, height = 5.6, 
                      sharex = False, sharey  = False)


titles = ["Execution with " + str(sample_HI[sample_HI['execution_id'] == e]['num_clusters'].values[0]) + " clusters for fitting and " + \
                    str(sample_HI[sample_HI['execution_id'] == e]['num_restarts'].values[0]) + " random restarts" \
                                              for e in list(sample_HI['execution_id'].unique())]  


    

for ax, execution in zip(g.axes.flat, sample_HI['execution_id'].unique()):
    
    g.map(sm.qqplot, data = sample_HI[ sample_HI['execution_id'] == execution] [['heterogeneity_index']].to_numpy(), line = 's', ax = ax, 
          # to avoid conflict between fmt and kwargs user warning: https://matplotlib.org/stable/_modules/matplotlib/axes/_base.html
         fmt = 'none')



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

p_values = [0] * len(sample_HI['execution_id'].unique())


for execution, i in zip(sample_HI['execution_id'].unique(), range(0,12)):

    sampleX = sample_HI[ (sample_HI['group'] == 'Ctrl') & (sample_HI['execution_id'] ==  execution)] [['execution_id', 'heterogeneity_index']].reset_index(drop = True) [['heterogeneity_index']].to_numpy()
    sampleY = sample_HI[ (sample_HI['group'] == 'Sunit') & (sample_HI['execution_id'] ==  execution)] [['execution_id', 'heterogeneity_index']].reset_index(drop = True) [['heterogeneity_index']].to_numpy()


    print("Execution with {num_clusters} clusters for fitting and {num_restarts} random restarts:".format(
            num_clusters = execution[ execution.find( execution[execution.find("o")], execution.find("o") + 1) + 1 : execution.find("r")].split("c")[0],
                num_restarts = execution[ execution.find( execution[execution.find("o")], execution.find("o") + 1) + 1 : execution.find("r")].split("c")[1]))


    # compare samples
    stat, p = mannwhitneyu(sampleX, sampleY) # alternative = None (default)
    
    p_values[i] = float(round(Decimal(p), 3))
    
    
    # The Mann-Whitney U statistic corresponding with sample x
    print('Statistics = %.3f, p = %.3f' % (stat, p))
    
    # interpretation
    alpha = 0.05
    if p > alpha:
    	print('Same distribution (fail to reject H0)\n')
    else:
    	print('Different distribution (reject H0)\n')         
    
    
    
    
    
    
    
    
    
#####  Samples' heterogeneity index value based on the clusters CCF  #####
    
    
# point is the mean value of heterogeneity index of the samples of that mice for that execution

g = sns.catplot(x = "mouse", y = "heterogeneity_index", hue = "group", palette=sns.color_palette("tab10"), \
      col = "execution_id", data = sample_HI, legend = False, \
          sharex = False, sharey = False, col_wrap= 3, kind = "point", join = False, height = 5.7, aspect = 1.8)
    

g.fig.set_size_inches(38, 22.5)
    
g.fig.subplots_adjust(hspace = 0.6)

g.fig.subplots_adjust(wspace = 0.43)
    


g.set_axis_labels("Mouse ID", "Sample TH Index")


titles = ["Execution with " + str(sample_HI[sample_HI['execution_id'] == e]['num_clusters'].values[0]) + " clusters for fitting and " + \
                    str(sample_HI[sample_HI['execution_id'] == e]['num_restarts'].values[0]) + " random restarts " \
                                              for e in list(sample_HI['execution_id'].unique())]    
    
    
    
for ax, title, p in zip(g.axes.flat, titles, range(0,12)):
    
    p_value = "(" + r'$\bf{{p = {p_value}}}$'.format(p_value = str(p_values[p])) + ")" 
    
    
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



group_legend = plt.gca().legend(bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=1.0, fontsize = 24.5, labels = ['Control', 'Sunitinib', '(' + r'$\bf{{\alpha = 0.05}}$' + ')'], title = "Group", title_fontsize = 26, shadow = True, \
        facecolor = 'white', borderpad = 0.8, labelspacing = 0.6, handletextpad = 0.7)    


# same color as data points
group_legend.legendHandles[0].set_color(plt.gca().get_children()[0].get_facecolor()[0:3])
group_legend.legendHandles[1].set_color(plt.gca().get_children()[1].get_facecolor()[0:3])
group_legend.legendHandles[2].set_visible(False)


# set markers' sizes
group_legend.legendHandles[0].set_markersize(23)
group_legend.legendHandles[1].set_markersize(23)


    
for legHan in group_legend.legendHandles[ : -1]:
    legHan.set_linewidth(6)

    
group_legend._legend_box.sep = 20


# save figure
figure_name = figure_name + '.png'
  
plt.savefig(figure_name, bbox_inches = 'tight')    

plt.show()

    
    
    
    
    
#--------------------------------------------------------------------------------------------------






# 3 - Using the mutations' CCF distribution of a cluster


os.chdir("../3 - Cluster mutations CCF distributions")


# Converting the CCF to percentage
allMiceOutput['cellular_prevalence'] = allMiceOutput['cellular_prevalence'] * 100



for e in list(allMiceOutput['execution_id'].unique()):
        
    
    execution = allMiceOutput[ allMiceOutput['execution_id'] == e]
    
    
    g = sns.displot(execution, x = "cellular_prevalence", hue = "group", col = "cluster_id", kind = "kde", cut = 0, \
                    palette=sns.color_palette("tab10", 2), legend = False, facet_kws={'sharey' : False, 'sharex' : False} ) 
  
    
        
    g.set_xlabels("Cluster Mutations Cellular Fraction").set_titles("Cluster {col_name}")  
        


    g.fig.tight_layout()   

    
    g.fig.subplots_adjust(hspace = 0.57)

    g.fig.subplots_adjust(wspace = 0.52)

    

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
        
    
    
    figure_name = "CCF mutations distribution per cluster with {n_c} clusters used while fitting and {n_r} random restarts of variational inference"  \
                     .format(n_c = execution['num_clusters'].values[0], n_r = execution['num_restarts'].values[0])


    plt.suptitle(figure_name, fontsize = 18, y = 1.13).set_weight("bold")


    
    group_legend = plt.gca().legend(handles = plt.gca().get_lines()[::-1], labels = ['Control', 'Sunitinib'], bbox_to_anchor=(1.08, 0.5), loc='center left', borderaxespad=1.0, fontsize = 16, title = "Group", title_fontsize = 17, shadow = True, \
            facecolor = 'white', borderpad = 0.6, labelspacing = 0.6, handletextpad = 0.5)    

        
    for legHan in group_legend.legendHandles:
        legHan.set_linewidth(4)

        
        
        
    group_legend._legend_box.sep = 15



    # save figure
    figure_name = figure_name + '.png'
      
    plt.savefig(figure_name, bbox_inches = 'tight')    
    
    plt.show()







#--------------------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------------------







#####  Visualize data with all mice and not using all CNAs  #####


os.chdir("../../Not all CNAs")



# 1 - Using the clusters average ccf                            


os.chdir("./1 - Clusters average CCF")


# Converting the CCF to percentage
clusterNotAll_ccf['cellular_prevalence'] = clusterNotAll_ccf['cellular_prevalence'] * 100




for e in list(clusterNotAll_ccf['execution_id'].unique()):
    
    
    execution = clusterNotAll_ccf[ clusterNotAll_ccf['execution_id'] == e]
    

        
    g = sns.catplot(x = "mouse", y = "cellular_prevalence", hue = "group", palette=sns.color_palette("tab10"), \
      col = "cluster_id", data = execution, legend = False, \
          sharex = False, sharey = False)


    # refer on the thesis that cancer cell fraction is in percentage
    g.set_axis_labels("Mouse ID", "Cancer Cell Fraction").\
        set_titles ("Cluster {col_name}")
    
    
    
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
        
        
        
        new_max = np.max([np.max(ax.collections[i].get_offsets(), axis = 0)[1] for i in range(0, 4)])
        
        new_min = np.min([np.min(ax.collections[i].get_offsets(), axis = 0)[1] for i in range(0, 4)])
        
        

        ax.set_ylim(new_min * 0.95, new_max * 1.05)
     
        
        for mp in ax.get_children()[0:4]:

             # increase data points size
             mp.set_sizes([50,50,50,50])
        
  
    
    g.fig.tight_layout()   
    
    
    g.fig.subplots_adjust(hspace = 0.57)
    
    g.fig.subplots_adjust(wspace = 0.47)
        



    figure_name = "Average cluster CCF with {n_c} clusters used while fitting and {n_r} random restarts of variational inference" \
                     .format(n_c = execution['num_clusters'].values[0], n_r = execution['num_restarts'].values[0])
    
    plt.suptitle(figure_name, fontsize = 18, y = 1.12).set_weight("bold")

    

    group_legend = plt.gca().legend(bbox_to_anchor=(1.10, 0.5), loc='center left', borderaxespad=1.0, fontsize = 15, labels = ['Control', 'Sunitinib'], title = "Group", title_fontsize = 16, shadow = True, \
            facecolor = 'white', borderpad = 0.6, labelspacing = 0.6, handletextpad = 0.2)    
    
    
  
    # same color as data points
    group_legend.legendHandles[0].set_color(plt.gca().get_children()[1].get_facecolor()[0:3])
    group_legend.legendHandles[1].set_color(plt.gca().get_children()[2].get_facecolor()[0:3])
    
    
    # set markers' sizes
    group_legend.legendHandles[0].set_sizes([20])
    group_legend.legendHandles[1].set_sizes([20])


        
    for legHan in group_legend.legendHandles:
        legHan.set_linewidth(5.5)
    
        
    group_legend._legend_box.sep = 13
    


    # save figure
    figure_name = figure_name + '.png'
      
    plt.savefig(figure_name, bbox_inches = 'tight')    

    plt.show()






#--------------------------------------------------------------------------------------------------






# 2 - Using the samples' heterogeneity index value based on the clusters CCF

os.chdir("../2 - Heterogeneity index value")



# Checking if the data is normally distributed so we can use other statistical measures like t-test


#####  Sample tumor heterogeneity index distribution by execution  #####


g = sns.FacetGrid(sampleNotAll_HI[['execution_id','heterogeneity_index']], col = 'execution_id', col_wrap = 3, aspect = 1.2, height = 5.6, 
                      sharex = False, sharey  = False)

g.map(sns.histplot, "heterogeneity_index", kde = True, legend = False)



g.fig.set_size_inches(24, 16)
    
g.fig.subplots_adjust(hspace = 0.67)
  
g.fig.subplots_adjust(wspace = 0.28) 




g.set_xlabels('Sample TH index')

titles = ["Execution with " + str(sampleNotAll_HI[sampleNotAll_HI['execution_id'] == e]['num_clusters'].values[0]) + " clusters for fitting and " + \
                    str(sampleNotAll_HI[sampleNotAll_HI['execution_id'] == e]['num_restarts'].values[0]) + " random restarts" \
                                              for e in list(sampleNotAll_HI['execution_id'].unique())]    
    



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


g = sns.FacetGrid(sampleNotAll_HI[['execution_id','heterogeneity_index']], col = 'execution_id', col_wrap = 3, aspect = 1.2, height = 5.6, 
                      sharex = False, sharey  = False)


titles = ["Execution with " + str(sampleNotAll_HI[sampleNotAll_HI['execution_id'] == e]['num_clusters'].values[0]) + " clusters for fitting and " + \
                    str(sampleNotAll_HI[sampleNotAll_HI['execution_id'] == e]['num_restarts'].values[0]) + " random restarts" \
                                              for e in list(sampleNotAll_HI['execution_id'].unique())]  


    

for ax, execution in zip(g.axes.flat, sampleNotAll_HI['execution_id'].unique()):
    
    g.map(sm.qqplot, data = sampleNotAll_HI[ sampleNotAll_HI['execution_id'] == execution] [['heterogeneity_index']].to_numpy(), line = 's', ax = ax, 
          # to avoid conflict between fmt and kwargs user warning: https://matplotlib.org/stable/_modules/matplotlib/axes/_base.html
         fmt = 'none')



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

p_values = [0] * len(sampleNotAll_HI['execution_id'].unique())


for execution, i in zip(sampleNotAll_HI['execution_id'].unique(), range(0,12)):

    sampleX = sampleNotAll_HI[ (sampleNotAll_HI['group'] == 'Ctrl') & (sampleNotAll_HI['execution_id'] ==  execution)] [['execution_id', 'heterogeneity_index']].reset_index(drop = True) [['heterogeneity_index']].to_numpy()
    sampleY = sampleNotAll_HI[ (sampleNotAll_HI['group'] == 'Sunit') & (sampleNotAll_HI['execution_id'] ==  execution)] [['execution_id', 'heterogeneity_index']].reset_index(drop = True) [['heterogeneity_index']].to_numpy()


    print("Execution with {num_clusters} clusters for fitting and {num_restarts} random restarts:".format(
            num_clusters = execution[ execution.find("s") + 1 : execution.find("r")].split("c")[0],
                num_restarts = execution[ execution.find("s") + 1 : execution.find("r")].split("c")[1]))


    # compare samples
    stat, p = mannwhitneyu(sampleX, sampleY) # alternative = None (default)
    
    p_values[i] = float(round(Decimal(p), 3))
    
    
    # The Mann-Whitney U statistic corresponding with sample x
    print('Statistics = %.3f, p = %.3f' % (stat, p))
    
    # interpretation
    alpha = 0.05
    if p > alpha:
    	print('Same distribution (fail to reject H0)\n')
    else:
    	print('Different distribution (reject H0)\n')         
    
    
    
    
    
    
    
    
    
#####  Samples' heterogeneity index value based on the clusters CCF  #####
    
    
# point is the mean value of heterogeneity index of the samples of that mice for that execution

g = sns.catplot(x = "mouse", y = "heterogeneity_index", hue = "group", palette=sns.color_palette("tab10"), \
      col = "execution_id", data = sampleNotAll_HI, legend = False, \
          sharex = False, sharey = False, col_wrap= 3, kind = "point", join = False, height = 5.7, aspect = 1.8)
    

g.fig.set_size_inches(38, 22.5)
    
g.fig.subplots_adjust(hspace = 0.6)

g.fig.subplots_adjust(wspace = 0.43)
    


g.set_axis_labels("Mouse ID", "Sample TH Index")


titles = ["Execution with " + str(sampleNotAll_HI[sampleNotAll_HI['execution_id'] == e]['num_clusters'].values[0]) + " clusters for fitting and " + \
                    str(sampleNotAll_HI[sampleNotAll_HI['execution_id'] == e]['num_restarts'].values[0]) + " random restarts " \
                                              for e in list(sampleNotAll_HI['execution_id'].unique())]    
    
    
    
for ax, title, p in zip(g.axes.flat, titles, range(0,12)):
    
    p_value = "(" + r'$\bf{{p = {p_value}}}$'.format(p_value = str(p_values[p])) + ")" 
    
    
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



group_legend = plt.gca().legend(bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=1.0, fontsize = 24.5, labels = ['Control', 'Sunitinib', '(' + r'$\bf{{\alpha = 0.05}}$' + ')'], title = "Group", title_fontsize = 26, shadow = True, \
        facecolor = 'white', borderpad = 0.8, labelspacing = 0.6, handletextpad = 0.7)    


# same color as data points
group_legend.legendHandles[0].set_color(plt.gca().get_children()[0].get_facecolor()[0:3])
group_legend.legendHandles[1].set_color(plt.gca().get_children()[1].get_facecolor()[0:3])
group_legend.legendHandles[2].set_visible(False)


# set markers' sizes
group_legend.legendHandles[0].set_markersize(23)
group_legend.legendHandles[1].set_markersize(23)


    
for legHan in group_legend.legendHandles[ : -1]:
    legHan.set_linewidth(6)

    
group_legend._legend_box.sep = 20


# save figure
figure_name = figure_name + '.png'
  
plt.savefig(figure_name, bbox_inches = 'tight')    

plt.show()






#--------------------------------------------------------------------------------------------------






# 3 - Using the mutations' CCF distribution of a cluster


os.chdir("../3 - Cluster mutations CCF distributions")


# Converting the CCF to percentage
allMiceOutputNotAll['cellular_prevalence'] = allMiceOutputNotAll['cellular_prevalence'] * 100



for e in list(allMiceOutputNotAll['execution_id'].unique()):
        
    
    execution = allMiceOutputNotAll[ allMiceOutputNotAll['execution_id'] == e]
    
    
    g = sns.displot(execution, x = "cellular_prevalence", hue = "group", col = "cluster_id", kind = "kde", cut = 0, \
                    palette=sns.color_palette("tab10", 2), legend = False, facet_kws={'sharey' : False, 'sharex' : False} ) 
  
    
        
    g.set_xlabels("Cluster Mutations Cellular Fraction").set_titles("Cluster {col_name}")  
        


    g.fig.tight_layout()   

    
    g.fig.subplots_adjust(hspace = 0.57)

    g.fig.subplots_adjust(wspace = 0.52)

    

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
        
    
    
    figure_name = "CCF mutations distribution per cluster with {n_c} clusters used while fitting and {n_r} random restarts of variational inference"  \
                     .format(n_c = execution['num_clusters'].values[0], n_r = execution['num_restarts'].values[0])


    plt.suptitle(figure_name, fontsize = 18, y = 1.13).set_weight("bold")


    
    group_legend = plt.gca().legend(handles = plt.gca().get_lines()[::-1], labels = ['Control', 'Sunitinib'], bbox_to_anchor=(1.08, 0.5), loc='center left', borderaxespad=1.0, fontsize = 16, title = "Group", title_fontsize = 17, shadow = True, \
            facecolor = 'white', borderpad = 0.6, labelspacing = 0.6, handletextpad = 0.5)    

        
    for legHan in group_legend.legendHandles:
        legHan.set_linewidth(4)

        
        
        
    group_legend._legend_box.sep = 15



    # save figure
    figure_name = figure_name + '.png'
      
    plt.savefig(figure_name, bbox_inches = 'tight')    
    
    plt.show()







#--------------------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------------------







#####  Visualize data with different mice and using all CNAs  #####


os.chdir("../../../differentMiceFiles/With all CNAs")



# 1 - Using the samples' heterogeneity index value based on the clusters CCF 


os.chdir("./1 - Heterogeneity index value")


# here we compare the samples' heterogeneity index using different mice files and not all mices in one file
# even though the 4 mice are all in the x axis

# apply to each of the elements of the column execution_id the lambda function
mice_sample_HI['execution_id'] = mice_sample_HI['execution_id'].apply(lambda e: e[4:])




# point is the mean value of heterogeneity index of the samples of that mice for that execution

g = sns.catplot(x = "mouse", y = "heterogeneity_index", hue = "group", palette=sns.color_palette("tab10"), \
      col = "execution_id", data = mice_sample_HI, legend = False, \
          sharex = False, sharey = False, col_wrap = 3,  kind = "point", join = False, height = 5.7, aspect = 1.8)
    


g.fig.subplots_adjust(hspace = 0.65)

g.fig.subplots_adjust(wspace = 0.4)
    


g.set_axis_labels("Mouse ID", "Sample TH Index")


titles = ["Execution with " + str(mice_sample_HI[mice_sample_HI['execution_id'] == e]['num_clusters'].values[0]) + " clusters for fitting and " + \
                    str(mice_sample_HI[mice_sample_HI['execution_id'] == e]['num_restarts'].values[0]) + " random restarts" \
                                              for e in list(mice_sample_HI['execution_id'].unique())]    
    
    

for ax, title in zip(g.axes.flat, titles):
    
    ax.set_title(title, fontsize = 24.3,  y = 1.12)

   
    xlabel = ax.get_xlabel()
    
    ylabel = ax.get_ylabel()
    
    
    if(xlabel == ''):
        xlabel = 'Mouse ID'

    if(ylabel == ''): # or Sample Heterogeneity Index variation if using kind = "point"
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


group_legend = plt.gca().legend(bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=1.0, fontsize = 23.4, labels = ['Control', 'Sunitinib'], title = "Group", title_fontsize = 24.4, shadow = True, \
        facecolor = 'white', borderpad = 0.8, labelspacing = 0.6, handletextpad = 0.7)    


# same color as data points
group_legend.legendHandles[0].set_color(plt.gca().get_children()[0].get_facecolor()[0:3])
group_legend.legendHandles[1].set_color(plt.gca().get_children()[1].get_facecolor()[0:3])

# set markers' sizes
group_legend.legendHandles[0].set_markersize(20)
group_legend.legendHandles[1].set_markersize(20)

    
for legHan in group_legend.legendHandles:
    legHan.set_linewidth(5.5)

    
group_legend._legend_box.sep = 15


# save figure
figure_name = figure_name + '.png'
  
plt.savefig(figure_name, bbox_inches = 'tight')    

plt.show()



    



#--------------------------------------------------------------------------------------------------







# We can not use the clusters average CCF nor the mutations' CCF distribution of a cluster plots, 
# because the mice are in different files and we can not compare the clusters distributions 
# since a cluster in a mouse has no relation with a cluster in other mouse


# 2 - Comparing the number of clusters between mice


os.chdir("../2 - Comparing the number of clusters")



g = sns.catplot(x = "mouse", y = "len_clusters", hue = "group",  \
          palette=sns.color_palette("tab10"), col = "execution_id", data = mice_clusters_info, legend = False, \
             kind = "bar", sharex = False, sharey = False, col_wrap= 3, height = 5.7, aspect = 1.8)

    
    
g.fig.subplots_adjust(hspace = 0.65)

g.fig.subplots_adjust(wspace = 0.4)



g.set_axis_labels("Mouse ID", "Number of Clones")



titles = ["Execution with " + str(mice_clusters_info[mice_clusters_info['execution_id'] == e]['num_clusters'].values[0]) + " clusters for fitting and " + \
                   str(mice_clusters_info[mice_clusters_info['execution_id'] == e]['num_restarts'].values[0]) + " random restarts" \
                                             for e in list(mice_clusters_info['execution_id'].unique())]
   
    
    
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
   
    k = 0
   

    for s in range(0,8): 
    
        if (s == 4):
            k = 0
    
    
        # gets left coordinate of each bar
        l_coord = str(ax.patches[s].get_x())
       
        w_h = str(ax.patches[s].get_width())
       
        coords = [l_coord, w_h] 
       
    
        r_coord = sum(Decimal(i) for i in coords)
       
    
        pos_shift = abs((r_coord - Decimal(l_coord)) / Decimal(2))
       
        coords[1] = pos_shift
       
    
       
        xlocs = ax.get_xticks()
    
        tick_coord = xlocs[k]
     
        
        k = k + 1
    
            
    
        if ((Decimal(int(tick_coord)) + Decimal(0.001)).compare(abs(Decimal(r_coord))) == 1):
            ax.patches[s].set_x(float(sum(Decimal(i) for i in coords)))
        
    
        else:
            ax.patches[s].set_x(float(Decimal(l_coord) - pos_shift ))
    
    
    
        ax.set_yticks([0, 1, 2, 3, 4, 5])

    
  
figure_name = "Comparison of the number of clones obtained for each mouse in each execution"
    
    
plt.suptitle(figure_name, fontsize = 27, y = 1.065).set_weight("bold")

    
group_legend = plt.gca().legend(bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=1.0, fontsize = 23.4, labels = ['Control', 'Sunitinib'], title = "Group", title_fontsize = 24.4, shadow = True, \
        facecolor = 'white', borderpad = 0.8, labelspacing = 0.6, handletextpad = 0.7)    


# same color as data points
group_legend.legendHandles[0].set_color(plt.gca().get_children()[0].get_facecolor()[0:3])
group_legend.legendHandles[1].set_color(plt.gca().get_children()[5].get_facecolor()[0:3])


# set markers' sizes
group_legend.legendHandles[0].set_markersize(20)
group_legend.legendHandles[1].set_markersize(20)


    
for legHan in group_legend.legendHandles:
    legHan.set_linewidth(6.5)

    
group_legend._legend_box.sep = 15



# save figure
figure_name = figure_name + '.png'
  
plt.savefig(figure_name, bbox_inches = 'tight')    

plt.show()







#--------------------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------------------







#####  Visualize data with different mice and not using all CNAs  #####


os.chdir("../../Not all CNAs")



# 1 - Using the samples' heterogeneity index value based on the clusters CCF 


os.chdir("./1 - Heterogeneity index value")


# here we compare the samples' heterogeneity index using different mice files and not all mices in one file
# even though the 4 mice are all in the x axis

# apply to each of the elements of the column execution_id the lambda function
miceNotAll_sample_HI['execution_id'] = miceNotAll_sample_HI['execution_id'].apply(lambda e: e[4:])




# point is the mean value of heterogeneity index of the samples of that mice for that execution

g = sns.catplot(x = "mouse", y = "heterogeneity_index", hue = "group", palette=sns.color_palette("tab10"), \
      col = "execution_id", data = miceNotAll_sample_HI, legend = False, \
          sharex = False, sharey = False, col_wrap = 3,  kind = "point", join = False, height = 5.7, aspect = 1.8)



g.fig.subplots_adjust(hspace = 0.65)

g.fig.subplots_adjust(wspace = 0.4)
    


g.set_axis_labels("Mouse ID", "Sample TH Index")


titles = ["Execution with " + str(miceNotAll_sample_HI[miceNotAll_sample_HI['execution_id'] == e]['num_clusters'].values[0]) + " clusters for fitting and " + \
                    str(miceNotAll_sample_HI[miceNotAll_sample_HI['execution_id'] == e]['num_restarts'].values[0]) + " random restarts" \
                                              for e in list(miceNotAll_sample_HI['execution_id'].unique())]    
    
    

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


group_legend = plt.gca().legend(bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=1.0, fontsize = 23.4, labels = ['Control', 'Sunitinib'], title = "Group", title_fontsize = 24.4, shadow = True, \
        facecolor = 'white', borderpad = 0.8, labelspacing = 0.6, handletextpad = 0.7)    


# same color as data points
group_legend.legendHandles[0].set_color(plt.gca().get_children()[0].get_facecolor()[0:3])
group_legend.legendHandles[1].set_color(plt.gca().get_children()[1].get_facecolor()[0:3])

# set markers' sizes
group_legend.legendHandles[0].set_markersize(20)
group_legend.legendHandles[1].set_markersize(20)

    
for legHan in group_legend.legendHandles:
    legHan.set_linewidth(5.5)

    
group_legend._legend_box.sep = 15


# save figure
figure_name = figure_name + '.png'
  
plt.savefig(figure_name, bbox_inches = 'tight')    

plt.show()







#--------------------------------------------------------------------------------------------------







# We can not use the clusters average CCF nor the mutations' CCF distribution of a cluster plots, 
# because the mice are in different files and we can not compare the clusters distributions 
# since a cluster in a mouse has no relation with a cluster in other mouse


# 2 - Comparing the number of clusters between mice


os.chdir("../2 - Comparing the number of clusters")



g = sns.catplot(x = "mouse", y = "len_clusters", hue = "group",  \
          palette=sns.color_palette("tab10"), col = "execution_id", data = mice_clusters_infoNotAll, legend = False, \
             kind = "bar", sharex = False, sharey = False, col_wrap= 3, height = 5.7, aspect = 1.8)

    
    
g.fig.subplots_adjust(hspace = 0.65)

g.fig.subplots_adjust(wspace = 0.4)



g.set_axis_labels("Mouse ID", "Number of Clones")



titles = ["Execution with " + str(mice_clusters_infoNotAll[ mice_clusters_infoNotAll['execution_id'] == e]['num_clusters'].values[0]) + " clusters for fitting and " + \
                   str(mice_clusters_infoNotAll[ mice_clusters_infoNotAll['execution_id'] == e]['num_restarts'].values[0]) + " random restarts" \
                                             for e in list(mice_clusters_infoNotAll['execution_id'].unique())]
   
    
    
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
   
    k = 0
   

    for s in range(0,8): 
    
        if (s == 4):
            k = 0
    
    
        # gets left coordinate of each bar
        l_coord = str(ax.patches[s].get_x())
       
        w_h = str(ax.patches[s].get_width())
       
        coords = [l_coord, w_h] 
       
    
        r_coord = sum(Decimal(i) for i in coords)
       
    
        pos_shift = abs((r_coord - Decimal(l_coord)) / Decimal(2))
       
        coords[1] = pos_shift
       
    
       
        xlocs = ax.get_xticks()
    
        tick_coord = xlocs[k]
     
        
        k = k + 1
    
            
    
        if ((Decimal(int(tick_coord)) + Decimal(0.001)).compare(abs(Decimal(r_coord))) == 1):
            ax.patches[s].set_x(float(sum(Decimal(i) for i in coords)))
        
    
        else:
            ax.patches[s].set_x(float(Decimal(l_coord) - pos_shift ))
    
    
    
        ax.set_yticks([0, 1, 2, 3, 4, 5])

    
  
figure_name = "Comparison of the number of clones obtained for each mouse in each execution"
    
    
plt.suptitle(figure_name, fontsize = 27, y = 1.065).set_weight("bold")

    
group_legend = plt.gca().legend(bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=1.0, fontsize = 23.4, labels = ['Control', 'Sunitinib'], title = "Group", title_fontsize = 24.4, shadow = True, \
        facecolor = 'white', borderpad = 0.8, labelspacing = 0.6, handletextpad = 0.7)    


# same color as data points
group_legend.legendHandles[0].set_color(plt.gca().get_children()[0].get_facecolor()[0:3])
group_legend.legendHandles[1].set_color(plt.gca().get_children()[5].get_facecolor()[0:3])


# set markers' sizes
group_legend.legendHandles[0].set_markersize(20)
group_legend.legendHandles[1].set_markersize(20)


    
for legHan in group_legend.legendHandles:
    legHan.set_linewidth(6.5)

    
group_legend._legend_box.sep = 15



# save figure
figure_name = figure_name + '.png'
  
plt.savefig(figure_name, bbox_inches = 'tight')    

plt.show()






