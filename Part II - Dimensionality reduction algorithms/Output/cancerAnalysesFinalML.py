# -*- coding: utf-8 -*-

"""
Created on Wed Sep 15 19:32:05 2021

@author: david
"""


# Use machine learning methods to assess the difference between the mice case study control and sunitinib groups
# To analyse the possible effects of the drug treatment on genes that may distinguish between both groups


#####  Libraries  #####


# Pandas for data management
import pandas as pd

# To access the needed dataframe objects with the input data
import os

# Matplotlib for additional customization
from matplotlib import pyplot as plt

# Reverse the matplotlib settings back to default (inclunding colors of graphs)  
import matplotlib as mat
mat.rc_file_defaults()

# Seaborn for plotting and styling
import seaborn as sns
sns.set_style("ticks")

# Principal Components Analyses
from sklearn.decomposition import PCA

# Truncated Singular Values Decomposition
from sklearn.decomposition import TruncatedSVD

# T-distributed Stochastic Neighbor Embedding
from sklearn.manifold import TSNE

# Multidimensional scaling
from sklearn.manifold import MDS

# Isomap Embedding (non-linear dimensionality reduction through Isometric Mapping)
from sklearn.manifold import Isomap

# Numpy for data management
import numpy as np

# Stats library to access the kernel density method
from scipy import stats


# Uncomment the two lines below when not using google colab

# Allow google colab file sharing
from google.colab import drive

# Connect to the drive folder root
drive.mount('/content/drive')




#--------------------------------------------------------------------------------------------------




#####  Input management methods  #####


# Define the data and target values matrices needed for the machine learning algorithms

"""no parameters""" # uses global variables

# returns: the data matrix and target values column needed for the machine learning algorithms to be applied
def dataTargetMatrices():
    # If we want to access the values distribution per mice group, the control and sunitinib groups
    miceSamplesGenesInfo.loc[:, 'MICE_GROUP'] = miceSamplesGenesInfo.index.str[3:-2]
    
    # If we want to access the values distribution per mice identifier
    miceSamplesGenesInfo.loc[:, 'MICE_ID'] = miceSamplesGenesInfo.index.str[:3]
    
    # If we want to access the values distribution per mice sample
    miceSamplesGenesInfo.loc[:, 'MICE_SAMPLE'] = miceSamplesGenesInfo.index
    

    # The mice samples input values to apply the machine learning algorithms  
    inputMice = miceSamplesGenesInfo.loc[:, miceSamplesGenesInfo.columns[:-3]].values
    
    # The target value for the mice samples machine learning algorithms
    miceExperimentGroup = list(miceSamplesGenesInfo['MICE_GROUP'])

    

    # The patients TCGA input values to apply the machine learning algorithms 
    inputPatientsTCGA = patientsGenesInfoTCGA.loc[: , patientsGenesInfoTCGA.columns[:-1]].values

    # The target value for the patients TCGA to apply the machine learning algorithms
    patientsTumorType = list(patientsGenesInfoTCGA['TUMOR_TYPE'])
    


    # Define control and sunitinib groups in the mice and patients TCGA samples genes dataframe 
    miceSamplesPatientsTCGA.loc[miceSamplesPatientsTCGA.index.str[3:-2] == 'Ctrl', 'TUMOR_TYPE'] = 'MICE_BRCA_CTRL'
    miceSamplesPatientsTCGA.loc[miceSamplesPatientsTCGA.index.str[3:-2] == 'Sunit', 'TUMOR_TYPE'] = 'MICE_BRCA_SUNIT'


    # The mice samples and patients TCGA input values to apply the machine learning algorithms
    inputMicePatientsTCGA = miceSamplesPatientsTCGA.loc[: , miceSamplesPatientsTCGA.columns[:-1]].values

    # The target values for the mice samples and patients TCGA when using both data together to apply the machine learning algorithms
    miceBothTumorType = list(miceSamplesPatientsTCGA['TUMOR_TYPE'])[:16]
    patientsBothTumorType = list(miceSamplesPatientsTCGA['TUMOR_TYPE'])[16:]


    return inputMice, inputPatientsTCGA, inputMicePatientsTCGA, miceExperimentGroup, patientsTumorType, miceBothTumorType, patientsBothTumorType




#--------------------------------------------------------------------------------------------------




#####  Mice samples data dimensionality reduction methods  #####


def micePCA2D():
  pca = PCA(n_components = 2)         
  micePCA = pca.fit_transform(inputMice)


  micePcaDF = pd.DataFrame(data = micePCA, columns = ['PC_1', 'PC_2'])

  # target value
  micePcaDF['MICE_GROUP'] = miceExperimentGroup


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  sns.scatterplot(data = micePcaDF, x = 'PC_1', y = 'PC_2', hue = 'MICE_GROUP', palette = ['b', 'r'], ax = ax, s = 60.0)

  group_legend = plt.gca().legend(bbox_to_anchor=(1.0, 0.5), loc='center left', borderaxespad=1.0, fontsize = 14,  \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)    


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([50.0]), group_legend.legendHandles))


  ax.set_xlabel('1' + '$^{st}$' + ' principal component', fontsize = 14.5, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' principal component', fontsize = 14.5, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14)
      

  # To then use to save the figures
  figure_name = "Mice samples applying PCA with 2 principal components"

  plt.suptitle(figure_name, fontsize = 15, y = 0.94).set_weight("bold")

  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def micePCA3D():
  pca = PCA(n_components = 3)         
  micePCA = pca.fit_transform(inputMice)


  micePcaDF = pd.DataFrame(data = micePCA, columns = ['PC_1', 'PC_2', 'PC_3'])

  # target value
  micePcaDF['MICE_GROUP'] = miceExperimentGroup


  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  fig.set_size_inches(16,10)


  list(map(lambda experimentGroup, color: ax.scatter(micePcaDF[micePcaDF['MICE_GROUP'] == experimentGroup]['PC_1'], 
                                                     micePcaDF[micePcaDF['MICE_GROUP'] == experimentGroup]['PC_2'],
                                                     micePcaDF[micePcaDF['MICE_GROUP'] == experimentGroup]['PC_3'],
                                                     c = color, s = 50.0), list(micePcaDF['MICE_GROUP'].unique()), ['b', 'r']))


  ax.set_xlabel('1' + '$^{st}$' + ' principal component', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' principal component', fontsize = 14.5, labelpad = 15)
  ax.set_zlabel('3' + '$^{rd}$' + ' principal component', fontsize = 14.5, labelpad = 15)

  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)


  group_legend = plt.gca().legend(bbox_to_anchor=(1.0, 0.5), loc='center left', borderaxespad=1.0, fontsize = 14,   \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)     
                  


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([50.0]), group_legend.legendHandles))


  # To then use to save the figures
  figure_name = "Mice samples applying PCA with 3 principal components"

  plt.suptitle(figure_name, fontsize = 15, y = 0.88).set_weight("bold")

  # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def miceTruncatedSVD2D():

  # Dimensionality reduction using truncated SVD

  svd = TruncatedSVD(n_components = 2)         
  miceSVD = svd.fit_transform(inputMice)


  miceSvdDF = pd.DataFrame(data = miceSVD, columns = ['SVD_1', 'SVD_2'])

  # target value
  miceSvdDF['MICE_GROUP'] = miceExperimentGroup


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  sns.scatterplot(data = miceSvdDF, x = 'SVD_1', y = 'SVD_2', hue = 'MICE_GROUP', palette = ['b', 'r'], ax = ax, s = 60.0)

  group_legend = plt.gca().legend(bbox_to_anchor=(1.0, 0.5), loc='center left', borderaxespad=1.0, fontsize = 14,  \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)    


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([50.0]), group_legend.legendHandles))


  ax.set_xlabel('1' + '$^{st}$' + ' truncated svd component', fontsize = 14.5, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' truncated svd component', fontsize = 14.5, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14)
      

  # To then use to save the figures
  figure_name = "Mice samples applying truncated SVD with 2 principal components"

  plt.suptitle(figure_name, fontsize = 15, y = 0.94).set_weight("bold")

  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()
    



def miceTruncatedSVD3D():

  svd = TruncatedSVD(n_components = 3)         
  miceSVD = svd.fit_transform(inputMice)


  miceSvdDF = pd.DataFrame(data = miceSVD, columns = ['SVD_1', 'SVD_2', 'SVD_3'])

  # target value
  miceSvdDF['MICE_GROUP'] = miceExperimentGroup


  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  fig.set_size_inches(16,10)



  list(map(lambda experimentGroup, color: ax.scatter(miceSvdDF[miceSvdDF['MICE_GROUP'] == experimentGroup]['SVD_1'], 
                                                     miceSvdDF[miceSvdDF['MICE_GROUP'] == experimentGroup]['SVD_2'],
                                                     miceSvdDF[miceSvdDF['MICE_GROUP'] == experimentGroup]['SVD_3'],
                                                     c = color, s = 50.0), list(miceSvdDF['MICE_GROUP'].unique()), ['b', 'r']))



  ax.set_xlabel('1' + '$^{st}$' + ' truncated svd component', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' truncated svd component', fontsize = 14.5, labelpad = 15)
  ax.set_zlabel('3' + '$^{rd}$' + ' truncated svd component', fontsize = 14.5, labelpad = 15)

  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)


  group_legend = plt.gca().legend(bbox_to_anchor=(1.0, 0.5), loc='center left', borderaxespad=1.0, fontsize = 14,   \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)     
                  


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([50.0]), group_legend.legendHandles))


  # To then use to save the figures
  figure_name = "Mice samples applying truncated SVD with 3 principal components"

  plt.suptitle(figure_name, fontsize = 15, y = 0.88).set_weight("bold")

  # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def miceTSNE2D():

  # T-distributed Stochastic Neighbor Embedding

  tsne = TSNE(n_components= 2, perplexity=10) 
  miceTSNE = tsne.fit_transform(inputMice)
      

  miceTsneDF = pd.DataFrame(data = miceTSNE, columns = ['TSNE_1', 'TSNE_2'])

  # target value
  miceTsneDF['MICE_GROUP'] = miceExperimentGroup


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  sns.scatterplot(data = miceTsneDF, x = 'TSNE_1', y = 'TSNE_2', hue = 'MICE_GROUP', palette = ['b', 'r'], ax = ax, s = 60.0)
                      

  group_legend = plt.gca().legend(bbox_to_anchor=(1.0, 0.5), loc='center left', borderaxespad=1.0, fontsize = 14,  \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)    


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([50.0]), group_legend.legendHandles))


  ax.set_xlabel('1' + '$^{st}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14)
      

  # To then use to save the figures
  figure_name = "Mice samples applying TSNE with 2 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.94).set_weight("bold")

  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()




def miceTSNE3D():

  # T-distributed Stochastic Neighbor Embedding

  tsne = TSNE(n_components= 3, perplexity=10) 
  miceTSNE = tsne.fit_transform(inputMice)
      

  miceTsneDF = pd.DataFrame(data = miceTSNE, columns = ['TSNE_1', 'TSNE_2', 'TSNE_3'])

  # target value
  miceTsneDF['MICE_GROUP'] = miceExperimentGroup



  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  fig.set_size_inches(16,10)



  list(map(lambda experimentGroup, color: ax.scatter(miceTsneDF[miceTsneDF['MICE_GROUP'] == experimentGroup]['TSNE_1'], 
                                                     miceTsneDF[miceTsneDF['MICE_GROUP'] == experimentGroup]['TSNE_2'],
                                                     miceTsneDF[miceTsneDF['MICE_GROUP'] == experimentGroup]['TSNE_3'],
                                                     c = color, s = 50.0), list(miceTsneDF['MICE_GROUP'].unique()), ['b', 'r']))



  ax.set_xlabel('1' + '$^{st}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 15)
  ax.set_zlabel('3' + '$^{rd}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 15)

  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)


  group_legend = plt.gca().legend(bbox_to_anchor=(1.0, 0.5), loc='center left', borderaxespad=1.0, fontsize = 14,   \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)     
                  


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([50.0]), group_legend.legendHandles))


  # To then use to save the figures
  figure_name = "Mice samples applying TSNE with 3 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.88).set_weight("bold")

  # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def micePcaTSNE2D():

  # T-distributed Stochastic Neighbor Embedding using PCA first to 16 dimensions
  pca = PCA(n_components = 16)         
  micePCA = pca.fit_transform(inputMice)

  tsne = TSNE(n_components= 2, perplexity=10) 
  miceTSNE = tsne.fit_transform(micePCA)
      
  miceTsneDF = pd.DataFrame(data = miceTSNE, columns = ['TSNE_1', 'TSNE_2'])

  # target value
  miceTsneDF['MICE_GROUP'] = miceExperimentGroup


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  sns.scatterplot(data = miceTsneDF, x = 'TSNE_1', y = 'TSNE_2', hue = 'MICE_GROUP', palette = ['b', 'r'], ax = ax, s = 60.0)
                      

  group_legend = plt.gca().legend(bbox_to_anchor=(1.0, 0.5), loc='center left', borderaxespad=1.0, fontsize = 14,  \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)    


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([50.0]), group_legend.legendHandles))


  ax.set_xlabel('1' + '$^{st}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14)
      

  # To then use to save the figures
  figure_name = "Mice samples applying PCA to 16 dimensions followed by TSNE with 2 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.94).set_weight("bold")

  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()




def micePcaTSNE3D():

  pca = PCA(n_components = 16)         
  micePCA = pca.fit_transform(inputMice)

  tsne = TSNE(n_components= 3, perplexity=10) 
  miceTSNE = tsne.fit_transform(micePCA)
      

  miceTsneDF = pd.DataFrame(data = miceTSNE, columns = ['TSNE_1', 'TSNE_2', 'TSNE_3'])

  # target value
  miceTsneDF['MICE_GROUP'] = miceExperimentGroup



  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  fig.set_size_inches(16,10)



  list(map(lambda experimentGroup, color: ax.scatter(miceTsneDF[miceTsneDF['MICE_GROUP'] == experimentGroup]['TSNE_1'], 
                                                     miceTsneDF[miceTsneDF['MICE_GROUP'] == experimentGroup]['TSNE_2'],
                                                     miceTsneDF[miceTsneDF['MICE_GROUP'] == experimentGroup]['TSNE_3'],
                                                     c = color, s = 50.0), list(miceTsneDF['MICE_GROUP'].unique()), ['b', 'r']))



  ax.set_xlabel('1' + '$^{st}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 15)
  ax.set_zlabel('3' + '$^{rd}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 15)

  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)


  group_legend = plt.gca().legend(bbox_to_anchor=(1.0, 0.5), loc='center left', borderaxespad=1.0, fontsize = 14,   \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)     
                  


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([50.0]), group_legend.legendHandles))


  # To then use to save the figures
  figure_name = "Mice samples applying PCA to 16 dimensions followed by TSNE with 3 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.88).set_weight("bold")

  # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()




def miceMDS2D():

  # Multidimensional scaling

  mds = MDS(n_components = 2)         
  miceMDS = mds.fit_transform(inputMice)


  miceMdsDF = pd.DataFrame(data = miceMDS, columns = ['MDS_1', 'MDS_2'])

  # target value
  miceMdsDF['MICE_GROUP'] = miceExperimentGroup


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  sns.scatterplot(data = miceMdsDF, x = 'MDS_1', y = 'MDS_2', hue = 'MICE_GROUP', palette = ['b', 'r'], ax = ax, s = 60.0)
                      

  group_legend = plt.gca().legend(bbox_to_anchor=(1.0, 0.5), loc='center left', borderaxespad=1.0, fontsize = 14,  \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)    


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([50.0]), group_legend.legendHandles))


  ax.set_xlabel('1' + '$^{st}$' + ' dimension of MDS', fontsize = 14.5, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of MDS', fontsize = 14.5, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14)
      

  # To then use to save the figures
  figure_name = "Mice samples applying MDS with 2 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.94).set_weight("bold")

  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()




def miceMDS3D():

  # Multidimensional scaling

  mds = MDS(n_components = 3)         
  miceMDS = mds.fit_transform(inputMice)


  miceMdsDF = pd.DataFrame(data = miceMDS, columns = ['MDS_1', 'MDS_2', 'MDS_3'])

  # target value
  miceMdsDF['MICE_GROUP'] = miceExperimentGroup



  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  fig.set_size_inches(16,10)



  list(map(lambda experimentGroup, color: ax.scatter(miceMdsDF[miceMdsDF['MICE_GROUP'] == experimentGroup]['MDS_1'], 
                                                     miceMdsDF[miceMdsDF['MICE_GROUP'] == experimentGroup]['MDS_2'],
                                                     miceMdsDF[miceMdsDF['MICE_GROUP'] == experimentGroup]['MDS_3'],
                                                     c = color, s = 50.0), list(miceMdsDF['MICE_GROUP'].unique()), ['b', 'r']))


  ax.set_xlabel('1' + '$^{st}$' + ' dimension of MDS', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of MDS', fontsize = 14.5, labelpad = 15)
  ax.set_zlabel('3' + '$^{rd}$' + ' dimension of MDS', fontsize = 14.5, labelpad = 15)

  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)


  group_legend = plt.gca().legend(bbox_to_anchor=(1.0, 0.5), loc='center left', borderaxespad=1.0, fontsize = 14,   \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)     
                  


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([50.0]), group_legend.legendHandles))


  # To then use to save the figures
  figure_name = "Mice samples applying MDS with 3 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.88).set_weight("bold")

  # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()




def micePcaMDS2D():

  # Multidimensional scaling

  pca = PCA(n_components = 16)         
  micePCA = pca.fit_transform(inputMice)


  mds = MDS(n_components = 2)         
  miceMDS = mds.fit_transform(micePCA)


  miceMdsDF = pd.DataFrame(data = miceMDS, columns = ['MDS_1', 'MDS_2'])

  # target value
  miceMdsDF['MICE_GROUP'] = miceExperimentGroup


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  sns.scatterplot(data = miceMdsDF, x = 'MDS_1', y = 'MDS_2', hue = 'MICE_GROUP', palette = ['b', 'r'], ax = ax, s = 60.0)
                      

  group_legend = plt.gca().legend(bbox_to_anchor=(1.0, 0.5), loc='center left', borderaxespad=1.0, fontsize = 14,  \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)    


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([50.0]), group_legend.legendHandles))


  ax.set_xlabel('1' + '$^{st}$' + ' dimension of MDS', fontsize = 14.5, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of MDS', fontsize = 14.5, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14)
      

  # To then use to save the figures
  figure_name = "Mice samples applying PCA to 16 dimensions followed by MDS with 2 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.94).set_weight("bold")

  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def micePcaMDS3D():

  # Multidimensional scaling

  pca = PCA(n_components = 16)         
  micePCA = pca.fit_transform(inputMice)

  mds = MDS(n_components = 3)         
  miceMDS = mds.fit_transform(micePCA)


  miceMdsDF = pd.DataFrame(data = miceMDS, columns = ['MDS_1', 'MDS_2', 'MDS_3'])

  # target value
  miceMdsDF['MICE_GROUP'] = miceExperimentGroup



  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  fig.set_size_inches(16,10)



  list(map(lambda experimentGroup, color: ax.scatter(miceMdsDF[miceMdsDF['MICE_GROUP'] == experimentGroup]['MDS_1'], 
                                                     miceMdsDF[miceMdsDF['MICE_GROUP'] == experimentGroup]['MDS_2'],
                                                     miceMdsDF[miceMdsDF['MICE_GROUP'] == experimentGroup]['MDS_3'],
                                                     c = color, s = 50.0), list(miceMdsDF['MICE_GROUP'].unique()), ['b', 'r']))


  ax.set_xlabel('1' + '$^{st}$' + ' dimension of MDS', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of MDS', fontsize = 14.5, labelpad = 15)
  ax.set_zlabel('3' + '$^{rd}$' + ' dimension of MDS', fontsize = 14.5, labelpad = 15)

  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)


  group_legend = plt.gca().legend(bbox_to_anchor=(1.0, 0.5), loc='center left', borderaxespad=1.0, fontsize = 14,   \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)     
                  


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([50.0]), group_legend.legendHandles))


  # To then use to save the figures
  figure_name = "Mice samples applying PCA to 16 dimensions followed by MDS with 3 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.88).set_weight("bold")

  # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()




def miceIsomap2D():

  # Non-linear dimensionality reduction through Isometric Mapping

  isomap = Isomap(n_components= 2)    
  miceIsomap = isomap.fit_transform(inputMice)    


  miceIsomapDF = pd.DataFrame(data = miceIsomap, columns = ['ISOMAP_1', 'ISOMAP_2'])

  # target value
  miceIsomapDF['MICE_GROUP'] = miceExperimentGroup


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  sns.scatterplot(data = miceIsomapDF, x = 'ISOMAP_1', y = 'ISOMAP_2', hue = 'MICE_GROUP', palette = ['b', 'r'], ax = ax, s = 60.0)
                      

  group_legend = plt.gca().legend(bbox_to_anchor=(1.0, 0.5), loc='center left', borderaxespad=1.0, fontsize = 14,  \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)    


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([50.0]), group_legend.legendHandles))


  ax.set_xlabel('1' + '$^{st}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14)
      

  # To then use to save the figures
  figure_name = "Mice samples applying Isomap with 2 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.94).set_weight("bold")

  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()




def miceIsomap3D():

  # Non-linear dimensionality reduction through Isometric Mapping

  isomap = Isomap(n_components = 3)    
  miceIsomap = isomap.fit_transform(inputMice)    


  miceIsomapDF = pd.DataFrame(data = miceIsomap, columns = ['ISOMAP_1', 'ISOMAP_2', 'ISOMAP_3'])

  # target value
  miceIsomapDF['MICE_GROUP'] = miceExperimentGroup



  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  fig.set_size_inches(16,10)



  list(map(lambda experimentGroup, color: ax.scatter(miceIsomapDF[miceIsomapDF['MICE_GROUP'] == experimentGroup]['ISOMAP_1'], 
                                                     miceIsomapDF[miceIsomapDF['MICE_GROUP'] == experimentGroup]['ISOMAP_2'],
                                                     miceIsomapDF[miceIsomapDF['MICE_GROUP'] == experimentGroup]['ISOMAP_3'],
                                                     c = color, s = 50.0), list(miceIsomapDF['MICE_GROUP'].unique()), ['b', 'r']))



  ax.set_xlabel('1' + '$^{st}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 15)
  ax.set_zlabel('3' + '$^{rd}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 15)

  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)


  group_legend = plt.gca().legend(bbox_to_anchor=(1.0, 0.5), loc='center left', borderaxespad=1.0, fontsize = 14,   \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)     
                  


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([50.0]), group_legend.legendHandles))


  # To then use to save the figures
  figure_name = "Mice samples applying Isomap with 3 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.88).set_weight("bold")

  # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()




# With 16 dimensions the Isomap applied after PCA would not have difference with applying the Isomap separately 
def micePcaIsomap2D():

  # Non-linear dimensionality reduction through Isometric Mapping
  pca = PCA(n_components = 6)         
  micePCA = pca.fit_transform(inputMice)

  isomap = Isomap(n_components= 2)    
  miceIsomap = isomap.fit_transform(micePCA)    


  miceIsomapDF = pd.DataFrame(data = miceIsomap, columns = ['ISOMAP_1', 'ISOMAP_2'])

  # target value
  miceIsomapDF['MICE_GROUP'] = miceExperimentGroup


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  sns.scatterplot(data = miceIsomapDF, x = 'ISOMAP_1', y = 'ISOMAP_2', hue = 'MICE_GROUP', palette = ['b', 'r'], ax = ax, s = 60.0)
                      

  group_legend = plt.gca().legend(bbox_to_anchor=(1.0, 0.5), loc='center left', borderaxespad=1.0, fontsize = 14,  \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)    


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([50.0]), group_legend.legendHandles))


  ax.set_xlabel('1' + '$^{st}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14)
      

  # To then use to save the figures
  figure_name = "Mice samples applying PCA to 6 dimensions followed by Isomap with 2 dimensions"  

  plt.suptitle(figure_name, fontsize = 15, y = 0.94).set_weight("bold")

  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()


   

def micePcaIsomap3D():
  # Non-linear dimensionality reduction through Isometric Mapping
  pca = PCA(n_components = 6)         
  micePCA = pca.fit_transform(inputMice)

  isomap = Isomap(n_components= 3)    
  miceIsomap = isomap.fit_transform(micePCA)    

  miceIsomapDF = pd.DataFrame(data = miceIsomap, columns = ['ISOMAP_1', 'ISOMAP_2', 'ISOMAP_3'])

  # target value
  miceIsomapDF['MICE_GROUP'] = miceExperimentGroup



  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  fig.set_size_inches(16,10)



  list(map(lambda experimentGroup, color: ax.scatter(miceIsomapDF[miceIsomapDF['MICE_GROUP'] == experimentGroup]['ISOMAP_1'], 
                                                     miceIsomapDF[miceIsomapDF['MICE_GROUP'] == experimentGroup]['ISOMAP_2'],
                                                     miceIsomapDF[miceIsomapDF['MICE_GROUP'] == experimentGroup]['ISOMAP_3'],
                                                     c = color, s = 50.0), list(miceIsomapDF['MICE_GROUP'].unique()), ['b', 'r']))



  ax.set_xlabel('1' + '$^{st}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 15)
  ax.set_zlabel('3' + '$^{rd}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 15)

  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)


  group_legend = plt.gca().legend(bbox_to_anchor=(1.0, 0.5), loc='center left', borderaxespad=1.0, fontsize = 14,   \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)     
                  


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([50.0]), group_legend.legendHandles))


  # To then use to save the figures
  figure_name = "Mice samples applying PCA to 6 dimensions followed by Isomap with 3 dimensions" 

  plt.suptitle(figure_name, fontsize = 15, y = 0.88).set_weight("bold")

    # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()




#--------------------------------------------------------------------------------------------------




####  Patients TCGA data dimensionality reduction methods  #####


def patientsPCA2D():
  pca = PCA(n_components = 2)         
  patientsPCA = pca.fit_transform(inputPatientsTCGA)

  patientsPcaDF = pd.DataFrame(data = patientsPCA, columns = ['PC_1', 'PC_2'])

  #target value
  patientsPcaDF['TUMOR_TYPE'] = patientsTumorType


  patientsPcaDF = patientsPcaDF.sort_values(['TUMOR_TYPE'])


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  g = sns.scatterplot(data = patientsPcaDF, x = 'PC_1', y = 'PC_2', hue = 'TUMOR_TYPE', ax = ax, s = 90.0)


  ax.set_xlabel('1' + '$^{st}$' + ' principal component', fontsize = 14.8, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' principal component', fontsize = 14.8, labelpad = 15)

  # default of axis is both
  ax.tick_params(labelsize = 14.4)
      
  legend_handles, _= g.get_legend_handles_labels()      

  tumor_legend = plt.gca().legend(bbox_to_anchor=(1.025, 0.5), loc='center left', borderaxespad=1.0, fontsize = 13.6, handles = legend_handles, \
      labels = list(patientsPcaDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14.3, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.2)    
                  


  figure_name = "Patients TCGA applying PCA with 2 principal components"

  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")

  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()


  return tumor_legend, patientsPCA




def patientsPcaDensity2D(patientsPCA):
  # Calculate the density of points, normalize it and encode it in the alpha channel of a colormap
  # Low density points are emphasized by darker color, while the high density points with the magma colormap

  patientsPCA = patientsPCA.T

  kernel = stats.gaussian_kde(patientsPCA)
  weights = kernel(patientsPCA)
  weights = weights/weights.max()


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  g = sns.scatterplot(x = patientsPCA[0], y = patientsPCA[1], hue = weights, palette="magma", ax = ax, s = 85.0)


  ax.set_xlabel('1' + '$^{st}$' + ' principal component', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' principal component', fontsize = 14.5, labelpad = 15)

  # default of axis is both
  ax.tick_params(labelsize = 14)


  g.get_legend().remove()

  sm = plt.cm.ScalarMappable(cmap="magma")
  sm.set_array([])

  cbar = g.figure.colorbar(sm, ticks = [0, 1])
  
  cbar.set_ticklabels(["Low", "High"])

  cbar.set_label(label='Density', size=14, weight='bold', rotation = 360, labelpad = -1, y = 0.5)
  cbar.ax.tick_params(labelsize=13.5)


  figure_name = "Patients TCGA applying PCA with 2 principal components - density zones"

  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")

  # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def patientsPcaKernel2D(patientsPCA):
  x = patientsPCA[:, 0]
  y = patientsPCA[:, 1]

  # Manually adjusted
  xmin, xmax = -8, 3
  ymin, ymax = -3, 6

  # Create meshgrid
  xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]


  positions = np.vstack([xx.ravel(), yy.ravel()])
  values = np.vstack([x, y])
  kernel = stats.gaussian_kde(values)
  f = np.reshape(kernel(positions).T, xx.shape)


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)
  
  ax.set_xlim(xmin, xmax)
  ax.set_ylim(ymin, ymax)
  cfset = ax.contourf(xx, yy, f, cmap='coolwarm')
  ax.imshow(np.rot90(f), cmap='coolwarm', extent=[xmin, xmax, ymin, ymax])
  cset = ax.contour(xx, yy, f, colors='k')
  ax.clabel(cset, inline=1, fontsize=12)
  ax.set_xlabel('1' + '$^{st}$' + ' principal component', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' principal component', fontsize = 14.5, labelpad = 15)

  # default of axis is both
  ax.tick_params(labelsize = 14)


  figure_name = "Patients TCGA applying PCA with 2 principal components - gaussian density estimation"
 
  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")

  # save figure
  figure_name = '3 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()


  return xx, yy, f




def patientsPCA3D(tumorLegend):

  pca = PCA(n_components = 3)         
  patientsPCA = pca.fit_transform(inputPatientsTCGA)


  patientsPcaDF = pd.DataFrame(data = patientsPCA, columns = ['PC_1', 'PC_2', 'PC_3'])

  #target value
  patientsPcaDF['TUMOR_TYPE'] = patientsTumorType


  patientsPcaDF = patientsPcaDF.sort_values(['TUMOR_TYPE'])


  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,10)


  list(map(lambda tumorType, legend: ax.scatter(patientsPcaDF[patientsPcaDF['TUMOR_TYPE'] == tumorType]['PC_1'], 
                                                patientsPcaDF[patientsPcaDF['TUMOR_TYPE'] == tumorType]['PC_2'], 
                                                patientsPcaDF[patientsPcaDF['TUMOR_TYPE'] == tumorType]['PC_3'], 
                                                s = 60.0, c = legend.get_facecolor()), 
                                                list(patientsPcaDF['TUMOR_TYPE'].unique()), tumorLegend.legendHandles))

  ax.set_xlabel('1' + '$^{st}$' + ' principal component', fontsize = 14.8, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' principal component', fontsize = 14.8, labelpad = 15)
  ax.set_zlabel('3' + '$^{rd}$' + ' principal component', fontsize = 14.8, labelpad = 15)

  ax.tick_params(labelsize = 14.3)
  ax.tick_params(labelsize = 14.3)
  ax.tick_params(labelsize = 14.3)


  group_legend = plt.gca().legend(bbox_to_anchor=(1.01, 0.5), loc='center left', borderaxespad=1.0, fontsize = 13.5,  \
      labels = list(patientsPcaDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14.5, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.6, handletextpad = 0.2)    
                  

  figure_name = "Patients TCGA applying PCA with 3 principal components"

  plt.suptitle(figure_name, fontsize = 15.3, y = 0.90).set_weight("bold")

  # save figure
  figure_name = '4 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def patientsPcaKernel3D(xx, yy, f):

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,10)

  surf = ax.plot_surface(xx, yy, f, rstride=1, cstride=1, cmap='coolwarm', edgecolor='none')
  ax.set_xlabel('1' + '$^{st}$' + ' principal component', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' principal component', fontsize = 14.5, labelpad = 15)
  ax.set_zlabel('Probability density function', fontsize = 14.5, labelpad = 28)
  

  ax.tick_params(axis = 'x', labelsize = 14)
  ax.tick_params(axis = 'y', labelsize = 14)
  ax.tick_params(axis = 'z', labelsize = 14, pad = 13.0)


  cbar = fig.colorbar(surf, shrink=0.5, aspect=5) # add color bar indicating the PDF

  cbar.ax.tick_params(labelsize=13)


  ax.view_init(10, 120)


  figure_name = "Patients TCGA applying PCA with 3 principal components - 2D probability density distribution"
 
  plt.suptitle(figure_name, fontsize = 15, y = 0.83).set_weight("bold")

  # save figure
  figure_name = '5 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def patientsTruncatedSVD2D():

  # Dimensionality reduction using truncated SVD

  svd = TruncatedSVD(n_components = 2)         
  patientsSVD = svd.fit_transform(inputPatientsTCGA)


  patientsSvdDF = pd.DataFrame(data = patientsSVD, columns = ['SVD_1', 'SVD_2'])

  # target value
  patientsSvdDF['TUMOR_TYPE'] = patientsTumorType


  patientsSvdDF = patientsSvdDF.sort_values(['TUMOR_TYPE'])


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  g = sns.scatterplot(data = patientsSvdDF, x = 'SVD_1', y = 'SVD_2', hue = 'TUMOR_TYPE', ax = ax, s = 90.0)


  ax.set_xlabel('1' + '$^{st}$' + ' truncated svd component', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' truncated svd component', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14.4)
      

  legend_handles, _= g.get_legend_handles_labels()      

  tumor_legend = plt.gca().legend(bbox_to_anchor=(1.025, 0.5), loc='center left', borderaxespad=1.0, fontsize = 13.6, handles = legend_handles, \
      labels = list(patientsSvdDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14.3, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.2)    
                  

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([55.0]), tumor_legend.legendHandles))


  figure_name = "Patients TCGA applying truncated SVD with 2 principal components"

  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")

  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()


  return tumor_legend, patientsSVD




def patientsTruncatedSvdDensity2D(patientsSVD):
  # Calculate the density of points, normalize it and encode it in the alpha channel of a colormap
  # Low density points are emphasized by darker color, while the high density points with the magma colormap

  patientsSVD = patientsSVD.T

  kernel = stats.gaussian_kde(patientsSVD)
  weights = kernel(patientsSVD)
  weights = weights/weights.max()


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  g = sns.scatterplot(x = patientsSVD[0], y = patientsSVD[1], hue = weights, palette="magma", ax = ax, s = 85.0)


  ax.set_xlabel('1' + '$^{st}$' + ' truncated svd component', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' truncated svd component', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14.4)


  g.get_legend().remove()

  sm = plt.cm.ScalarMappable(cmap="magma")
  sm.set_array([])

  cbar = g.figure.colorbar(sm, ticks = [0, 1])
  
  cbar.set_ticklabels(["Low", "High"])

  cbar.set_label(label='Density', size=14, weight='bold', rotation = 360, labelpad = -1, y = 0.5)
  cbar.ax.tick_params(labelsize=13.5)
  

  figure_name = "Patients TCGA applying truncated SVD with 2 principal components - density zones"

  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")

  # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def patientsTruncatedSvdKernel2D(patientsSVD):
  x = patientsSVD[:, 0]
  y = patientsSVD[:, 1]

  # Manually adjusted
  xmin, xmax = -4, 8
  ymin, ymax = -4, 4

  # Create meshgrid
  xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]


  positions = np.vstack([xx.ravel(), yy.ravel()])
  values = np.vstack([x, y])
  kernel = stats.gaussian_kde(values)
  f = np.reshape(kernel(positions).T, xx.shape)


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)
  
  ax.set_xlim(xmin, xmax)
  ax.set_ylim(ymin, ymax)
  cfset = ax.contourf(xx, yy, f, cmap='coolwarm')
  ax.imshow(np.rot90(f), cmap='coolwarm', extent=[xmin, xmax, ymin, ymax])
  cset = ax.contour(xx, yy, f, colors='k')
  ax.clabel(cset, inline=1, fontsize=12)
  ax.set_xlabel('1' + '$^{st}$' + ' truncated svd component', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' truncated svd component', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14.4)


  figure_name = "Patients TCGA applying truncated SVD with 2 principal components - gaussian density estimation"
 
  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")

  # save figure
  figure_name = '3 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()


  return xx, yy, f




def patientsTruncatedSVD3D(tumorLegend):

  # Dimensionality reduction using truncated SVD

  svd = TruncatedSVD(n_components = 3)         
  patientsSVD = svd.fit_transform(inputPatientsTCGA)

  patientsSvdDF = pd.DataFrame(data = patientsSVD, columns = ['SVD_1', 'SVD_2', 'SVD_3'])

  #target value
  patientsSvdDF['TUMOR_TYPE'] = patientsTumorType


  patientsSvdDF = patientsSvdDF.sort_values(['TUMOR_TYPE'])


  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,10)


  list(map(lambda tumorType, legend: ax.scatter(patientsSvdDF[patientsSvdDF['TUMOR_TYPE'] == tumorType]['SVD_1'], 
                                                patientsSvdDF[patientsSvdDF['TUMOR_TYPE'] == tumorType]['SVD_2'], 
                                                patientsSvdDF[patientsSvdDF['TUMOR_TYPE'] == tumorType]['SVD_3'], 
                                                s = 60.0, c = legend.get_facecolor()), 
                                                list(patientsSvdDF['TUMOR_TYPE'].unique()), tumorLegend.legendHandles))

  ax.set_xlabel('1' + '$^{st}$' + ' truncated svd component', fontsize = 14.8, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' truncated svd component', fontsize = 14.8, labelpad = 15)
  ax.set_zlabel('3' + '$^{rd}$' + ' truncated svd component', fontsize = 14.8, labelpad = 15)

  ax.tick_params(labelsize = 14.3)
  ax.tick_params(labelsize = 14.3)
  ax.tick_params(labelsize = 14.3)


  group_legend = plt.gca().legend(bbox_to_anchor=(1.01, 0.5), loc='center left', borderaxespad=1.0, fontsize = 13.5,  \
      labels = list(patientsSvdDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14.5, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.6, handletextpad = 0.2)    
                  

  figure_name = "Patients TCGA applying truncated SVD with 3 principal components"

  plt.suptitle(figure_name, fontsize = 15, y = 0.90).set_weight("bold")

  # save figure
  figure_name = '4 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def patientsTruncatedSvdKernel3D(xx, yy, f):

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,10)

  surf = ax.plot_surface(xx, yy, f, rstride=1, cstride=1, cmap='coolwarm', edgecolor='none')
  ax.set_xlabel('1' + '$^{st}$' + ' truncated svd component', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' truncated svd component', fontsize = 14.5, labelpad = 15)
  ax.set_zlabel('Probability density function', fontsize = 14.5, labelpad = 28)
  

  ax.tick_params(axis = 'x', labelsize = 14)
  ax.tick_params(axis = 'y', labelsize = 14)
  ax.tick_params(axis = 'z', labelsize = 14, pad = 13.0)


  cbar = fig.colorbar(surf, shrink=0.5, aspect=5) # add color bar indicating the PDF

  cbar.ax.tick_params(labelsize=13)


  ax.view_init(10, 120)


  figure_name = "Patients TCGA applying truncated SVD with 3 principal components - 2D probability density distribution"
 
  plt.suptitle(figure_name, fontsize = 15, y = 0.83).set_weight("bold")

  # save figure
  figure_name = '5 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()




# TSNE alone with big data takes a lot of time and memory
def patientsPca25TSNE2D():

  # Using PCA first, followed by tSNE (2D) 25 dimensions
  # T-distributed Stochastic Neighbor Embedding

  pca = PCA(n_components = 25)         
  patientsPCA = pca.fit_transform(inputPatientsTCGA)

  tsne = TSNE(n_components = 2, perplexity=40, n_iter=300)         
  patientsTSNE = tsne.fit_transform(patientsPCA)


  patientsTsneDF = pd.DataFrame(data = patientsTSNE, columns = ['TSNE_1', 'TSNE_2'])

  # target value
  patientsTsneDF['TUMOR_TYPE'] = patientsTumorType


  patientsTsneDF = patientsTsneDF.sort_values(['TUMOR_TYPE'])


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  g = sns.scatterplot(data = patientsTsneDF, x = 'TSNE_1', y = 'TSNE_2', hue = 'TUMOR_TYPE', ax = ax, s = 90.0)


  ax.set_xlabel('1' + '$^{st}$' + ' dimension of TSNE embedded space', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of TSNE embedded space', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14.4)
      

  legend_handles, _= g.get_legend_handles_labels()      

  tumor_legend = plt.gca().legend(bbox_to_anchor=(1.10, 0.5), loc='center left', borderaxespad=1.0, fontsize = 13.4, handles = legend_handles, \
      labels = list(patientsTsneDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.2)    
                  

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([55.0]), tumor_legend.legendHandles))

  figure_name = "Patients TCGA applying PCA first to 25 components followed by TSNE with 2 dimensions"

  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")

  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()


  return tumor_legend, patientsTSNE




def patientsPca25TsneDensity2D(patientsTSNE):

  # Calculate the density of points, normalize it and encode it in the alpha channel of a colormap
  # Low density points are emphasized by darker color, while the high density points with the magma colormap

  patientsTSNE = patientsTSNE.T

  kernel = stats.gaussian_kde(patientsTSNE)
  weights = kernel(patientsTSNE)
  weights = weights/weights.max()


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  g = sns.scatterplot(x = patientsTSNE[0], y = patientsTSNE[1], hue = weights, palette="magma", ax = ax, s = 85.0)


  ax.set_xlabel('1' + '$^{st}$' + ' dimension of TSNE embedded space', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of TSNE embedded space', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14.3)

  
  g.get_legend().remove()

  sm = plt.cm.ScalarMappable(cmap="magma")
  sm.set_array([])

  cbar = g.figure.colorbar(sm, ticks = [0, 1])
  
  cbar.set_ticklabels(["Low", "High"])

  cbar.set_label(label='Density', size=14, weight='bold', rotation = 360, labelpad = -1, y = 0.5)
  cbar.ax.tick_params(labelsize=13.5)


  figure_name = "Patients TCGA applying PCA first to 25 components followed by TSNE with 2 dimensions - density zones"

  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")


  # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()




def patientsPca25TsneKernel2D(patientsTSNE): 

  x = patientsTSNE[:, 0]
  y = patientsTSNE[:, 1]

  # Manually adjusted
  xmin, xmax = -13, 13 
  ymin, ymax = -13, 13

  # Create meshgrid
  xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]


  positions = np.vstack([xx.ravel(), yy.ravel()])
  values = np.vstack([x, y])
  kernel = stats.gaussian_kde(values)
  f = np.reshape(kernel(positions).T, xx.shape)


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)
  
  ax.set_xlim(xmin, xmax)
  ax.set_ylim(ymin, ymax)
  cfset = ax.contourf(xx, yy, f, cmap='coolwarm')
  ax.imshow(np.rot90(f), cmap='coolwarm', extent=[xmin, xmax, ymin, ymax])
  cset = ax.contour(xx, yy, f, colors='k')
  ax.clabel(cset, inline=1, fontsize=12)
  ax.set_xlabel('1' + '$^{st}$' + ' dimension of TSNE embedded space', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of TSNE embedded space', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14.3)


  figure_name = "Patients TCGA applying PCA first to 25 components followed by TSNE with 2 dimensions - gaussian density estimation"
 
  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")

  # save figure
  figure_name = '3 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()


  return xx, yy, f




def patientsPca25TSNE3D(tumorLegend):

  # Using PCA first, followed by tSNE (3D) 25 dimensions
  # T-distributed Stochastic Neighbor Embedding

  pca = PCA(n_components = 25)         
  patientsPCA = pca.fit_transform(inputPatientsTCGA)

  tsne = TSNE(n_components = 3, perplexity=40, n_iter=300)         
  patientsTSNE = tsne.fit_transform(patientsPCA)


  patientsTsneDF = pd.DataFrame(data = patientsTSNE, columns = ['TSNE_1', 'TSNE_2', 'TSNE_3'])

  #target value
  patientsTsneDF['TUMOR_TYPE'] = patientsTumorType


  patientsTsneDF = patientsTsneDF.sort_values(['TUMOR_TYPE'])


  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,10)


  list(map(lambda tumorType, legend: ax.scatter(patientsTsneDF[patientsTsneDF['TUMOR_TYPE'] == tumorType]['TSNE_1'], 
                                                patientsTsneDF[patientsTsneDF['TUMOR_TYPE'] == tumorType]['TSNE_2'], 
                                                patientsTsneDF[patientsTsneDF['TUMOR_TYPE'] == tumorType]['TSNE_3'], 
                                                s = 50.0, c = legend.get_facecolor()), 
                                                list(patientsTsneDF['TUMOR_TYPE'].unique()), tumorLegend.legendHandles))

  ax.set_xlabel('1' + '$^{st}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 14.5)
  ax.set_zlabel('3' + '$^{rd}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 14.5)

  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)


  group_legend = plt.gca().legend(bbox_to_anchor=(1.01, 0.5), loc='center left', borderaxespad=1.0, fontsize = 13.5,  \
      labels = list(patientsTsneDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14.5, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.6, handletextpad = 0.2)    
                  

  figure_name = "Patients TCGA applying PCA first to 25 components and TSNE with 3 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.90).set_weight("bold")


  # save figure
  figure_name = '4 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()




def patientsPca25TsneKernel3D(xx, yy, f):

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,12)

  surf = ax.plot_surface(xx, yy, f, rstride=1, cstride=1, cmap='coolwarm', edgecolor='none')
  ax.set_xlabel('1' + '$^{st}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 14.5)
  ax.set_zlabel('Probability density function', fontsize = 14.5, labelpad = 28)


  ax.tick_params(axis = 'x', labelsize = 14)
  ax.tick_params(axis = 'y', labelsize = 14, labelleft = False, labelright = True, pad = 13.0)
  ax.tick_params(axis = 'z', labelsize = 14, pad = 13.0)


  cbar = fig.colorbar(surf, shrink=0.5, aspect=5) # add color bar indicating the PDF

  cbar.ax.tick_params(labelsize=13)

  ax.view_init(elev=21)


  figure_name = "Patients TCGA applying PCA first to 25 components followed by TSNE with 3 dimensions - 2D probability density distribution"

  plt.suptitle(figure_name, fontsize = 15, y = 0.86).set_weight("bold")

  # save figure
  figure_name = '5 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




# MDS alone with big data takes a lot of time and memory
# MDS is not really efficient in really big datasets

def patientsPca25MDS2D(): 
  
  # Using PCA first to 25 dimensions, followed by MDS (2D)
  # Multidimensional scaling

  pca = PCA(n_components = 25)         
  patientsPCA = pca.fit_transform(inputPatientsTCGA)

  mds = MDS(n_components = 2)         
  patientsMDS = mds.fit_transform(patientsPCA)


  patientsMdsDF = pd.DataFrame(data = patientsMDS, columns = ['MDS_1', 'MDS_2'])


  #target value
  patientsMdsDF['TUMOR_TYPE'] = patientsTumorType

  patientsMdsDF = patientsMdsDF.sort_values(['TUMOR_TYPE'])


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  g = sns.scatterplot(data = patientsMdsDF, x = 'MDS_1', y = 'MDS_2', hue = 'TUMOR_TYPE', ax = ax, s = 90.0)


  ax.set_xlabel('1' + '$^{st}$' + ' dimension of MDS', fontsize = 14.8, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of MDS', fontsize = 14.8, labelpad = 15)

  # default of axis is both
  ax.tick_params(labelsize = 14.3)
      
  legend_handles, _= g.get_legend_handles_labels()      

  tumor_legend = plt.gca().legend(bbox_to_anchor=(1.1, 0.5), loc='center left', borderaxespad=1.0, fontsize = 13.6, handles = legend_handles, \
      labels = list(patientsMdsDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14.3, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.2)    
                  

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([55.0]), tumor_legend.legendHandles))

  figure_name = "Patients TCGA applying PCA first to 25 components followed by MDS with 2 dimensions"

  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")


  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()


  return tumor_legend, patientsMDS




def patientsPca25MdsDensity2D(patientsMDS):

  # Calculate the density of points, normalize it and encode it in the alpha channel of a colormap
  # Low density points are emphasized by darker color, while the high density points with the magma colormap

  patientsMDS = patientsMDS.T

  kernel = stats.gaussian_kde(patientsMDS)
  weights = kernel(patientsMDS)
  weights = weights/weights.max()


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  g = sns.scatterplot(x = patientsMDS[0], y = patientsMDS[1], hue = weights, palette="magma", ax = ax, s = 85.0)


  ax.set_xlabel('1' + '$^{st}$' + ' dimension of MDS', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of MDS', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14.3)

  
  g.get_legend().remove()

  sm = plt.cm.ScalarMappable(cmap="magma")
  sm.set_array([])

  cbar = g.figure.colorbar(sm, ticks = [0, 1])
  
  cbar.set_ticklabels(["Low", "High"])

  cbar.set_label(label='Density', size=14, weight='bold', rotation = 360, labelpad = -1, y = 0.5)
  cbar.ax.tick_params(labelsize=13.5)


  figure_name = "Patients TCGA applying PCA first to 25 components followed by MDS with 2 dimensions - density zones"

  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")

  # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def patientsPca25MdsKernel2D(patientsMDS):

  x = patientsMDS[:, 0]
  y = patientsMDS[:, 1]

  # Manually adjusted
  xmin, xmax = -13, 13 
  ymin, ymax = -13, 13

  # Create meshgrid
  xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]


  positions = np.vstack([xx.ravel(), yy.ravel()])
  values = np.vstack([x, y])
  kernel = stats.gaussian_kde(values)
  f = np.reshape(kernel(positions).T, xx.shape)


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)
  
  ax.set_xlim(xmin, xmax)
  ax.set_ylim(ymin, ymax)
  cfset = ax.contourf(xx, yy, f, cmap='coolwarm')
  ax.imshow(np.rot90(f), cmap='coolwarm', extent=[xmin, xmax, ymin, ymax])
  cset = ax.contour(xx, yy, f, colors='k')
  ax.clabel(cset, inline=1, fontsize=12)
  ax.set_xlabel('1' + '$^{st}$' + ' dimension of MDS', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of MDS', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14.3)


  figure_name = "Patients TCGA applying PCA first to 25 components followed by MDS with 2 dimensions - gaussian density estimation"
 
  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")

  # save figure
  figure_name = '3 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()


  return xx, yy, f
 



def patientsPca25MDS3D(tumorLegend):

  # Using PCA first to 25 dimensions, followed by MDS (3D)
  # Multidimensional scaling

  pca = PCA(n_components = 25)         
  patientsPCA = pca.fit_transform(inputPatientsTCGA)

  mds = MDS(n_components = 3)         
  patientsMDS = mds.fit_transform(patientsPCA)


  patientsMdsDF = pd.DataFrame(data = patientsMDS, columns = ['MDS_1', 'MDS_2', 'MDS_3'])

  #target value
  patientsMdsDF['TUMOR_TYPE'] = patientsTumorType

  patientsMdsDF = patientsMdsDF.sort_values(['TUMOR_TYPE'])


  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,10)


  list(map(lambda tumorType, legend: ax.scatter(patientsMdsDF[patientsMdsDF['TUMOR_TYPE'] == tumorType]['MDS_1'], 
                                                patientsMdsDF[patientsMdsDF['TUMOR_TYPE'] == tumorType]['MDS_2'], 
                                                patientsMdsDF[patientsMdsDF['TUMOR_TYPE'] == tumorType]['MDS_3'], 
                                                s = 60.0, c = legend.get_facecolor()), 
                                                list(patientsMdsDF['TUMOR_TYPE'].unique()), tumorLegend.legendHandles))

  ax.set_xlabel('1' + '$^{st}$' + ' dimension of MDS', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of MDS', fontsize = 14.8, labelpad = 14.5)
  ax.set_zlabel('3' + '$^{rd}$' + ' dimension of MDS', fontsize = 14.8, labelpad = 14.5)

  ax.tick_params(labelsize = 14.3)
  ax.tick_params(labelsize = 14.3)
  ax.tick_params(labelsize = 14.3)


  group_legend = plt.gca().legend(bbox_to_anchor=(1.01, 0.5), loc='center left', borderaxespad=1.0, fontsize = 13.8,  \
      labels = list(patientsMdsDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14.7, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.6, handletextpad = 0.2)    
                  

  figure_name = "Patients TCGA applying PCA first to 25 components followed by MDS with 3 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.90).set_weight("bold")


  # save figure
  figure_name = '4 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')  


  plt.show()




def patientsPca25MdsKernel3D(xx, yy, f):
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,10)

  surf = ax.plot_surface(xx, yy, f, rstride=1, cstride=1, cmap='coolwarm', edgecolor='none')
  ax.set_xlabel('1' + '$^{st}$' + ' dimension of MDS', fontsize = 14.5, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of MDS', fontsize = 14.5, labelpad = 14.5)
  ax.set_zlabel('Probability density function', fontsize = 14.5, labelpad = 28)


  ax.tick_params(axis = 'x', labelsize = 14)
  ax.tick_params(axis = 'y', labelsize = 14) #, labelleft = False, labelright = True, pad = 13.0)
  ax.tick_params(axis = 'z', labelsize = 14, pad = 13.0)


  cbar = fig.colorbar(surf, shrink=0.5, aspect=5) # add color bar indicating the PDF

  cbar.ax.tick_params(labelsize=13)

  # ax.view_init(elev=21)


  figure_name = "Patients TCGA applying PCA first to 25 components followed by MDS with 3 dimensions - 2D probability density distribution"

  plt.suptitle(figure_name, fontsize = 15, y = 0.90).set_weight("bold")

  # save figure
  figure_name = '5 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def patientsPca25Isomap2D():

  # Using PCA first to 25 dimensions, followed by Isomap (2D)
  # Non-linear dimensionality reduction through Isometric Mapping

  pca = PCA(n_components = 25)         
  patientsPCA = pca.fit_transform(inputPatientsTCGA)

  isomap = Isomap(n_components= 2)    
  patientsIsomap = isomap.fit_transform(patientsPCA)    


  patientsIsomapDF = pd.DataFrame(data = patientsIsomap, columns = ['ISOMAP_1', 'ISOMAP_2'])


  #target value
  patientsIsomapDF['TUMOR_TYPE'] = patientsTumorType

  patientsIsomapDF = patientsIsomapDF.sort_values(['TUMOR_TYPE'])


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  g = sns.scatterplot(data = patientsIsomapDF, x = 'ISOMAP_1', y = 'ISOMAP_2', hue = 'TUMOR_TYPE', ax = ax, s = 90.0)


  ax.set_xlabel('1' + '$^{st}$' + ' Isomap manifold dimension', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' Isomap manifold dimension', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14.4)
      
  legend_handles, _= g.get_legend_handles_labels()      

  tumor_legend = plt.gca().legend(bbox_to_anchor=(1.07, 0.5), loc='center left', borderaxespad=1.0, fontsize = 13.6, handles = legend_handles, \
      labels = list(patientsIsomapDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14.3, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.2)    
                  

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([55.0]), tumor_legend.legendHandles))

  figure_name = "Patients TCGA applying PCA first to 25 components followed by Isomap with 2 dimensions"


  plt.suptitle(figure_name, fontsize = 15, y = 0.94).set_weight("bold")

  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()


  return tumor_legend, patientsIsomap




def patientsPca25IsomapDensity2D(patientsIsomap):

  # Calculate the density of points, normalize it and encode it in the alpha channel of a colormap
  # Low density points are emphasized by darker color, while the high density points with the magma colormap

  patientsIsomap = patientsIsomap.T

  kernel = stats.gaussian_kde(patientsIsomap)
  weights = kernel(patientsIsomap)
  weights = weights/weights.max()


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  g = sns.scatterplot(x = patientsIsomap[0], y = patientsIsomap[1], hue = weights, palette="magma", ax = ax, s = 85.0)


  ax.set_xlabel('1' + '$^{st}$' + ' Isomap manifold dimension', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' Isomap manifold dimension', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14.3)

  
  g.get_legend().remove()

  sm = plt.cm.ScalarMappable(cmap="magma")
  sm.set_array([])

  cbar = g.figure.colorbar(sm, ticks = [0, 1])
  
  cbar.set_ticklabels(["Low", "High"])

  cbar.set_label(label='Density', size=14, weight='bold', rotation = 360, labelpad = -1, y = 0.5)
  cbar.ax.tick_params(labelsize=13.5)


  figure_name = "Patients TCGA applying PCA first to 25 components followed by Isomap with 2 dimensions - density zones"

  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")


  # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def patientsPca25IsomapKernel2D(patientsIsomap):

  x = patientsIsomap[:, 0]
  y = patientsIsomap[:, 1]

  # Manually adjusted
  xmin, xmax = -10, 10  
  ymin, ymax = -17, 8

  # Create meshgrid
  xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]


  positions = np.vstack([xx.ravel(), yy.ravel()])
  values = np.vstack([x, y])
  kernel = stats.gaussian_kde(values)
  f = np.reshape(kernel(positions).T, xx.shape)


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)
  
  ax.set_xlim(xmin, xmax)
  ax.set_ylim(ymin, ymax)
  cfset = ax.contourf(xx, yy, f, cmap='coolwarm')
  ax.imshow(np.rot90(f), cmap='coolwarm', extent=[xmin, xmax, ymin, ymax])
  cset = ax.contour(xx, yy, f, colors='k')
  ax.clabel(cset, inline=1, fontsize=12)
  ax.set_xlabel('1' + '$^{st}$' + ' Isomap manifold dimension', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' Isomap manifold dimension', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14.05)


  figure_name = "Patients TCGA applying PCA first to 25 components followed by Isomap with 2 dimensions - gaussian density estimation"
 
  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")

  # save figure
  figure_name = '3 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()


  return xx, yy, f




def patientsPca25Isomap3D(tumorLegend):

  # Using PCA first to 25 dimensions, followed by Isomap (3D)
  # Non-linear dimensionality reduction through Isometric Mapping

  pca = PCA(n_components = 25)         
  patientsPCA = pca.fit_transform(inputPatientsTCGA)

  isomap = Isomap(n_components= 3)    
  patientsIsomap = isomap.fit_transform(patientsPCA)    


  patientsIsomapDF = pd.DataFrame(data = patientsIsomap, columns = ['ISOMAP_1', 'ISOMAP_2', 'ISOMAP_3'])


  #target value
  patientsIsomapDF['TUMOR_TYPE'] = patientsTumorType

  patientsIsomapDF = patientsIsomapDF.sort_values(['TUMOR_TYPE'])


  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,10)


  list(map(lambda tumorType, legend: ax.scatter(patientsIsomapDF[patientsIsomapDF['TUMOR_TYPE'] == tumorType]['ISOMAP_1'], 
                                                patientsIsomapDF[patientsIsomapDF['TUMOR_TYPE'] == tumorType]['ISOMAP_2'], 
                                                patientsIsomapDF[patientsIsomapDF['TUMOR_TYPE'] == tumorType]['ISOMAP_3'], 
                                                s = 60.0, c = legend.get_facecolor()), 
                                                list(patientsIsomapDF['TUMOR_TYPE'].unique()), tumorLegend.legendHandles))

  ax.set_xlabel('1' + '$^{st}$' + ' Isomap manifold dimension', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' Isomap manifold dimension', fontsize = 14.8, labelpad = 14.5)
  ax.set_zlabel('3' + '$^{rd}$' + ' Isomap manifold dimension', fontsize = 14.8, labelpad = 14.5)

  ax.tick_params(labelsize = 14.3)
  ax.tick_params(labelsize = 14.3)
  ax.tick_params(labelsize = 14.3)


  group_legend = plt.gca().legend(bbox_to_anchor=(1.015, 0.6), loc='center left', borderaxespad=1.0, fontsize = 13.9,  \
      labels = list(patientsIsomapDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14.5, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.6, handletextpad = 0.2)    
                  

  figure_name = "Patients TCGA applying PCA first to 25 components followed by Isomap with 3 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.90).set_weight("bold")


  # save figure
  figure_name = '4 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def patientsPca25IsomapKernel3D(xx, yy, f):
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,10)

  surf = ax.plot_surface(xx, yy, f, rstride=1, cstride=1, cmap='coolwarm', edgecolor='none')
  ax.set_xlabel('1' + '$^{st}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 14.5)
  ax.set_zlabel('Probability density function', fontsize = 14.5, labelpad = 28)


  ax.tick_params(axis = 'x', labelsize = 14)
  ax.tick_params(axis = 'y', labelsize = 14) #, labelleft = False, labelright = True, pad = 13.0)
  ax.tick_params(axis = 'z', labelsize = 14, pad = 13.0)


  cbar = fig.colorbar(surf, shrink=0.5, aspect=5) # add color bar indicating the PDF

  cbar.ax.tick_params(labelsize=13)

  # ax.view_init(elev=21)


  figure_name = "Patients TCGA applying PCA first to 25 components followed by Isomap with 3 dimensions - 2D probability density distribution"

  plt.suptitle(figure_name, fontsize = 15, y = 0.88).set_weight("bold")


  # save figure
  figure_name = '5 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()



  
#--------------------------------------------------------------------------------------------------




#### Mice and Patients TCGA data dimensionality reduction methods  #####


def micePatientsPCA2D():
  pca = PCA(n_components = 2)         
  micePatientsPCA = pca.fit_transform(inputMicePatientsTCGA)


  micePCA = micePatientsPCA[:16, :]
  patientsPCA = micePatientsPCA[16:, :]


  micePcaDF = pd.DataFrame(data = micePCA, columns = ['PC_1', 'PC_2'])

  patientsPcaDF = pd.DataFrame(data = patientsPCA, columns = ['PC_1', 'PC_2'])



  #target value
  micePcaDF['TUMOR_TYPE'] = miceBothTumorType

  #target value
  patientsPcaDF['TUMOR_TYPE'] = patientsBothTumorType


  patientsPcaDF = patientsPcaDF.sort_values(['TUMOR_TYPE'])


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  g = sns.scatterplot(data = patientsPcaDF, x = 'PC_1', y = 'PC_2', hue = 'TUMOR_TYPE', ax = ax, s = 112.0, alpha = 0.5)


  g1 = sns.scatterplot(data = micePcaDF, x = 'PC_1', y = 'PC_2', hue = 'TUMOR_TYPE', palette = ['b', 'r'], ax = ax, s = 112.0)



  legend_handles_g1, _= g1.get_legend_handles_labels()  


  group_legend = plt.gca().legend(loc='best', borderaxespad=1.0, fontsize = 14.3, handles = legend_handles_g1, \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14.8, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)    
                  


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([55.0]), group_legend.legendHandles))


  plt.gca().add_artist(group_legend)



  ax.set_xlabel('1' + '$^{st}$' + ' principal component', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' principal component', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14.3)
      

                      
  legend_handles, _= g.get_legend_handles_labels()        


  tumor_legend = plt.gca().legend(bbox_to_anchor=(1.03, 0.5), loc='center left', borderaxespad=1.0, fontsize = 14, handles = legend_handles, \
      labels = list(patientsPcaDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14.6, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.2)    

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([55.0]), tumor_legend.legendHandles))

  # setting legend alpha 
  list(map(lambda legendHandle: legendHandle.set_alpha(0.5), tumor_legend.legendHandles))



  figure_name = "Mice and patients TCGA samples applying PCA with 2 principal components"

  plt.suptitle(figure_name, fontsize = 15.3, y = 0.94).set_weight("bold")


  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()


  return legend_handles, tumor_legend, micePatientsPCA




def micePatientsPcaDensity2D(micePatientsPCA):
  # Calculate the density of points, normalize it and encode it in the alpha channel of a colormap
  # Low density points are emphasized by darker color, while the high density points with the magma colormap

  micePatientsPCA = micePatientsPCA.T

  kernel = stats.gaussian_kde(micePatientsPCA)
  weights = kernel(micePatientsPCA)
  weights = weights/weights.max()


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  g = sns.scatterplot(x = micePatientsPCA[0], y = micePatientsPCA[1], hue = weights, palette="magma", ax = ax, s = 85.0)


  ax.set_xlabel('1' + '$^{st}$' + ' principal component', fontsize = 14.8, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' principal component', fontsize = 14.8, labelpad = 15)

  # default of axis is both
  ax.tick_params(labelsize = 14.4)

  
  g.get_legend().remove()

  sm = plt.cm.ScalarMappable(cmap="magma")
  sm.set_array([])

  cbar = g.figure.colorbar(sm, ticks = [0, 1])
  
  cbar.set_ticklabels(["Low", "High"])

  cbar.set_label(label='Density', size=14, weight='bold', rotation = 360, labelpad = -1, y = 0.5)
  cbar.ax.tick_params(labelsize=13.5)


  figure_name = "Mice and patients TCGA samples applying PCA with 2 principal components - density zones"

  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")


  # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def micePatientsPcaKernel2D(micePatientsPCA):

  x = micePatientsPCA[:, 0]
  y = micePatientsPCA[:, 1]


  # Manually adjusted
  xmin, xmax = -8, 3
  ymin, ymax = -3, 6

  # Create meshgrid
  xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]


  positions = np.vstack([xx.ravel(), yy.ravel()])
  values = np.vstack([x, y])
  kernel = stats.gaussian_kde(values)
  f = np.reshape(kernel(positions).T, xx.shape)


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)
  
  ax.set_xlim(xmin, xmax)
  ax.set_ylim(ymin, ymax)
  cfset = ax.contourf(xx, yy, f, cmap='coolwarm')
  ax.imshow(np.rot90(f), cmap='coolwarm', extent=[xmin, xmax, ymin, ymax])
  cset = ax.contour(xx, yy, f, colors='k')
  ax.clabel(cset, inline=1, fontsize=12)
  ax.set_xlabel('1' + '$^{st}$' + ' principal component', fontsize = 14.8, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' principal component', fontsize = 14.8, labelpad = 15)

  # default of axis is both
  ax.tick_params(labelsize = 14)


  figure_name = "Mice and patients TCGA samples applying PCA with 2 principal components - gaussian density estimation"


  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")

  # save figure
  figure_name = '3 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()


  return xx, yy, f




def micePatientsPCA3D(legendHandles, tumorLegend):
  pca = PCA(n_components = 3)         
  micePatientsPCA = pca.fit_transform(inputMicePatientsTCGA)


  micePCA = micePatientsPCA[:16, :]
  patientsPCA = micePatientsPCA[16:, :]


  micePcaDF = pd.DataFrame(data = micePCA, columns = ['PC_1', 'PC_2', 'PC_3'])

  patientsPcaDF = pd.DataFrame(data = patientsPCA, columns = ['PC_1', 'PC_2', 'PC_3'])



  #target value
  micePcaDF['TUMOR_TYPE'] = miceBothTumorType

  #target value
  patientsPcaDF['TUMOR_TYPE'] = patientsBothTumorType


  patientsPcaDF = patientsPcaDF.sort_values(['TUMOR_TYPE'])



  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,12)


  list(map(lambda tumorType, legend: ax.scatter(patientsPcaDF[patientsPcaDF['TUMOR_TYPE'] == tumorType]['PC_1'], 
                                                patientsPcaDF[patientsPcaDF['TUMOR_TYPE'] == tumorType]['PC_2'], 
                                                patientsPcaDF[patientsPcaDF['TUMOR_TYPE'] == tumorType]['PC_3'], 
                                                s = 60.0, c = legend.get_facecolor(), alpha = 0.5, zorder = 0.001),
                                                list(patientsPcaDF['TUMOR_TYPE'].unique()), tumorLegend.legendHandles))


  kwargs={'markersize': np.sqrt(60.0)}

  # to plot the mice points in front of the tumor points
  list(map(lambda markerColor, tumorType: ax.plot(micePcaDF[micePcaDF['TUMOR_TYPE'] == tumorType]['PC_1'], 
                                                  micePcaDF[micePcaDF['TUMOR_TYPE'] == tumorType]['PC_2'], 
                                                  micePcaDF[micePcaDF['TUMOR_TYPE'] == tumorType]['PC_3'],
                                                  markerColor, zorder = 1000, **kwargs),
                                                  ['ro', 'bo'], list(micePcaDF['TUMOR_TYPE'].unique())[::-1]))



  ax.set_xlabel('1' + '$^{st}$' + ' principal component', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' principal component', fontsize = 14.5, labelpad = 15)
  ax.set_zlabel('3' + '$^{rd}$' + ' principal component', fontsize = 14.5, labelpad = 15)

  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)


  ax.view_init(60, 270)



  group_legend = plt.gca().legend(bbox_to_anchor=(0.85, -0.08), loc='lower right', borderaxespad=1.0, fontsize = 14,   \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)     
                  

  plt.gca().add_artist(group_legend)



  tumor_legend = plt.gca().legend(bbox_to_anchor=(0.88, 0.528), loc='center left', borderaxespad=1.0, fontsize = 13.4, handles = legendHandles, \
      labels = list(patientsPcaDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.7, handletextpad = 0.2)    


  # setting legend alpha 
  list(map(lambda legendHandle: legendHandle.set_alpha(0.5), tumor_legend.legendHandles))



  figure_name = "Mice and patients TCGA samples applying PCA with 3 principal components"

  plt.suptitle(figure_name, fontsize = 15, y = 0.87).set_weight("bold")


  # save figure
  figure_name = '4 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def micePatientsPcaKernel3D(xx, yy, f):

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,10)

  surf = ax.plot_surface(xx, yy, f, rstride=1, cstride=1, cmap='coolwarm', edgecolor='none')
  ax.set_xlabel('1' + '$^{st}$' + ' principal component', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' principal component', fontsize = 14.5, labelpad = 15)
  ax.set_zlabel('Probability density function', fontsize = 14.5, labelpad = 28)
  

  ax.tick_params(axis = 'x', labelsize = 14)
  ax.tick_params(axis = 'y', labelsize = 14)
  ax.tick_params(axis = 'z', labelsize = 14, pad = 13.0)


  cbar = fig.colorbar(surf, shrink=0.5, aspect=5) # add color bar indicating the PDF

  cbar.ax.tick_params(labelsize=13)


  ax.view_init(10, 120)


  figure_name = "Mice and patients TCGA samples applying PCA with 3 principal components - 2D probability density distribution"
 
  plt.suptitle(figure_name, fontsize = 15, y = 0.83).set_weight("bold")

  # save figure
  figure_name = '5 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def micePatientsTruncatedSVD2D():

  # Dimensionality reduction using truncated SVD

  svd = TruncatedSVD(n_components = 2)         
  micePatientsSVD = svd.fit_transform(inputMicePatientsTCGA)

  miceSVD = micePatientsSVD[:16, :]
  patientsSVD = micePatientsSVD[16:, :]


  miceSvdDF = pd.DataFrame(data = miceSVD, columns = ['SVD_1', 'SVD_2'])

  patientsSvdDF = pd.DataFrame(data = patientsSVD,  columns = ['SVD_1', 'SVD_2'])



  #target value
  miceSvdDF['TUMOR_TYPE'] = miceBothTumorType

  #target value
  patientsSvdDF['TUMOR_TYPE'] = patientsBothTumorType


  patientsSvdDF = patientsSvdDF.sort_values(['TUMOR_TYPE'])



  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  g = sns.scatterplot(data = patientsSvdDF, x = 'SVD_1', y = 'SVD_2', hue = 'TUMOR_TYPE', ax = ax, s = 112.0, alpha = 0.5)


  g1 = sns.scatterplot(data = miceSvdDF, x = 'SVD_1', y = 'SVD_2', hue = 'TUMOR_TYPE', palette = ['b', 'r'], ax = ax, s = 112.0)



  legend_handles_g1, _= g1.get_legend_handles_labels()  


  group_legend = plt.gca().legend(loc='best', borderaxespad=1.0, fontsize = 14.3, handles = legend_handles_g1, \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14.6, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)    
                  


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([55.0]), group_legend.legendHandles))


  plt.gca().add_artist(group_legend)



  ax.set_xlabel('1' + '$^{st}$' + ' truncated svd component', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' truncated svd component', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14.4)
      

                      
  legend_handles, _= g.get_legend_handles_labels()        


  tumor_legend = plt.gca().legend(bbox_to_anchor=(1.06, 0.5), loc='center left', borderaxespad=1.0, fontsize = 13.8, handles = legend_handles, \
      labels = list(patientsSvdDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14.4, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.2)    

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([55.0]), tumor_legend.legendHandles))

  # setting legend alpha 
  list(map(lambda legendHandle: legendHandle.set_alpha(0.5), tumor_legend.legendHandles))



  figure_name = "Mice and patients TCGA samples applying truncated SVD with 2 principal components"

  plt.suptitle(figure_name, fontsize = 15, y = 0.94).set_weight("bold")

  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()


  return legend_handles, tumor_legend, micePatientsSVD




def micePatientsTruncatedSvdDensity2D(micePatientsSVD):
  # Calculate the density of points, normalize it and encode it in the alpha channel of a colormap
  # Low density points are emphasized by darker color, while the high density points with the magma colormap

  micePatientsSVD = micePatientsSVD.T

  kernel = stats.gaussian_kde(micePatientsSVD)
  weights = kernel(micePatientsSVD)
  weights = weights/weights.max()


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  g = sns.scatterplot(x = micePatientsSVD[0], y = micePatientsSVD[1], hue = weights, palette="magma", ax = ax, s = 85.0)


  ax.set_xlabel('1' + '$^{st}$' + ' truncated svd component', fontsize = 14.8, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' truncated svd component', fontsize = 14.8, labelpad = 15)

  # default of axis is both
  ax.tick_params(labelsize = 14.4)


  g.get_legend().remove()

  sm = plt.cm.ScalarMappable(cmap="magma")
  sm.set_array([])

  cbar = g.figure.colorbar(sm, ticks = [0, 1])
  
  cbar.set_ticklabels(["Low", "High"])

  cbar.set_label(label='Density', size=14, weight='bold', rotation = 360, labelpad = -1, y = 0.5)
  cbar.ax.tick_params(labelsize=13.5)


  figure_name = "Mice and patients TCGA samples applying truncated SVD with 2 principal components - density zones"

  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")


  # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def micePatientsTruncatedSvdKernel2D(micePatientsSVD):

  x = micePatientsSVD[:, 0]
  y = micePatientsSVD[:, 1]


  # Manually adjusted
  xmin, xmax = -4, 8
  ymin, ymax = -4, 4

  # Create meshgrid
  xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]


  positions = np.vstack([xx.ravel(), yy.ravel()])
  values = np.vstack([x, y])
  kernel = stats.gaussian_kde(values)
  f = np.reshape(kernel(positions).T, xx.shape)


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)
  
  ax.set_xlim(xmin, xmax)
  ax.set_ylim(ymin, ymax)
  cfset = ax.contourf(xx, yy, f, cmap='coolwarm')
  ax.imshow(np.rot90(f), cmap='coolwarm', extent=[xmin, xmax, ymin, ymax])
  cset = ax.contour(xx, yy, f, colors='k')
  ax.clabel(cset, inline=1, fontsize=12)
  ax.set_xlabel('1' + '$^{st}$' + ' truncated svd component', fontsize = 14.8, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' truncated svd component', fontsize = 14.8, labelpad = 15)

  # default of axis is both
  ax.tick_params(labelsize = 14)


  figure_name = "Mice and patients TCGA samples applying truncated SVD with 2 principal components - gaussian density estimation"


  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")

  # save figure
  figure_name = '3 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()


  return xx, yy, f




def micePatientsTruncatedSVD3D(legendHandles, tumorLegend):
 
  svd = TruncatedSVD(n_components = 3)         
  micePatientsSVD = svd.fit_transform(inputMicePatientsTCGA)

  miceSVD = micePatientsSVD[:16, :]
  patientsSVD = micePatientsSVD[16:, :]


  miceSvdDF = pd.DataFrame(data = miceSVD, columns = ['SVD_1', 'SVD_2', 'SVD_3'])

  patientsSvdDF = pd.DataFrame(data = patientsSVD,  columns = ['SVD_1', 'SVD_2', 'SVD_3'])



  #target value
  miceSvdDF['TUMOR_TYPE'] = miceBothTumorType

  #target value
  patientsSvdDF['TUMOR_TYPE'] = patientsBothTumorType


  patientsSvdDF = patientsSvdDF.sort_values(['TUMOR_TYPE'])


  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,12)


  list(map(lambda tumorType, legend: ax.scatter(patientsSvdDF[patientsSvdDF['TUMOR_TYPE'] == tumorType]['SVD_1'], 
                                                patientsSvdDF[patientsSvdDF['TUMOR_TYPE'] == tumorType]['SVD_2'], 
                                                patientsSvdDF[patientsSvdDF['TUMOR_TYPE'] == tumorType]['SVD_3'], 
                                                s = 60.0, c = legend.get_facecolor(), alpha = 0.5, zorder = 0.001),
                                                list(patientsSvdDF['TUMOR_TYPE'].unique()), tumorLegend.legendHandles))


  kwargs={'markersize': np.sqrt(60.0)}

  # to plot the mice points in front of the tumor points
  list(map(lambda markerColor, tumorType: ax.plot(miceSvdDF[miceSvdDF['TUMOR_TYPE'] == tumorType]['SVD_1'], 
                                                  miceSvdDF[miceSvdDF['TUMOR_TYPE'] == tumorType]['SVD_2'], 
                                                  miceSvdDF[miceSvdDF['TUMOR_TYPE'] == tumorType]['SVD_3'],
                                                  markerColor, zorder = 1000, **kwargs),
                                                  ['ro', 'bo'], list(miceSvdDF['TUMOR_TYPE'].unique())[::-1]))



  ax.set_xlabel('1' + '$^{st}$' + ' truncated svd component', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' truncated svd component', fontsize = 14.5, labelpad = 15)
  ax.set_zlabel('3' + '$^{rd}$' + ' truncated svd component', fontsize = 14.5, labelpad = 15)

  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)


  ax.view_init(60, 270)



  group_legend = plt.gca().legend(bbox_to_anchor=(0.85, -0.08), loc='lower right', borderaxespad=1.0, fontsize = 14.2,   \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14.4, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)     
                  

  plt.gca().add_artist(group_legend)



  tumor_legend = plt.gca().legend(bbox_to_anchor=(0.92, 0.528), loc='center left', borderaxespad=1.0, fontsize = 13.2, handles = legendHandles, \
      labels = list(patientsSvdDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.7, handletextpad = 0.2)    


  # setting legend alpha 
  list(map(lambda legendHandle: legendHandle.set_alpha(0.5), tumor_legend.legendHandles))



  figure_name = "Mice and patients TCGA samples applying truncated SVD with 3 principal components"

  plt.suptitle(figure_name, fontsize = 15, y = 0.89).set_weight("bold")

  # save figure
  figure_name = '4 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def micePatientsTruncatedSvdKernel3D(xx, yy, f):

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,10)

  surf = ax.plot_surface(xx, yy, f, rstride=1, cstride=1, cmap='coolwarm', edgecolor='none')
  ax.set_xlabel('1' + '$^{st}$' + ' truncated svd component', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' truncated svd component', fontsize = 14.5, labelpad = 15)
  ax.set_zlabel('Probability density function', fontsize = 14.5, labelpad = 28)
  

  ax.tick_params(axis = 'x', labelsize = 14)
  ax.tick_params(axis = 'y', labelsize = 14)
  ax.tick_params(axis = 'z', labelsize = 14, pad = 13.0)


  cbar = fig.colorbar(surf, shrink=0.5, aspect=5) # add color bar indicating the PDF

  cbar.ax.tick_params(labelsize=13)


  ax.view_init(10, 120)


  figure_name = "Mice and patients TCGA samples applying truncated SVD with 3 principal components - 2D probability density distribution"
 
  plt.suptitle(figure_name, fontsize = 15, y = 0.83).set_weight("bold")

  # save figure
  figure_name = '5 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def micePatientsPca25TSNE2D():

  # Using PCA first, followed by tSNE (2D) 25 dimensions
  # T-distributed Stochastic Neighbor Embedding

  pca = PCA(n_components = 25)         
  micePatientsPCA = pca.fit_transform(inputMicePatientsTCGA)

  tsne = TSNE(n_components = 2, perplexity=40, n_iter=300)         
  micePatientsTSNE = tsne.fit_transform(micePatientsPCA)



  miceTSNE = micePatientsTSNE[:16, :]
  patientsTSNE = micePatientsTSNE[16:, :]


  miceTsneDF = pd.DataFrame(data = miceTSNE, columns = ['TSNE_1', 'TSNE_2'])

  patientsTsneDF = pd.DataFrame(data = patientsTSNE,  columns = ['TSNE_1', 'TSNE_2'])



  #target value
  miceTsneDF['TUMOR_TYPE'] = miceBothTumorType

  #target value
  patientsTsneDF['TUMOR_TYPE'] = patientsBothTumorType


  patientsTsneDF = patientsTsneDF.sort_values(['TUMOR_TYPE'])



  fig, ax = plt.subplots()
  fig.set_size_inches(12.0, 8.80)


  g = sns.scatterplot(data = patientsTsneDF, x = 'TSNE_1', y = 'TSNE_2', hue = 'TUMOR_TYPE', ax = ax, s = 112.0, alpha = 0.5)


  # All the mice points are in the same place, on top of one another
  miceTsneDF = miceTsneDF.sort_values(['TUMOR_TYPE'], ascending = False)

  g1 = sns.scatterplot(data = miceTsneDF, x = 'TSNE_1', y = 'TSNE_2', hue = 'TUMOR_TYPE', palette = ['r', 'b'], ax = ax, s = 112.0)



  legend_handles_g1, _= g1.get_legend_handles_labels()  


  group_legend = plt.gca().legend(loc='best', borderaxespad=1.0, fontsize = 14.2, handles = legend_handles_g1, \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14.4, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)    
                  


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([55.0]), group_legend.legendHandles))


  plt.gca().add_artist(group_legend)



  ax.set_xlabel('1' + '$^{st}$' + ' dimension of TSNE embedded space', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of TSNE embedded space', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14.3)
      

                      
  legend_handles, _= g.get_legend_handles_labels()        


  tumor_legend = plt.gca().legend(bbox_to_anchor=(1.19, 0.5), loc='center left', borderaxespad=1.0, fontsize = 13.8, handles = legend_handles, \
      labels = list(patientsTsneDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14.3, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.2)    

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([55.0]), tumor_legend.legendHandles))

  # setting legend alpha 
  list(map(lambda legendHandle: legendHandle.set_alpha(0.5), tumor_legend.legendHandles))



  figure_name = "Mice and patients TCGA samples applying PCA first to 25 components followed by TSNE with 2 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.94).set_weight("bold")


  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()


  return legend_handles, tumor_legend, micePatientsTSNE




def micePatientsPca25TsneDensity2D(micePatientsTSNE):

  # Calculate the density of points, normalize it and encode it in the alpha channel of a colormap
  # Low density points are emphasized by darker color, while the high density points with the magma colormap

  micePatientsTSNE = micePatientsTSNE.T

  kernel = stats.gaussian_kde(micePatientsTSNE)
  weights = kernel(micePatientsTSNE)
  weights = weights/weights.max()


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  g = sns.scatterplot(x = micePatientsTSNE[0], y = micePatientsTSNE[1], hue = weights, palette="magma", ax = ax, s = 85.0)


  ax.set_xlabel('1' + '$^{st}$' + ' dimension of TSNE embedded space', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of TSNE embedded space', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14)


  g.get_legend().remove()

  sm = plt.cm.ScalarMappable(cmap="magma")
  sm.set_array([])

  cbar = g.figure.colorbar(sm, ticks = [0, 1])
  
  cbar.set_ticklabels(["Low", "High"])

  cbar.set_label(label='Density', size=14, weight='bold', rotation = 360, labelpad = -1, y = 0.5)
  cbar.ax.tick_params(labelsize=13.5)


  figure_name = "Mice and patients TCGA samples applying PCA first to 25 components followed by TSNE with 2 dimensions - density zones"

  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")


  # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def micePatientsPca25TsneKernel2D(micePatientsTSNE): 

  x = micePatientsTSNE[:, 0]
  y = micePatientsTSNE[:, 1]

  # Manually adjusted
  xmin, xmax = -13, 13 
  ymin, ymax = -13, 13

  # Create meshgrid
  xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]


  positions = np.vstack([xx.ravel(), yy.ravel()])
  values = np.vstack([x, y])
  kernel = stats.gaussian_kde(values)
  f = np.reshape(kernel(positions).T, xx.shape)


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)
  
  ax.set_xlim(xmin, xmax)
  ax.set_ylim(ymin, ymax)
  cfset = ax.contourf(xx, yy, f, cmap='coolwarm')
  ax.imshow(np.rot90(f), cmap='coolwarm', extent=[xmin, xmax, ymin, ymax])
  cset = ax.contour(xx, yy, f, colors='k')
  ax.clabel(cset, inline=1, fontsize=12)
  ax.set_xlabel('1' + '$^{st}$' + ' dimension of TSNE embedded space', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of TSNE embedded space', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14.1)


  figure_name = "Mice and patients TCGA samples applying PCA first to 25 components followed by TSNE with 2 dimensions - gaussian density estimation"
 
  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")

  # save figure
  figure_name = '3 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    

  plt.show()


  return xx, yy, f




def micePatientsPca25TSNE3D(legendHandles, tumorLegend):
 
  # Using PCA first, followed by tSNE (3D) 25 dimensions
  # T-distributed Stochastic Neighbor Embedding

  pca = PCA(n_components = 25)         
  micePatientsPCA = pca.fit_transform(inputMicePatientsTCGA)

  tsne = TSNE(n_components = 3, perplexity=40, n_iter=300)         
  micePatientsTSNE = tsne.fit_transform(micePatientsPCA)



  miceTSNE = micePatientsTSNE[:16, :]
  patientsTSNE = micePatientsTSNE[16:, :]


  miceTsneDF = pd.DataFrame(data = miceTSNE, columns = ['TSNE_1', 'TSNE_2', 'TSNE_3'])

  patientsTsneDF = pd.DataFrame(data = patientsTSNE,  columns = ['TSNE_1', 'TSNE_2', 'TSNE_3'])



  #target value
  miceTsneDF['TUMOR_TYPE'] = miceBothTumorType

  #target value
  patientsTsneDF['TUMOR_TYPE'] = patientsBothTumorType


  patientsTsneDF = patientsTsneDF.sort_values(['TUMOR_TYPE'])



  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,12)


  list(map(lambda tumorType, legend: ax.scatter(patientsTsneDF[patientsTsneDF['TUMOR_TYPE'] == tumorType]['TSNE_1'], 
                                                patientsTsneDF[patientsTsneDF['TUMOR_TYPE'] == tumorType]['TSNE_2'], 
                                                patientsTsneDF[patientsTsneDF['TUMOR_TYPE'] == tumorType]['TSNE_3'], 
                                                s = 60.0, c = legend.get_facecolor(), alpha = 0.5, zorder = 0.001),
                                                list(patientsTsneDF['TUMOR_TYPE'].unique()), tumorLegend.legendHandles))


  kwargs={'markersize': np.sqrt(60.0)}

  # to plot the mice points in front of the tumor points
  list(map(lambda markerColor, tumorType: ax.plot(miceTsneDF[miceTsneDF['TUMOR_TYPE'] == tumorType]['TSNE_1'], 
                                                  miceTsneDF[miceTsneDF['TUMOR_TYPE'] == tumorType]['TSNE_2'], 
                                                  miceTsneDF[miceTsneDF['TUMOR_TYPE'] == tumorType]['TSNE_3'],
                                                  markerColor, zorder = 1000, **kwargs),
                                                  ['ro', 'bo'], list(miceTsneDF['TUMOR_TYPE'].unique())[::-1]))



  ax.set_xlabel('1' + '$^{st}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 15)
  ax.set_zlabel('3' + '$^{rd}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 15)

  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)


  # ax.view_init(60, 270)



  group_legend = plt.gca().legend(bbox_to_anchor=(1.00, -0.08), loc='lower right', borderaxespad=1.0, fontsize = 14,   \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)     
                  

  plt.gca().add_artist(group_legend)



  tumor_legend = plt.gca().legend(bbox_to_anchor=(1.02, 0.528), loc='center left', borderaxespad=1.0, fontsize = 13.5, handles = legendHandles, \
      labels = list(patientsTsneDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14.3, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.7, handletextpad = 0.2)    


  # setting legend alpha 
  list(map(lambda legendHandle: legendHandle.set_alpha(0.5), tumor_legend.legendHandles))



  figure_name = "Mice and patients TCGA samples applying PCA first to 25 components followed by TSNE with 3 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.90).set_weight("bold")


  # save figure
  figure_name = '4 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def micePatientsPca25TsneKernel3D(xx, yy, f):

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,12)

  surf = ax.plot_surface(xx, yy, f, rstride=1, cstride=1, cmap='coolwarm', edgecolor='none')
  ax.set_xlabel('1' + '$^{st}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of TSNE embedded space', fontsize = 14.5, labelpad = 14.5)
  ax.set_zlabel('Probability density function', fontsize = 14.5, labelpad = 28)


  ax.tick_params(axis = 'x', labelsize = 14)
  ax.tick_params(axis = 'y', labelsize = 14, labelleft = False, labelright = True, pad = 13.0)
  ax.tick_params(axis = 'z', labelsize = 14, pad = 13.0)


  cbar = fig.colorbar(surf, shrink=0.5, aspect=5) # add color bar indicating the PDF

  cbar.ax.tick_params(labelsize=13)

  ax.view_init(elev=21)


  figure_name = "Mice and patients TCGA samples applying PCA first to 25 components followed by TSNE with 3 dimensions - 2D probability density distribution"

  plt.suptitle(figure_name, fontsize = 15, y = 0.86).set_weight("bold")


  # save figure
  figure_name = '5 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()






def micePatientsPca25MDS2D():

  # Using PCA first to 25 dimensions, followed by MDS (2D)
  # Multidimensional scaling

  pca = PCA(n_components = 25)         
  micePatientsPCA = pca.fit_transform(inputMicePatientsTCGA)

  mds = MDS(n_components = 2)         
  micePatientsMDS = mds.fit_transform(micePatientsPCA)


  miceMDS = micePatientsMDS[:16, :]
  patientsMDS = micePatientsMDS[16:, :]


  miceMdsDF = pd.DataFrame(data = miceMDS, columns = ['MDS_1', 'MDS_2'])

  patientsMdsDF = pd.DataFrame(data = patientsMDS,  columns = ['MDS_1', 'MDS_2'])



  #target value
  miceMdsDF['TUMOR_TYPE'] = miceBothTumorType

  #target value
  patientsMdsDF['TUMOR_TYPE'] = patientsBothTumorType


  patientsMdsDF = patientsMdsDF.sort_values(['TUMOR_TYPE'])



  fig, ax = plt.subplots()
  fig.set_size_inches(12.0, 8.80)


  g = sns.scatterplot(data = patientsMdsDF, x = 'MDS_1', y = 'MDS_2', hue = 'TUMOR_TYPE', ax = ax, s = 112.0, alpha = 0.5)


  g1 = sns.scatterplot(data = miceMdsDF, x = 'MDS_1', y = 'MDS_2', hue = 'TUMOR_TYPE', palette = ['b', 'r'], ax = ax, s = 112.0)



  legend_handles_g1, _= g1.get_legend_handles_labels()  


  group_legend = plt.gca().legend(loc='best', borderaxespad=1.0, fontsize = 14.2, handles = legend_handles_g1, \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14.4, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)    
                  


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([55.0]), group_legend.legendHandles))


  plt.gca().add_artist(group_legend)



  ax.set_xlabel('1' + '$^{st}$' + ' dimension of MDS', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of MDS', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14.1)
      

                      
  legend_handles, _= g.get_legend_handles_labels()        


  tumor_legend = plt.gca().legend(bbox_to_anchor=(1.19, 0.5), loc='center left', borderaxespad=1.0, fontsize = 13.8, handles = legend_handles, \
      labels = list(patientsMdsDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14.3, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.2)    

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([55.0]), tumor_legend.legendHandles))

  # setting legend alpha 
  list(map(lambda legendHandle: legendHandle.set_alpha(0.5), tumor_legend.legendHandles))



  figure_name = "Mice and patients TCGA samples applying PCA first to 25 components followed by MDS with 2 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.94).set_weight("bold")


  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()


  return legend_handles, tumor_legend, micePatientsMDS




def micePatientsPca25MdsDensity2D(micePatientsMDS):


  # Calculate the density of points, normalize it and encode it in the alpha channel of a colormap
  # Low density points are emphasized by darker color, while the high density points with the magma colormap

  micePatientsMDS = micePatientsMDS.T

  kernel = stats.gaussian_kde(micePatientsMDS)
  weights = kernel(micePatientsMDS)
  weights = weights/weights.max()


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  g = sns.scatterplot(x = micePatientsMDS[0], y = micePatientsMDS[1], hue = weights, palette="magma", ax = ax, s = 85.0)


  ax.set_xlabel('1' + '$^{st}$' + ' dimension of MDS', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of MDS', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14)


  g.get_legend().remove()

  sm = plt.cm.ScalarMappable(cmap="magma")
  sm.set_array([])

  cbar = g.figure.colorbar(sm, ticks = [0, 1])
  
  cbar.set_ticklabels(["Low", "High"])

  cbar.set_label(label='Density', size=14, weight='bold', rotation = 360, labelpad = -1, y = 0.5)
  cbar.ax.tick_params(labelsize=13.5)


  figure_name = "Mice and patients TCGA samples applying PCA first to 25 components followed by MDS with 2 dimensions - density zones"

  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")


  # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def micePatientsPca25MdsKernel2D(micePatientsMDS):

  x = micePatientsMDS[:, 0]
  y = micePatientsMDS[:, 1]

  # Manually adjusted
  xmin, xmax = -13, 13 
  ymin, ymax = -13, 13

  # Create meshgrid
  xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]


  positions = np.vstack([xx.ravel(), yy.ravel()])
  values = np.vstack([x, y])
  kernel = stats.gaussian_kde(values)
  f = np.reshape(kernel(positions).T, xx.shape)


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)
  
  ax.set_xlim(xmin, xmax)
  ax.set_ylim(ymin, ymax)
  cfset = ax.contourf(xx, yy, f, cmap='coolwarm')
  ax.imshow(np.rot90(f), cmap='coolwarm', extent=[xmin, xmax, ymin, ymax])
  cset = ax.contour(xx, yy, f, colors='k')
  ax.clabel(cset, inline=1, fontsize=12)
  ax.set_xlabel('1' + '$^{st}$' + ' dimension of MDS', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of MDS', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14)


  figure_name = "Mice and patients TCGA samples applying PCA first to 25 components followed by MDS with 2 dimensions - gaussian density estimation"
 
  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")


  # save figure
  figure_name = '3 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()


  return xx, yy, f




def micePatientsPca25MDS3D(legendHandles, tumorLegend):

  # Using PCA first to 25 dimensions, followed by MDS (3D)
  # Multidimensional scaling

  pca = PCA(n_components = 25)         
  micePatientsPCA = pca.fit_transform(inputMicePatientsTCGA)

  mds = MDS(n_components = 3)         
  micePatientsMDS = mds.fit_transform(micePatientsPCA)


  miceMDS = micePatientsMDS[:16, :]
  patientsMDS = micePatientsMDS[16:, :]


  miceMdsDF = pd.DataFrame(data = miceMDS, columns = ['MDS_1', 'MDS_2', 'MDS_3'])

  patientsMdsDF = pd.DataFrame(data = patientsMDS,  columns = ['MDS_1', 'MDS_2', 'MDS_3'])


  #target value
  miceMdsDF['TUMOR_TYPE'] = miceBothTumorType

  #target value
  patientsMdsDF['TUMOR_TYPE'] = patientsBothTumorType


  patientsMdsDF = patientsMdsDF.sort_values(['TUMOR_TYPE'])



  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,12)


  list(map(lambda tumorType, legend: ax.scatter(patientsMdsDF[patientsMdsDF['TUMOR_TYPE'] == tumorType]['MDS_1'], 
                                                patientsMdsDF[patientsMdsDF['TUMOR_TYPE'] == tumorType]['MDS_2'], 
                                                patientsMdsDF[patientsMdsDF['TUMOR_TYPE'] == tumorType]['MDS_3'], 
                                                s = 60.0, c = legend.get_facecolor(), alpha = 0.5, zorder = 0.001),
                                                list(patientsMdsDF['TUMOR_TYPE'].unique()), tumorLegend.legendHandles))


  kwargs={'markersize': np.sqrt(60.0)}

  # to plot the mice points in front of the tumor points
  list(map(lambda markerColor, tumorType: ax.plot(miceMdsDF[miceMdsDF['TUMOR_TYPE'] == tumorType]['MDS_1'], 
                                                  miceMdsDF[miceMdsDF['TUMOR_TYPE'] == tumorType]['MDS_2'], 
                                                  miceMdsDF[miceMdsDF['TUMOR_TYPE'] == tumorType]['MDS_3'],
                                                  markerColor, zorder = 1000, **kwargs),
                                                  ['ro', 'bo'], list(miceMdsDF['TUMOR_TYPE'].unique())[::-1]))


  ax.set_xlabel('1' + '$^{st}$' + ' dimension of MDS', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of MDS', fontsize = 14.5, labelpad = 15)
  ax.set_zlabel('3' + '$^{rd}$' + ' dimension of MDS', fontsize = 14.5, labelpad = 15)

  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)


  ax.view_init(60, 270)



  group_legend = plt.gca().legend(bbox_to_anchor=(0.98, -0.08), loc='lower right', borderaxespad=1.0, fontsize = 14.2,   \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14.4, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)     
                  

  plt.gca().add_artist(group_legend)



  tumor_legend = plt.gca().legend(bbox_to_anchor=(1.00, 0.528), loc='center left', borderaxespad=1.0, fontsize = 13.5, handles = legendHandles, \
      labels = list(patientsMdsDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14.3, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.7, handletextpad = 0.2)    


  # setting legend alpha 
  list(map(lambda legendHandle: legendHandle.set_alpha(0.5), tumor_legend.legendHandles))



  figure_name = "Mice and patients TCGA samples applying PCA first to 25 components followed by MDS with 3 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.87).set_weight("bold")

  # save figure
  figure_name = '4 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()

 


def micePatientsPca25MdsKernel3D(xx, yy, f):

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,10)

  surf = ax.plot_surface(xx, yy, f, rstride=1, cstride=1, cmap='coolwarm', edgecolor='none')
  ax.set_xlabel('1' + '$^{st}$' + ' dimension of MDS', fontsize = 14.5, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' dimension of MDS', fontsize = 14.5, labelpad = 14.5)
  ax.set_zlabel('Probability density function', fontsize = 14.5, labelpad = 28)


  ax.tick_params(axis = 'x', labelsize = 14)
  ax.tick_params(axis = 'y', labelsize = 14) #, labelleft = False, labelright = True, pad = 13.0)
  ax.tick_params(axis = 'z', labelsize = 14, pad = 13.0)


  cbar = fig.colorbar(surf, shrink=0.5, aspect=5) # add color bar indicating the PDF

  cbar.ax.tick_params(labelsize=13)

  # ax.view_init(elev=21)


  figure_name = "Mice and patients TCGA samples applying PCA first to 25 components followed by MDS with 3 dimensions - 2D probability density distribution"

  plt.suptitle(figure_name, fontsize = 15, y = 0.90).set_weight("bold")


  # save figure
  figure_name = '5 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def micePatientsPca25Isomap2D():

  # Using PCA first to 25 dimensions, followed by Isomap (2D)
  # Non-linear dimensionality reduction through Isometric Mapping


  pca = PCA(n_components = 25)         
  micePatientsPCA = pca.fit_transform(inputMicePatientsTCGA)

  isomap = Isomap(n_components= 2)         
  micePatientsIsomap = isomap.fit_transform(micePatientsPCA)


  miceIsomap = micePatientsIsomap[:16, :]
  patientsIsomap = micePatientsIsomap[16:, :]


  miceIsomapDF = pd.DataFrame(data = miceIsomap, columns = ['ISOMAP_1', 'ISOMAP_2'])

  patientsIsomapDF = pd.DataFrame(data = patientsIsomap,  columns = ['ISOMAP_1', 'ISOMAP_2'])



  #target value
  miceIsomapDF['TUMOR_TYPE'] = miceBothTumorType

  #target value
  patientsIsomapDF['TUMOR_TYPE'] = patientsBothTumorType


  patientsIsomapDF = patientsIsomapDF.sort_values(['TUMOR_TYPE'])



  fig, ax = plt.subplots()
  fig.set_size_inches(12.2, 8.80)


  g = sns.scatterplot(data = patientsIsomapDF, x = 'ISOMAP_1', y = 'ISOMAP_2', hue = 'TUMOR_TYPE', ax = ax, s = 112.0, alpha = 0.5)


  # All the mice points are in the same place, on top of one another
  miceIsomapDF = miceIsomapDF.sort_values(['TUMOR_TYPE'], ascending = False)


  g1 = sns.scatterplot(data = miceIsomapDF, x = 'ISOMAP_1', y = 'ISOMAP_2', hue = 'TUMOR_TYPE', palette = ['r', 'b'], ax = ax, s = 112.0)



  legend_handles_g1, _= g1.get_legend_handles_labels()  


  group_legend = plt.gca().legend(loc='best', borderaxespad=1.0, fontsize = 14.2, handles = legend_handles_g1, \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14.4, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)    
                  


  # setting legend marker color
  list(map(lambda legendHandle, colors: legendHandle.set_color(colors), group_legend.legendHandles, ['b', 'r']))

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([55.0]), group_legend.legendHandles))


  plt.gca().add_artist(group_legend)



  ax.set_xlabel('1' + '$^{st}$' + ' Isomap manifold dimension', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' Isomap manifold dimension', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14.2)
      

                      
  legend_handles, _= g.get_legend_handles_labels()        


  tumor_legend = plt.gca().legend(bbox_to_anchor=(1.19, 0.5), loc='center left', borderaxespad=1.0, fontsize = 13.8, handles = legend_handles, \
      labels = list(patientsIsomapDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14.3, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.2)    

  # setting legend marker size
  list(map(lambda legendHandle: legendHandle.set_sizes([55.0]), tumor_legend.legendHandles))

  # setting legend alpha 
  list(map(lambda legendHandle: legendHandle.set_alpha(0.5), tumor_legend.legendHandles))



  figure_name = "Mice and patients TCGA samples applying PCA first to 25 components followed by Isomap with 2 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.94).set_weight("bold")


  # save figure
  figure_name = '1 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()


  return legend_handles, tumor_legend, micePatientsIsomap




def micePatientsPca25IsomapDensity2D(micePatientsIsomap):

  # Calculate the density of points, normalize it and encode it in the alpha channel of a colormap
  # Low density points are emphasized by darker color, while the high density points with the magma colormap

  micePatientsIsomap = micePatientsIsomap.T

  kernel = stats.gaussian_kde(micePatientsIsomap)
  weights = kernel(micePatientsIsomap)
  weights = weights/weights.max()


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)


  g = sns.scatterplot(x = micePatientsIsomap[0], y = micePatientsIsomap[1], hue = weights, palette="magma", ax = ax, s = 85.0)


  ax.set_xlabel('1' + '$^{st}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14)


  g.get_legend().remove()

  sm = plt.cm.ScalarMappable(cmap="magma")
  sm.set_array([])

  cbar = g.figure.colorbar(sm, ticks = [0, 1])
  
  cbar.set_ticklabels(["Low", "High"])

  cbar.set_label(label='Density', size=14, weight='bold', rotation = 360, labelpad = -1, y = 0.5)
  cbar.ax.tick_params(labelsize=13.5)


  figure_name = "Mice and patients TCGA samples applying PCA first to 25 components followed by Isomap with 2 dimensions - density zones"

  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")


  # save figure
  figure_name = '2 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def micePatientsPca25IsomapKernel2D(micePatientsIsomap):

  x = micePatientsIsomap[:, 0]
  y = micePatientsIsomap[:, 1]

  # Manually adjusted
  xmin, xmax = -10, 10  
  ymin, ymax = -17, 8

  # Create meshgrid
  xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]


  positions = np.vstack([xx.ravel(), yy.ravel()])
  values = np.vstack([x, y])
  kernel = stats.gaussian_kde(values)
  f = np.reshape(kernel(positions).T, xx.shape)


  fig, ax = plt.subplots()
  fig.set_size_inches(11.7, 8.27)
  
  ax.set_xlim(xmin, xmax)
  ax.set_ylim(ymin, ymax)
  cfset = ax.contourf(xx, yy, f, cmap='coolwarm')
  ax.imshow(np.rot90(f), cmap='coolwarm', extent=[xmin, xmax, ymin, ymax])
  cset = ax.contour(xx, yy, f, colors='k')
  ax.clabel(cset, inline=1, fontsize=12)
  ax.set_xlabel('1' + '$^{st}$' + ' Isomap manifold dimension', fontsize = 14.8, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' Isomap manifold dimension', fontsize = 14.8, labelpad = 14.5)

  # default of axis is both
  ax.tick_params(labelsize = 14)


  figure_name = "Mice and patients TCGA samples applying PCA first to 25 components followed by Isomap with 2 dimensions - gaussian density estimation"
 
  plt.suptitle(figure_name, fontsize = 16, y = 0.945).set_weight("bold")

  # save figure
  figure_name = '3 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()


  return xx, yy, f




def micePatientsPca25Isomap3D(legendHandles, tumorLegend):

  # Using PCA first to 25 dimensions, followed by Isomap (3D)
  # Non-linear dimensionality reduction through Isometric Mapping


  pca = PCA(n_components = 25)         
  micePatientsPCA = pca.fit_transform(inputMicePatientsTCGA)

  isomap = Isomap(n_components= 3)         
  micePatientsIsomap = isomap.fit_transform(micePatientsPCA)


  miceIsomap = micePatientsIsomap[:16, :]
  patientsIsomap = micePatientsIsomap[16:, :]


  miceIsomapDF = pd.DataFrame(data = miceIsomap, columns = ['ISOMAP_1', 'ISOMAP_2', 'ISOMAP_3'])

  patientsIsomapDF = pd.DataFrame(data = patientsIsomap,  columns = ['ISOMAP_1', 'ISOMAP_2', 'ISOMAP_3'])



  #target value
  miceIsomapDF['TUMOR_TYPE'] = miceBothTumorType

  #target value
  patientsIsomapDF['TUMOR_TYPE'] = patientsBothTumorType


  patientsIsomapDF = patientsIsomapDF.sort_values(['TUMOR_TYPE'])



  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,12)


  list(map(lambda tumorType, legend: ax.scatter(patientsIsomapDF[patientsIsomapDF['TUMOR_TYPE'] == tumorType]['ISOMAP_1'], 
                                                patientsIsomapDF[patientsIsomapDF['TUMOR_TYPE'] == tumorType]['ISOMAP_2'], 
                                                patientsIsomapDF[patientsIsomapDF['TUMOR_TYPE'] == tumorType]['ISOMAP_3'], 
                                                s = 60.0, c = legend.get_facecolor(), alpha = 0.5, zorder = 0.001),
                                                list(patientsIsomapDF['TUMOR_TYPE'].unique()), tumorLegend.legendHandles))


  kwargs={'markersize': np.sqrt(60.0)}

  # to plot the mice points in front of the tumor points
  list(map(lambda markerColor, tumorType: ax.plot(miceIsomapDF[miceIsomapDF['TUMOR_TYPE'] == tumorType]['ISOMAP_1'], 
                                                  miceIsomapDF[miceIsomapDF['TUMOR_TYPE'] == tumorType]['ISOMAP_2'], 
                                                  miceIsomapDF[miceIsomapDF['TUMOR_TYPE'] == tumorType]['ISOMAP_3'],
                                                  markerColor, zorder = 1000, **kwargs),
                                                  ['ro', 'bo'], list(miceIsomapDF['TUMOR_TYPE'].unique())[::-1]))


  ax.set_xlabel('1' + '$^{st}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 15)
  ax.set_ylabel('2' + '$^{nd}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 15)
  ax.set_zlabel('3' + '$^{rd}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 15)

  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)
  ax.tick_params(labelsize = 14)


  ax.view_init(60, 270)



  group_legend = plt.gca().legend(bbox_to_anchor=(1.00, -0.08), loc='lower right', borderaxespad=1.0, fontsize = 14,   \
      labels = ["Control", "Sunitinib"], title = r"$\bf{Experiment \ group}$", title_fontsize = 14, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.85, handletextpad = 0.01)     
                  

  plt.gca().add_artist(group_legend)



  tumor_legend = plt.gca().legend(bbox_to_anchor=(1.03, 0.528), loc='center left', borderaxespad=1.0, fontsize = 13.5, handles = legendHandles, \
      labels = list(patientsIsomapDF['TUMOR_TYPE'].unique()), title = r"$\bf{Tumor \ type}$", title_fontsize = 14.3, shadow = True, \
          facecolor = 'white', borderpad = 0.7, labelspacing = 0.7, handletextpad = 0.2)    


  # setting legend alpha 
  list(map(lambda legendHandle: legendHandle.set_alpha(0.5), tumor_legend.legendHandles))



  figure_name = "Mice and patients TCGA samples applying PCA first to 25 components followed by Isomap with 3 dimensions"

  plt.suptitle(figure_name, fontsize = 15, y = 0.90).set_weight("bold")


  # save figure
  figure_name = '4 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




def micePatientsPca25IsomapKernel3D(xx, yy, f):

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  fig.set_size_inches(16,10)

  surf = ax.plot_surface(xx, yy, f, rstride=1, cstride=1, cmap='coolwarm', edgecolor='none')
  ax.set_xlabel('1' + '$^{st}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 14.5)
  ax.set_ylabel('2' + '$^{nd}$' + ' Isomap manifold dimension', fontsize = 14.5, labelpad = 14.5)
  ax.set_zlabel('Probability density function', fontsize = 14.5, labelpad = 28)


  ax.tick_params(axis = 'x', labelsize = 14)
  ax.tick_params(axis = 'y', labelsize = 14) #, labelleft = False, labelright = True, pad = 13.0)
  ax.tick_params(axis = 'z', labelsize = 14, pad = 13.0)


  cbar = fig.colorbar(surf, shrink=0.5, aspect=5) # add color bar indicating the PDF

  cbar.ax.tick_params(labelsize=13)

  # ax.view_init(elev=21)


  figure_name = "Mice and patients TCGA samples applying PCA first to 25 components followed by Isomap with 3 dimensions - 2D probability density distribution"

  plt.suptitle(figure_name, fontsize = 15, y = 0.89).set_weight("bold")

  
  # save figure
  figure_name = '5 - ' + figure_name + '.png'
    
  plt.savefig(figure_name, bbox_inches = 'tight')    


  plt.show()




#--------------------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------------------




#####  Load the mice samples and patients TCGA dataframes needed to input to the machine learning methods  #####


# Current directory: Dados 2 º Parte/


os.chdir("./drive/My Drive/Dados 2 º Parte/InputTablesFilesML")


# Mice samples genes
miceSamplesGenesInfo = pd.read_pickle("miceSamplesGenesInfo.pkl")

# Patients TCGA genes
patientsGenesInfoTCGA = pd.read_pickle("patientsGenesInfoTCGA.pkl")

# Mice samples and patients TCGA genes together
miceSamplesPatientsTCGA = pd.read_pickle("miceSamplesPatientsTCGA.pkl")


# Defining data and target values arrays for the mice samples and patient TCGA genes information 
inputMice, inputPatientsTCGA, inputMicePatientsTCGA, miceExperimentGroup, patientsTumorType, miceBothTumorType, patientsBothTumorType = dataTargetMatrices()




#--------------------------------------------------------------------------------------------------




# Change directory to mice samples plots output
os.chdir("../miceSamplesPlotsOutput")




#####  Apply dimensionality reduction methods using the mice samples genes to input to the machine learning methods  #####



os.chdir("./1 - PCA")


# Using PCA only (2D)
micePCA2D()

# Using PCA only (3D)
micePCA3D()



os.chdir("../2 - TruncatedSVD")


# Using Truncated SVD only (2D)
miceTruncatedSVD2D()

# Using Truncated SVD only (3D)
miceTruncatedSVD3D()



os.chdir("../3 - TSNE")


# Using TSNE only (2D)
miceTSNE2D()

# Using TSNE only (3D)
miceTSNE3D()



os.chdir("../4 - PCA&TSNE")


# Using PCA to 16 dimensions followed by TSNE (2D)
micePcaTSNE2D()

# Using PCA to 16 dimensions followed by TSNE (3D)
micePcaTSNE3D()



os.chdir("../5 - MDS")


# Using MDS only (2D)
miceMDS2D()

# Using MDS only (3D)
miceMDS3D()



os.chdir("../6 - PCA&MDS")


# Using PCA to 16 dimensions followed by MDS (2D)
micePcaMDS2D()

# Using PCA to 16 dimensions followed by MDS (3D)
micePcaMDS3D()



os.chdir("../7 - Isomap")


# Using Isomap only (2D)
miceIsomap2D()

# Using Isomap only (3D)
miceIsomap3D()



os.chdir("../8 - PCA&Isomap")

# Using PCA to 16 dimensions followed by Isomap (2D)
micePcaIsomap2D()

# Using PCA to 16 dimensions followed by Isomap (3D)
micePcaIsomap3D()




#--------------------------------------------------------------------------------------------------




# Change directory to patients TCGA plots output
os.chdir("../../patientsTCGAPlotsOutput")






#####  Apply dimensionality reduction methods using patients TCGA genes to input to the machine learning methods  #####



os.chdir("./1 - PCA")


# Using PCA only (2D)
tumor_legend, patientsPCA = patientsPCA2D()

# Density zones of the PCA points
patientsPcaDensity2D(patientsPCA)

# Gaussian Kernel density estimation 2D plot of the PCA points 
xx, yy, f = patientsPcaKernel2D(patientsPCA)

# Using PCA only (3D)
patientsPCA3D(tumor_legend)

# Gaussian Kernel density estimation 3D plot of the PCA points
patientsPcaKernel3D(xx, yy, f)



os.chdir("../2 - TruncatedSVD")


# Using Truncated SVD only (2D)
tumor_legend, patientsSVD = patientsTruncatedSVD2D()

# Density zones of the truncated SVD points
patientsTruncatedSvdDensity2D(patientsSVD)

# Gaussian Kernel density estimation 2D plot of the truncated SVD points  
xx, yy, f = patientsTruncatedSvdKernel2D(patientsSVD)

# Using Truncated SVD only (3D)
patientsTruncatedSVD3D(tumor_legend)

# Gaussian Kernel density estimation 3D plot of the truncated SVD points
patientsTruncatedSvdKernel3D(xx, yy, f)




os.chdir("../3 - PCA&TSNE")


# Using PCA to 25 dimensions followed by TSNE (2D)
tumor_legend, patientsTSNE = patientsPca25TSNE2D()

# Density zones of the TSNE points
patientsPca25TsneDensity2D(patientsTSNE)

# Gaussian Kernel density estimation 2D plot of the TSNE points 
xx, yy, f = patientsPca25TsneKernel2D(patientsTSNE)

# Using PCA to 25 dimensions followed by TSNE (3D)
patientsPca25TSNE3D(tumor_legend)

# Gaussian Kernel density estimation 3D plot of the TSNE points
patientsPca25TsneKernel3D(xx, yy, f)




os.chdir("../4 - PCA&MDS")


# Using PCA to 25 dimensions followed by MDS (2D)
tumor_legend, patientsMDS = patientsPca25MDS2D()

# Density zones of the MDS points
patientsPca25MdsDensity2D(patientsMDS)

# Gaussian Kernel density estimation 2D plot of the MDS points 
xx, yy, f = patientsPca25MdsKernel2D(patientsMDS)

# Using PCA to 25 dimensions followed by MDS (3D)
patientsPca25MDS3D(tumor_legend)

# Gaussian Kernel density estimation 3D plot of the MDS points
patientsPca25MdsKernel3D(xx, yy, f)




os.chdir("../5 - PCA&Isomap")


# Using PCA to 25 dimensions followed by Isomap (2D)
tumor_legend, patientsIsomap = patientsPca25Isomap2D()

# Density zones of the isomap points
patientsPca25IsomapDensity2D(patientsIsomap)

# Gaussian Kernel density estimation 2D plot of the Isomap points 
xx, yy, f = patientsPca25IsomapKernel2D(patientsIsomap)

# Using PCA to 25 dimensions followed by Isomap (3D)
patientsPca25Isomap3D(tumor_legend)

# Gaussian Kernel density estimation 3D plot of the Isomap points
patientsPca25IsomapKernel3D(xx, yy, f)




#--------------------------------------------------------------------------------------------------




# Change directory to mice and patients TCGA samples plots output

os.chdir("../../micePatientsTCGASamplesPlotsOutput")




#####  Apply dimensionality reduction methods using mice samples and patients TCGA genes to input to the machine learning methods  #####



os.chdir("./1 - PCA")


# Using PCA only (2D)
legend_handles, tumor_legend, micePatientsPCA = micePatientsPCA2D()

# Density zones of the PCA points
micePatientsPcaDensity2D(micePatientsPCA)

# Gaussian Kernel density estimation 2D plot of the PCA points 
xx, yy, f = micePatientsPcaKernel2D(micePatientsPCA)

# Using PCA only (3D)
micePatientsPCA3D(legend_handles, tumor_legend)

# Gaussian Kernel density estimation 3D plot of the PCA points
micePatientsPcaKernel3D(xx, yy, f)




os.chdir("../2 - TruncatedSVD")


# Using Truncated SVD only (2D)
legend_handles, tumor_legend, micePatientsSVD = micePatientsTruncatedSVD2D()

# Density zones of the truncated SVD points
micePatientsTruncatedSvdDensity2D(micePatientsSVD)

# Gaussian Kernel density estimation 2D plot of the truncated SVD points  
xx, yy, f = micePatientsTruncatedSvdKernel2D(micePatientsSVD)

# Using Truncated SVD only (3D)
micePatientsTruncatedSVD3D(legend_handles, tumor_legend)

# Gaussian Kernel density estimation 3D plot of the truncated SVD points
micePatientsTruncatedSvdKernel3D(xx, yy, f)




os.chdir("../3 - PCA&TSNE")


# Using PCA to 25 dimensions followed by TSNE (2D)
legend_handles, tumor_legend, micePatientsTSNE = micePatientsPca25TSNE2D()

# Density zones of the TSNE points
micePatientsPca25TsneDensity2D(micePatientsTSNE)

# Gaussian Kernel density estimation 2D plot of the TSNE points 
xx, yy, f = micePatientsPca25TsneKernel2D(micePatientsTSNE)

# Using PCA to 25 dimensions followed by TSNE (3D)
micePatientsPca25TSNE3D(legend_handles, tumor_legend)

# Gaussian Kernel density estimation 3D plot of the TSNE points
micePatientsPca25TsneKernel3D(xx, yy, f)




os.chdir("../4 - PCA&MDS")


# Using PCA to 25 dimensions followed by MDS (2D)
legend_handles, tumor_legend, micePatientsMDS = micePatientsPca25MDS2D()

# Density zones of the MDS points
micePatientsPca25MdsDensity2D(micePatientsMDS)

# Gaussian Kernel density estimation 2D plot of the MDS points 
xx, yy, f = micePatientsPca25MdsKernel2D(micePatientsMDS)

# Using PCA to 25 dimensions followed by MDS (3D)
micePatientsPca25MDS3D(legend_handles, tumor_legend)

# Gaussian Kernel density estimation 3D plot of the MDS points
micePatientsPca25MdsKernel3D(xx, yy, f)




os.chdir("../5 - PCA&Isomap")


# Using PCA to 25 dimensions followed by Isomap (2D)
legend_handles, tumor_legend, micePatientsIsomap = micePatientsPca25Isomap2D()

# Density zones of the isomap points
micePatientsPca25IsomapDensity2D(micePatientsIsomap)

# Gaussian Kernel density estimation 2D plot of the Isomap points 
xx, yy, f = micePatientsPca25IsomapKernel2D(micePatientsIsomap)

# Using PCA to 25 dimensions followed by Isomap (3D)
micePatientsPca25Isomap3D(legend_handles, tumor_legend)

# Gaussian Kernel density estimation 3D plot of the Isomap points
micePatientsPca25IsomapKernel3D(xx, yy, f)






