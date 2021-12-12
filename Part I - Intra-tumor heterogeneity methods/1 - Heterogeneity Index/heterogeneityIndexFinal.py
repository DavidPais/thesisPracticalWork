#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 23:34:53 2021

@author: david
"""


#####  Libraries  #####


# Numpy for data management
import numpy as np

# Pandas also for data management
import pandas as pd

# Seaborn for plotting and styling
import seaborn as sns
sns.set_style("darkgrid")

# Matplotlib for additional customization
from matplotlib import pyplot as plt

# To open system files
import os

# For additional customization
import matplotlib as mat

# Reverse the matplotlib settings back to default (inclunding colors of graphs)  
mat.rc_file_defaults()

# Convert from a list of lists to a list
import itertools as it




#--------------------------------------------------------------------------------------------------




#####  Methods  #####


# Get allele reference and mutation counts

"""alleleInfo""" # contains information about alleles, incluindg the reference allele and alternative allele read counts 
                 # e.g. "0/1:15,10:0.350:6:4:0.600:546,373:3:12", in this case refCounts is 15 and mutCounts is 10

# returns: a tuple with the reference allele reads count and the mutated allele reads count    
def readCounts(alleleInfo):
    refCounts = int(alleleInfo[alleleInfo.find(':') + 1 : alleleInfo.find(',')])
    mutCounts = int(alleleInfo[alleleInfo.find(',') + 1 : alleleInfo.find(':', alleleInfo.find(','))])
    
    return (refCounts, mutCounts)


# Calculate variant allele frequencies of mutations

"""readsInfo""" # contains the pairs of refCounts and mutCounts returned by the readCounts function, e.g. (15, 10)

# returns: the variant allele frequency in percentage, e.g. 40%
def calculateVAFs(readsInfo):
    return (readsInfo[1] / (readsInfo[0] + readsInfo[1])) * 100


# Calculate Shannon's Index

"""binProb""" # probability of a mutation belonging to a certain bin

# returns: part of the shannon's index formula result
def shannonIndex(binProb):
    if binProb == 0.0:
        return 0
    
    return binProb * np.log(binProb)

    
# Calculate reads total coverage (refCounts + mutCounts)

"""readsInfo""" # contains the pairs of refCounts and mutCounts returned by the readCounts function, e.g. (15, 10)

# returns: the total reads coverage (sum of the reference count and of the mutation count reads)
def totalCoverage(readsInfo):
    return (readsInfo[0] + readsInfo[1])




#--------------------------------------------------------------------------------------------------




#####  Load SNVs and CNAs input data  #####


# Current directory: Dados/variantes_somaticas/Heterogeneity Index


# print(os.getcwd())
# can also change the python directory above in spyder accordingly


# change when needed (e.g. when changing from using all SNVs to only SNVs in non-abnormal copy number regions)
# os.chdir("../../variantes_somaticas/excel_format")

# when using all SNVs
os.chdir("../excel_format")


genomeInfo = []

for filename in sorted(os.listdir(os.getcwd())):
    with open(os.path.join(os.getcwd(), filename), 'r'):
        data = pd.read_csv(filename, delimiter=';')

        data = data[data['FILTER'] == 'PASS']
        
        # will not be included because of copy number regions identification method's error
        data = data[data['CHROM'] != "chrX"]
        
        # reset indexes to start on 0 and increment by 1
        data = data.reset_index()
        
        data = data[['CHROM', 'POS', 'TUMOR']]
        
        genomeInfo.append(data)



# when it equals 1 it's because we have to remove the SNVs in CNA regions with copy number != 2, or also with major copy 2 and minor copy 0
noSNVsWithCNAs = 0



if(noSNVsWithCNAs == 1):
    # this variable is only needed when the option of using only SNVs in regions with normal copy number (major copy number = 1 and minor copy number = 1) is used
    # when it equals one, we do not use the CNAs with major copy number = 2 and minor copy number = 0
    notAllCNAs = 0




# has to be done only after removing the mutations in copy number alteration zones if removing CNAs, if using all SNVs, can be done here
if(noSNVsWithCNAs == 0): 
    tumorInfo = [mouseGenome['TUMOR'].to_list() for mouseGenome in genomeInfo]


else: 

    os.chdir("../../CNVS/excel_format_germline")

    cnasInfo = []
    
        
    for filename in sorted(os.listdir(os.getcwd())):
        with open(os.path.join(os.getcwd(), filename), 'r'):
            data = pd.read_csv(filename, delimiter='\t')
            
            # remove chrY since is breast cancer (women only have X chromosome)
            data = data[:-1]
            
            
            if(notAllCNAs == 0):
                # for the chromosome X CNAs method errors
                data = data.drop( data[ ((data['cn1'] == 1) & (data['cn2'] == 1)) | ((data['cn1'] == 0) & (data['cn2'] == 0))].index )
            
            # CNAs with major copy number 2 and minor copy number 0 will not be used (counted as CNAs, but as an error)
            else:
                data = data.drop( data[ ((data['cn1'] == 1) & (data['cn2'] == 1)) | ((data['cn1'] == 0) & (data['cn2'] == 0)) |
                                              ((data['cn1'] == 2) & (data['cn2'] == 0)) ].index )

            # drop rows with nan values
            data = data.dropna() 
        
            # reset indexes to start on 0 and increment by 1
            data = data.reset_index(drop = True)
        
        
            data = data.rename(columns = {data.columns[0]: 'CHROM', data.columns[1]: 'START POS', data.columns[2]: 'END POS'})
                         # data.columns[8]: 'CN1', data.columns[9]: 'CN2'})
      
         
            data = data[['CHROM', 'START POS', 'END POS']] # 'CN1', 'CN2']]
           
            cnasInfo.append(data)



    # Remove SNVs that are in the regions with copy major copy number != 1 and minor copy number != 1   
    
    
    i = 0
        
    while (i < len(genomeInfo)):
        j = 0
        
        while j < len(cnasInfo[i]):
          
            # SNVs that are intersected by a CNA in the same sample
            genomeInfo[i] = genomeInfo[i].drop(genomeInfo[i][(genomeInfo[i].CHROM == cnasInfo[i].iloc[j]['CHROM']) & 
                              (cnasInfo[i].iloc[j]['START POS'] <= genomeInfo[i].POS) & 
                                  (genomeInfo[i].POS<= cnasInfo[i].iloc[j]['END POS'])].index)  
            
            j+=1
            
            
        i+=1
        
    
    # if using only mut in areas without CNAs it goes here tumorInfo
    tumorInfo = [mouseGenome['TUMOR'].to_list() for mouseGenome in genomeInfo]





#--------------------------------------------------------------------------------------------------



#####  Calculate Heterogeneity Indexes  #####



reads = [list(map(readCounts, x)) for x in tumorInfo]

vafs = [list(map(calculateVAFs, read)) for read in reads]



# take only the counts and not bin_edges (using 10 bins)
counts = [list(np.histogram(vaf, bins = 10)[0]) for vaf in vafs]
 


pi = [[countValue/len(tumorInfo[counts.index(rCount)]) for countValue in rCount] for rCount in counts]


heterogeneityIndices = [-sum(list(map(shannonIndex, mousePi))) for mousePi in pi]





if(noSNVsWithCNAs == 0): 
    # if including all SNVs
    os.chdir("../Heterogeneity Index/Images")

else: 
    # if not including SNVs in CNA regions with abnormal copy number
    os.chdir("../../variantes_somaticas/Heterogeneity Index/Images")
    



#####  Dataframe of heterogeneity indexes to be used for plotting  #####


miceHI = pd.DataFrame(index=np.arange(16), columns=np.arange(4))


# heterogeneity index value, mouse, region, group type in columns

miceHI.columns = ['TH index', 'Mouse', 'Region', 'Group']

groups = ['Ctrl', 'Sunit']

mice = ['Mouse 49', 'Mouse 55', 'Mouse 61', 'Mouse 62']

regions = ['Region 1', 'Region 2', 'Region 3', 'Region 4']


# graphic analysis of distributions of data 
c_names = []


i = 0

j = 1

k = 0


for r in miceHI.index:
    miceHI.at[r, 'TH index'] = heterogeneityIndices[i]
    
    c_names.append(' ')
    

    if(i < 8):
        miceHI.at[r, 'Group'] = groups[0]
         
    else:
        miceHI.at[r, 'Group'] = groups[1]
    
    
    if(4 * j - i == 0): 
        j+=1

    
    miceHI.at[r, 'Mouse'] = mice[j - 1]

    c_names[i] = mice[j - 1] [0 : mice[j - 1].find(' ')][0] + mice[j - 1] [mice[j - 1].find(' ') + 1 : len(mice[j - 1])]


    if(k == 4): 
        k = 0
    
    miceHI.at[r, 'Region'] = regions[k]

    c_names[i] = c_names[i] + regions[k] [0 : regions[k].find(' ')][0] +  regions[k] [regions[k].find(' ') + 1 : len(regions[k])]
    

    k+=1
    
    i+=1
    
    



#--------------------------------------------------------------------------------------------------





#####  Data analysis using all samples (here is not very important if we use all SNVs or only those in regions with copy number == 2)  #####



tumorInfoAll  = []

# copying by value and not by reference
tumorInfoAll.extend(tumorInfo)

tumorInfoAll = list(it.chain.from_iterable(tumorInfoAll))




# Data's reads total coverage analysis
samplesTotalCoverage = list(map(totalCoverage, map(readCounts, tumorInfoAll)))

tot_coverage = pd.DataFrame(samplesTotalCoverage, columns = ['Total coverage'])



# Data's variant allele frequency analysis

samplesTotalVAF = list(map(calculateVAFs, map(readCounts, tumorInfoAll)))

mut_allele_frequency = pd.DataFrame(samplesTotalVAF, columns = ['Variant Allele Frequency'])







#####  Plotting data (Genome sequencing data reads total coverage distribution using all samples)  ##### 



g = sns. displot(data = tot_coverage, x = "Total coverage", kind = "kde", cut = 0)


data_points = g.fig.get_children()[1].get_lines()[0].get_xydata()


max_xy = data_points[np.argmax(data_points[:, 1])]


maxDensity = "Total coverage with \nmax density is {max_d}".format(max_d = round(max_xy[0], 2))


g.fig.get_children()[1].annotate(maxDensity, fontsize = 14, xy=(max_xy [0], max_xy[1]),  xycoords='data',
            xytext=(1.2, 0.952), textcoords='axes fraction',
            arrowprops=dict(facecolor='black', shrink=0.05),
            horizontalalignment='right', verticalalignment='center')



figure_name = "Genome sequencing data reads total coverage distribution using all samples"



for ax in g.axes.flat:
    
    
    xlabel = ax.get_xlabel()
    
    ylabel = ax.get_ylabel()
    
                
    ax.set_xlabel(xlabel, fontsize = 14, labelpad = 8)
    ax.set_ylabel(ylabel, fontsize = 14, labelpad = 8)
    
    # default of axis is both
    ax.tick_params(labelsize = 13)




plt.suptitle(figure_name, fontsize = 15, y = 1.08).set_weight("bold")


# save figure

figure_name = figure_name + '.png'
  
plt.savefig(figure_name, bbox_inches = 'tight')    

plt.show()







#####  Plotting data (Variant Allele Frequency of mutations using all samples)  ##### 



g = sns.displot(mut_allele_frequency, x = 'Variant Allele Frequency', kde = True)



figure_name = "Variant Allele Frequency distribution using all samples"



for ax in g.axes.flat:
    
    
    xlabel = ax.get_xlabel()
    
    ylabel = ax.get_ylabel()
    
                
    ax.set_xlabel(xlabel, fontsize = 14, labelpad = 8)
    ax.set_ylabel(ylabel, fontsize = 14, labelpad = 8)
    
    # default of axis is both
    ax.tick_params(labelsize = 13)




plt.suptitle(figure_name, fontsize = 15, y = 1.08).set_weight("bold")


# save figure

figure_name = figure_name + '.png'
  
plt.savefig(figure_name, bbox_inches = 'tight')    

plt.show()





#--------------------------------------------------------------------------------------------------





#####  Data analysis by sample (here is not important if we use all SNVs or only those in regions with copy number == 2)  #####



# Data's reads total coverage analysis
tumorInfoCoverage = []

# Data's variant allele frequency analysis
tumorInfoVAF = []


# copying by value and not by reference
tumorInfoCoverage.extend(tumorInfo)

# copying by value and not by reference
tumorInfoVAF.extend(tumorInfo)


samplesInfoCoverage = []

samplesInfoVAF = []


for s in range (0, len(tumorInfo)):
    tumorInfoCoverage[s] = list(map(totalCoverage, map(readCounts, tumorInfoCoverage[s])))
    tumorInfoVAF[s] = list(map(calculateVAFs, map(readCounts, tumorInfoVAF[s])))
    
    
    for s_size in range(0, len (tumorInfo[s])):
        aux_cov = []
        aux_vaf = []
        
        aux_cov.append(tumorInfoCoverage[s][s_size])
        aux_vaf.append(tumorInfoVAF[s][s_size])
         
        aux_cov.append(c_names[s])
        aux_vaf.append(c_names[s])
        
        tumorInfoCoverage[s][s_size] = aux_cov
        tumorInfoVAF[s][s_size] = aux_vaf


    samplesInfoCoverage.append(pd.DataFrame(tumorInfoCoverage[s], columns = ['Total coverage', 'Sample']))
    samplesInfoVAF.append(pd.DataFrame(tumorInfoVAF[s], columns = ['Variant Allele Frequency', 'Sample']))
    

samplesInfoCoverage = pd.concat(samplesInfoCoverage)
samplesInfoVAF = pd.concat(samplesInfoVAF)






#####  Plotting data (Genome sequencing data reads total coverage distribution by sample)  ##### 

    
g = sns.displot(samplesInfoCoverage, x = 'Total coverage', kind = "kde", cut = 0, col = "Sample", col_wrap = 4, aspect = 1.2, \
                    height = 5.6, legend = False, facet_kws={'sharey' : False, 'sharex' : False})
    
g.fig.set_size_inches(21, 12)
    
g.fig.subplots_adjust(hspace = 1.0)
  
g.fig.subplots_adjust(wspace = 0.45) 



g.set_xlabels("Total coverage").set_titles("Sample {col_name}")  



for ax in g.fig.get_children()[1:] :
    
        
    data_points = ax.get_lines()[0].get_xydata()
    
    
    max_xy = data_points[np.argmax(data_points[:, 1])]
    
    
    
    maxDensity = "Total coverage with \nmax density is {max_d}".format(max_d = round(max_xy[0], 2))
    

    ax.annotate(maxDensity, fontsize = 16, xy=(max_xy[0], max_xy[1]),  xycoords='data',
                xytext=(1.05, 0.952), textcoords='axes fraction',
                arrowprops=dict(facecolor='black', shrink=0.05),
                horizontalalignment='right', verticalalignment='center')



for ax in g.axes.flat:
    
    
    xlabel = ax.get_xlabel()
    
    ylabel = ax.get_ylabel()
    
    
    if(xlabel == ''):
        xlabel = 'Total coverage'


    if(ylabel == ''):
        ylabel = 'Density'


                
    ax.set_xlabel(xlabel, fontsize = 15, labelpad = 8)
    ax.set_ylabel(ylabel, fontsize = 15, labelpad = 8)

    ax.set_title(ax.get_title(), fontsize = 16, y = 1.20).set_weight("bold")
    
    
    # default of axis is both
    ax.tick_params(labelsize = 14)




figure_name = "Genome sequencing data reads total coverage distribution by sample"


plt.suptitle(figure_name, fontsize = 18, y = 1.10).set_weight("bold")


# save figure

figure_name = figure_name + '.png'
  
plt.savefig(figure_name, bbox_inches = 'tight')    

plt.show()







#####  Plotting data (Variant Allele Frequency by sample)  ##### 


g = sns.displot(samplesInfoVAF, x = 'Variant Allele Frequency', kde = True, col = "Sample", col_wrap = 4, aspect = 1.4, \
                    height = 5.6, legend = False, facet_kws={'sharey' : False, 'sharex' : False})
    
g.fig.set_size_inches(21, 12)
    
g.fig.subplots_adjust(hspace = 0.8)
  
g.fig.subplots_adjust(wspace = 0.25) 



g.set_xlabels("Variant Allele Frequency").set_titles("Sample {col_name}")  




for ax in g.axes.flat:
    
    
    xlabel = ax.get_xlabel()
    
    ylabel = ax.get_ylabel()
    
    
    if(xlabel == ''):
        xlabel = 'Variant Allele Frequency'


    if(ylabel == ''):
        ylabel = 'Count'


                
    ax.set_xlabel(xlabel, fontsize = 15, labelpad = 8)
    ax.set_ylabel(ylabel, fontsize = 15, labelpad = 8)

    ax.set_title(ax.get_title(), fontsize = 16, y = 1.04).set_weight("bold")
    
    
    # default of axis is both
    ax.tick_params(labelsize = 14)



figure_name = "Variant Allele Frequency distribution by sample"


plt.suptitle(figure_name, fontsize = 18, y = 1.08).set_weight("bold")


# save figure

figure_name = figure_name + '.png'
  
plt.savefig(figure_name, bbox_inches = 'tight')    

plt.show()







#--------------------------------------------------------------------------------------------------






#####  Plotting data (Heterogeneity Index)  #####


# Graphics with all SNVs or only with SNVs not intersected by CNA regions with abnormal copy number 
    


# all SNVs will be used 
if(noSNVsWithCNAs == 0):
    figure_name = "TH index comparison between experiment groups (with all SNVs)"
     
else:    
    
    if(notAllCNAs == 0):
        figure_name = "TH index comparison between experiment groups (with not all SNVs and all CNAs)"    
    else:
        figure_name = "TH index comparison between experiment groups (with not all SNVs and not all CNAs)"
        
    
    
g = sns.catplot(x = "Mouse", y = "TH index", hue = "Group", palette=sns.color_palette("tab10", 2), \
    data = miceHI, legend = False, sharex = False, sharey = False, aspect = 1.2) 


g.set_axis_labels("Mouse", "Tumor Heterogeneity index")
    
    
for ax in g.axes.flat:
    
    xlabel = ax.get_xlabel()
    
    ylabel = ax.get_ylabel()
    
    ax.set_xlabel(xlabel, fontsize = 15, labelpad = 13)
    ax.set_ylabel(ylabel, fontsize = 15, labelpad = 13)


    # default of axis is both
    ax.tick_params(labelsize = 13)
    
        
    for mp in ax.get_children()[0:4]:

         # increase data points size
         mp.set_sizes([42,42,42,42])
    


group_legend = plt.gca().legend(bbox_to_anchor=(1.03, 0.5), loc='center left', borderaxespad=1.0, fontsize = 13, labels = ['Control', 'Sunitinib'], title = "Group", title_fontsize = 14, shadow = True, \
        facecolor = 'white', borderpad = 0.7, labelspacing = 0.6, handletextpad = 0.2)    


# same color as data points
group_legend.legendHandles[0].set_color(plt.gca().get_children()[1].get_facecolor()[0:3])
group_legend.legendHandles[1].set_color(plt.gca().get_children()[2].get_facecolor()[0:3])


# set markers' sizes
group_legend.legendHandles[0].set_sizes([10])
group_legend.legendHandles[1].set_sizes([10])


    
for legHan in group_legend.legendHandles:
    legHan.set_linewidth(4.5)
    
        
group_legend._legend_box.sep = 10

    

plt.suptitle(figure_name, fontsize = 16, y = 1.10).set_weight("bold")



# save figure

figure_name = figure_name + '.png'
  
plt.savefig(figure_name, bbox_inches = 'tight')    

plt.show()
    





