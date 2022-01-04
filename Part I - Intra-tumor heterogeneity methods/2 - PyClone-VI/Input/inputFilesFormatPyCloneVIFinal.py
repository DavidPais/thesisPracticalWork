#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 18:07:23 2021

@author: david
"""


#####  Libraries  #####


# To open system files
import os

# To access all files of a specific type without using ifs
import glob

# Pandas for data management
import pandas as pd

# To test run times
# import time

# start_time = time.time()
# print("--- %s seconds ---" % (time.time() - start_time))







#####  Constants  #####


# Change the values of the constants accordingly


# Here CNAs are always included, since PyClone-VI belongs to the category of clonal inference algorithms that need CNAs information as input

# if CNAS_ERRORS = 0, all CNAs are included and CNAs with major copy number of 2 and minor copy number of 0 are not considered sequencing errors
# if CNAS_ERRORS = 1, the CNAs with major copy number of 2 and minor copy number of 0 are considered sequencing errors and are excluded

CNAS_ERRORS = 0




#--------------------------------------------------------------------------------------------------




#####  Methods  #####


def loadCNAs(filename):
    """
    Load each mouse sample CNAs to a dataframe.

    Parameters
            filename (str): the name of each csv file with each mouse sample CNAs to be load into a dataframe.

    Returns
            cnasAttributes (dataframe): each mouse sample CNAs information.
    """

    # skipfooter = 1, removes last line of each dataframe, removes chrY since the data is about breast cancer cell line (women only have the X chromosome)
    # engine = "python" to be able to use skipfooter
    # delimiter = ';' on windows csv and '\t' on linux csv
    cnasAttributes = pd.read_csv(filename, delimiter=';',  
                                 dtype = {'CHROM': 'category', 'START_POS': 'uint32', 'END_POS': 'uint32', 'MAJOR_COPY' : pd.UInt8Dtype(), 'MINOR_COPY' : pd.UInt8Dtype()},
                                 names = ['CHROM', 'START_POS', 'END_POS', 'MAJOR_COPY', 'MINOR_COPY'], header = 0,
                                 usecols = [0, 1, 2, 8, 9], skipfooter = 1, engine = "python")

    # drop rows with nan values
    cnasAttributes = cnasAttributes.dropna() 
        
    # change data type to reduce memory usage
    # because of nan or empty values, the first data type used was UInt8Dtype()
    cnasAttributes = cnasAttributes.astype({'MAJOR_COPY' : 'uint8', 'MINOR_COPY' : 'uint8'})
    

    # eliminate CNAs sequencing errors in chromosome X
    cnasAttributes = cnasAttributes.drop(cnasAttributes[(cnasAttributes['MAJOR_COPY'] == 0) & (cnasAttributes['MINOR_COPY'] == 0)].index)
    
    # CNAs with major copy number 2 and minor copy number 0 will not be used (considered as possible sequencing errors)
    if CNAS_ERRORS == 1:
        cnasAttributes = cnasAttributes.drop(cnasAttributes[((cnasAttributes['MAJOR_COPY'] == 2) & (cnasAttributes['MINOR_COPY'] == 0))].index)
    
    
    # reset indexes to start on 0 and increment by 1
    cnasAttributes.index = range(len(cnasAttributes.index))
  

    return cnasAttributes




def getSamplesIDs():
    """
    Get the each of the mice samples IDs.

    Returns  
            samples (list): The identifiers of the mice samples
    """
    
    groups = ['Ctrl', 'Sunit']
    
    miceID = ['Mouse 49', 'Mouse 55', 'Mouse 61', 'Mouse 62']
    
    regionsID = ['Region 1', 'Region 2', 'Region 3', 'Region 4']


    # To identify each mouse sample with the format M'xx'G'y'R'z', where xx = 49, 55, 61 or 62, y = Ctrl or Sunit, and z = 1, 2, 3 or 4 
    samples = [miceID[int(i/4)][0] + miceID[int(i/4)][-2:] + groups[int(i/8)] +  regionsID[i % 4][0] + regionsID[i % 4][-1:]  for i in range(0,16)]
    
    
    return samples




def loadSNVs(filename, sampleID): 
    """
    Load each mouse sample SNVs to a dataframe.
  
    Parameters  
            filename (str): the name of each csv file with a mouse sample SNVs to be load into a dataframe.
            sampleID (str): the identifier of the SNVs sample

    Returns
            genomeInfo (dataframe): each mouse sample SNVs information.
    """

    genomeInfo = pd.read_csv(filename, delimiter=';',  dtype = {'CHROM': 'category', 'POS': 'uint32', 'ALT' : 'category', 'FILTER': 'category', 'TUMOR': str},
                             names = ['CHROM', 'POS', 'ALT', 'FILTER', 'TUMOR'], header = 0, usecols = [0, 1, 4, 6, 9])

    # SNVs in chr X will not be included because of CNAs caller methods errors (with some CNAs with major and minor copy number = 0 in chr X) 
    genomeInfo = genomeInfo[(genomeInfo['FILTER'] == 'PASS') & (genomeInfo['CHROM'] != 'chrX')]
  
    
    # PyClone-VI input columns major_cn, minor_cn, normal_cn and tumour_content (cellular fraction) default values
    genomeInfo['normal_cn'], genomeInfo['major_cn'], genomeInfo['minor_cn'], genomeInfo['tumour_content'] = [2, 1, 1, 1.0]
    
    # change data type to reduce memory usage
    genomeInfo = genomeInfo.astype ({'normal_cn' : 'uint8', 'major_cn' : 'uint8', 'minor_cn' : 'uint8', 'tumour_content' : 'float16'})        
          

    # Get the reference and alternative alleles read counts
    # e.g. "0/1:15,10:0.350:6:4:0.600:546,373:3:12", in this case refCounts is 15 and mutCounts is 10

    genomeInfo['ref_counts'], genomeInfo['alt_counts'] = \
        [genomeInfo['TUMOR'].str.extract(r':(.*?),', expand = False).astype('uint16'),
         genomeInfo['TUMOR'].str.extract(r',(.*?):', expand = False).astype('uint16')]


    # Mutation identifier for each of the SNVs (fatest way with join string for dataframes with less than 100 rows)
    genomeInfo['mutation_id'] = [':'.join(snvID) for snvID in \
                                 zip(genomeInfo['CHROM'].map(str.upper), genomeInfo['POS'].map(str), genomeInfo['ALT'].map(str))] 
    

    # Sample identifier for each of the SNVs
    genomeInfo['sample_id'] = sampleID
        
    genomeInfo['sample_id'] = genomeInfo['sample_id'].astype('category')
        
    
    # remove ALT, TUMOR and FILTER columns
    genomeInfo.drop(columns = ['ALT', 'TUMOR', 'FILTER'], axis = 1, inplace = True)
        
        
    # reset indexes to start on 0 and increment by 1
    genomeInfo.index = range(len(genomeInfo.index))
  
    
    return genomeInfo




def intersectSNVsCNAs(sampleSNVsInfo, sampleCNAsInfo):
    """
    Find which SNVs are intersected by which CNAs in each mouse sample.

    Parameters
            sampleSNVsInfo (dataframe): the SNVs of a specific sample.
            sampleCNAsInfo (dataframe): the CNAs of a specific sample.
                
    Returns
            snvsIntersectedByCNAs (dataframe): the SNVs intersected by CNAs and their respective minor and major copy numbers in each mouse sample.
    """

    # Intersect the SNVs and CNAs of a same mouse sample (how = 'left' in the case a SNV is not intersected by any CNA, all the SNVs are kept)
    snvsIntersectedByCNAs = pd.merge(sampleSNVsInfo, sampleCNAsInfo, on = ['CHROM'], how = 'left')
    
    
    # Get the SNVs that are intersected by CNAs
    snvsIntersectedByCNAs = snvsIntersectedByCNAs[((snvsIntersectedByCNAs['START_POS']) <= (snvsIntersectedByCNAs['POS'])) & \
                                                  ((snvsIntersectedByCNAs['POS']) <= (snvsIntersectedByCNAs['END_POS']))] 
                                         
    # The minor and major copy number values of each SNV are the minor and major copy number values of the CNAs that intersected them
    snvsIntersectedByCNAs['major_cn'], snvsIntersectedByCNAs['minor_cn'] = [snvsIntersectedByCNAs['MAJOR_COPY'], snvsIntersectedByCNAs['MINOR_COPY']]


    # drops CHROM, POS, START_POS, END_POS, MAJOR_COPY and MINOR_COPY columns, and order the remaining columns according to the PyClone-VI input files format
    snvsIntersectedByCNAs = snvsIntersectedByCNAs[['mutation_id', 'sample_id', 'ref_counts', 'alt_counts', 'normal_cn', 'major_cn', 'minor_cn', 'tumour_content']]


    # reset indexes to start on 0 and increment by 1
    snvsIntersectedByCNAs.index = range(len(snvsIntersectedByCNAs.index))


    return snvsIntersectedByCNAs




def getSNVsToInsertSamples(mouseSample, allSNVsIdentifiers):
    """
    Get the SNVs that are not in all samples (either using all mice samples, or the four samples of each mouse).

    Parameters
            mouseSample (dataframe): each sample snvs information in the set of mice used.
            allSNVsIdentifiers (series): the identifiers of all the SNVs that exist in the set of mice used.
  
    Returns
            allSamplesSNVs (dataframe): a dataframe with the SNVs that will be add to all samples
    """

    # keep = False removes all duplicates, so if a SNV exists in that sample, is removed and only the SNVs that are not in that sample are returned
    allSamplesSNVs = pd.DataFrame(pd.concat([allSNVsIdentifiers, mouseSample['mutation_id']], ignore_index = True).drop_duplicates(keep = False))

    
    allSamplesSNVs['sample_id'], allSamplesSNVs['ref_counts'], allSamplesSNVs['alt_counts'], \
        allSamplesSNVs['normal_cn'], allSamplesSNVs['major_cn'], allSamplesSNVs['minor_cn'], allSamplesSNVs['tumour_content'] = \
            [mouseSample['sample_id'].cat.categories[0], 0, 0, 2, 1, 1, 1.0]
            
    allSamplesSNVs = allSamplesSNVs.astype ({'mutation_id' : str, 'sample_id' : 'category', 'ref_counts' : 'uint16', 'alt_counts' : 'uint16', \
                                             'normal_cn' : 'uint8', 'major_cn' : 'uint8', 'minor_cn' : 'uint8', 'tumour_content' : 'float16'})        
            

    return allSamplesSNVs

    


def putSNVsAllSamples(miceSetSamples, snvsToPutInAllSamples):
    """
    Insert all the SNVs in all the mice set used samples.

    Parameters 
            miceSetSamples (dataframe): The SNVs that are in the mice set used samples.
            snvsToPutInAllSamples (dataframe): The SNVs that are to be inserted into all the samples of the mice set used.

    Returns
            miceSetWithAllSNVs (dataframe): All the SNVs in all samples of the mice set used.
    """

    # Insert the snvs in all samples
    miceSetWithAllSNVs = pd.concat([miceSetSamples, snvsToPutInAllSamples], ignore_index = True)

    
    # Order based on sample, chromosome, and position of each SNV
    miceSetWithAllSNVs['CHROM'], miceSetWithAllSNVs['POS'] = \
       [miceSetWithAllSNVs['mutation_id'].str.extract(r'R(.*?):', expand = False).astype('uint8'), \
        miceSetWithAllSNVs['mutation_id'].str.extract(r':(.*?):', expand = False).astype('uint32')]
           
    miceSetWithAllSNVs.sort_values(['sample_id','CHROM', 'POS'], inplace = True)

    miceSetWithAllSNVs.drop(columns = ['CHROM', 'POS'], axis = 1, inplace = True)


    # change data type to reduce memory usage
    miceSetWithAllSNVs = miceSetWithAllSNVs.astype ({'mutation_id' : str, 'sample_id' : 'category', 'ref_counts' : 'uint16', 'alt_counts' : 'uint16', \
                                                     'normal_cn' : 'uint8', 'major_cn' : 'uint8', 'minor_cn' : 'uint8', 'tumour_content' : 'float16'})        
             

    # reset indexes to start on 0 and increment by 1
    miceSetWithAllSNVs.index = range(len(miceSetWithAllSNVs.index))


    return miceSetWithAllSNVs







#--------------------------------------------------------------------------------------------------




#####  Load SNVs and CNAs input data  #####


# Current directory : Dados/CNVS


os.chdir("./excel_format_germline")


# load CNAs
cnasInfo = list(map(loadCNAs, sorted(glob.iglob("*.csv", recursive=False))))




os.chdir("../../variantes_somaticas/excel_format")


# load SNVs
snvsInfo = list(map(loadSNVs, sorted(glob.iglob("*.csv", recursive=False)), getSamplesIDs()))




# Get the minor and major copy number informations of the SNVs intersected by CNAs in each mouse sample
snvsInfo = list(map(intersectSNVsCNAs, snvsInfo, cnasInfo))







#--------------------------------------------------------------------------------------------------




#####  Get the SNVs dataframes that will be input to the PyClone-VI algorithm  #####


# With all mice samples together (given that all have the same cellular line inserted into them)
allMice = pd.concat(snvsInfo, ignore_index = True)


# Find which SNVs are not in all samples for each of the samples of the mice set used and insert them in each sample dataframe
# PyClone-VI will remove all SNVs that are not in all samples of the mice set used, so we put them with ref_counts and alt_counts = 0
# With one fix argument, all the identifiers of the SNVs that exist

snvsInAllSamples = pd.concat( \
   list(map(lambda mouseSample: getSNVsToInsertSamples(mouseSample, allMice['mutation_id'].drop_duplicates()), snvsInfo)), ignore_index = True)


# Put SNVs in all samples of the mice set used
allMice = putSNVsAllSamples(allMice, snvsInAllSamples)



   
# For the four samples of each of the mice
differentMice = [snvsInfo[mouseSamples : mouseSamples + 4] for mouseSamples in range(0, len(snvsInfo), 4)]


# For each mice four samples, get the SNVs that are not in all of the four samples of a mouse
snvsInAllSamplesMouse = [pd.concat( \
                            [getSNVsToInsertSamples(mouseSample, pd.concat(mouse)['mutation_id'].drop_duplicates()) for mouseSample in mouse]) \
                             for mouse in differentMice]


# Put SNVs in all samples of the mice set used
differentMice = [putSNVsAllSamples( \
                    pd.concat(mouseSamples), snvsToInsertMouse) for mouseSamples, snvsToInsertMouse in zip(differentMice, snvsInAllSamplesMouse)]





    
    
#####  Save the SNVs dataframes information inside files that will be input to the PyClone-VI algorithm  #####


os.chdir("../../CNVS/inputFilesPyCloneVI")


# A file for all the mice samples together

if CNAS_ERRORS == 0: # line terminator avoids empty lines between rows in the tsv file
    allMice.to_csv('allMice_snvsInfo.tsv', sep='\t', header = True, index=False, line_terminator = '\n')
    
else: # with CNAS_ERRORS = 1, change the name of the tsv file to open    
    allMice.to_csv('allMice_snvsInfoWithoutCNAsErrors.tsv', sep='\t', header = True, index=False, line_terminator = '\n')




# A file for the samples of each mouse

for mouseSamplesSNVs in differentMice:
    
    if CNAS_ERRORS == 0: 
        mouseSamplesSNVs.to_csv(mouseSamplesSNVs['sample_id'][0][0:3] + "_snvsInfo.tsv", sep='\t', header = True, index=False, line_terminator = '\n')

    else: 
        mouseSamplesSNVs.to_csv(mouseSamplesSNVs['sample_id'][0][0:3] + "_snvsInfoWithoutCNAsErrors.tsv", sep='\t', header = True, index=False, line_terminator = '\n')




# remove variables to save space
del snvsInAllSamples, snvsInAllSamplesMouse




