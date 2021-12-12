# -*- coding: utf-8 -*-

"""
Created on Sat Jul 31 17:10:04 2021

@author: david
"""

# Process the input matrices needed in order to apply the machine learning methods (mice samples and TCGA patients with copy number and somatic mutation informations)
# The goal is to access the differences between the sunitinib and control groups, to see if the treatment has any effect on the mice, also using TCGA as a manifold


#####  Libraries  #####


# Pandas for data management
import pandas as pd

# To open system files
import os

# To access all files of a specific type without using ifs
import glob

# To count the number of genes with more than one occurence
from collections import Counter

# To iterate through a list of lists
import itertools 

# To check the indexes of multiple occurences of an element in an array
import numpy as np

# To load number of columns of a file before slicing it based on its number of columns
import csv

# To test run times
# import time

# start_time = time.time()
# print("--- %s seconds ---" % (time.time() - start_time))




#--------------------------------------------------------------------------------------------------




#####  Methods  #####


# Load the genes information table with the genes chromosome, start position, end position, gene identifier and gene name

"""filename""" # The name of the file containing the gene information

# returns: The genes on the X and Y chromosomes and the genes information table
def loadGenesInfo(filename):
    genesInfo = pd.read_csv(filename, delimiter='\t', # chromosome is category when loading because of the X and Y chromosomes
        dtype = {'CHROM': 'category', 'START_POS': 'uint32', 'END_POS': 'uint32', 'GENE_ID': str, 'GENE_NAME': str},
            names = ['CHROM', 'START_POS', 'END_POS', 'GENE_ID', 'GENE_NAME'], header = 0) # using all columns
           

    # eliminate genes with more than one entry, when different gene ids/symbols correspond to the same gene name
    multipleIdsGenes = list(map(lambda geneOccurences: geneOccurences[0] if geneOccurences[1] > 1 else True, Counter(genesInfo['GENE_NAME'].tolist()).items()))

    genesInfo = genesInfo[~(genesInfo.GENE_NAME.isin(multipleIdsGenes))]


    # we will not look at X and Y chromosomes (variable needed to remove these from the TCGA CNAs and SNVs)
    genesChrXY = list(genesInfo.loc[np.in1d(genesInfo['CHROM'], ['X', 'Y']), 'GENE_NAME'].unique())
               
    # eliminate genes of chromosomes X and Y
    genesInfo = genesInfo[~(genesInfo.GENE_NAME.isin(genesChrXY))]
 
 
    # convert the chromosome number column from string to int8 to numerically order the dataframe based on it and compare with the CNAs and SNVs dataframes
    genesInfo['CHROM'] = genesInfo['CHROM'].astype('uint8')

    # sort dataframe by using different columns, first by chromosome, then by start position and finally by end position
    genesInfo.sort_values(['CHROM', 'START_POS', 'END_POS'], inplace = True)
    

    # For the gene position column based on the chromosome number, start and end positions of the gene
    genesInfo['GENE_POSITION'] = 'CHR' + genesInfo['CHROM'].astype(str) + ':' + genesInfo['START_POS'].astype(str) + ':' + genesInfo['END_POS'].astype(str)
        
    
    # rearrange columns
    genesInfo = genesInfo.iloc[:, [0, 1, 2, 5, 3, 4]]

    # reset indexes to start on 0 and increment by 1
    genesInfo.index = range(len(genesInfo.index)) 


    return genesChrXY, genesInfo




# Load the CNAs germline information and process which CNAs intersect each of the genes

"""filename, numberFiles""" # The name of the mouse file (sample) accessed in the current iteration
                            # The current file number to access the list positions corresponding to each mice, group, and region

# returns: the CNAs germline information and CNAs annotations with the genes which are intersected by CNAs
def loadMiceCNAs(filename, numberFiles):
    sampleCNAsGermline = pd.read_csv(filename, delimiter=';', # MAJOR_COPY and MINOR_COPY are pd.UInt8Dtype() because they have empty cell values, and so they can not be uint8
        dtype = {'CHROM': 'category', 'START_POS': 'uint32', 'END_POS': 'uint32', 'LOG2': 'float32', 'BAF' : 'float32', 'CI_UPPER': 'float32', 'CI_LOWER': 'float32',
                 'COPY_NUMBER': 'uint8', 'MAJOR_COPY': pd.UInt8Dtype(), 'MINOR_COPY': pd.UInt8Dtype(), 'DEPTH' : 'float32'}, \
            names = ['CHROM', 'START_POS', 'END_POS', 'LOG2', 'BAF', 'CI_UPPER', 'CI_LOWER', 'COPY_NUMBER', 'MAJOR_COPY', 'MINOR_COPY', 'DEPTH'], header = 0,
                usecols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]) # not using the probes and weight columns


    # drop any rows with nan values (removes Y chromosome) 
    # after removing the empty cells, MAJOR_COPY and MINOR_COPY can be of type uint8
    sampleCNAsGermline.dropna(axis=0, how='any', inplace=True)
    
    # change data types to reduce memory usage   
    sampleCNAsGermline['MAJOR_COPY'] = sampleCNAsGermline['MAJOR_COPY'].astype('uint8')
    sampleCNAsGermline['MINOR_COPY'] = sampleCNAsGermline['MINOR_COPY'].astype('uint8')
  
    
    # remove copy number regions columns errors with major copy 0 and minor copy 0 (removes X chromosome)  
    sampleCNAsGermline = sampleCNAsGermline[~((sampleCNAsGermline['MAJOR_COPY'] == 0) & (sampleCNAsGermline['MINOR_COPY'] == 0))]



    # Process the copy alteration ID and sample ID columns, and which CNAs intersect which genes


    # for the copy alteration identifier column (to have a unique identifier/key of the CNA)
    sampleCNAsGermline['COPY_ALTERATION_ID'] = [':'.join(cnaID) for cnaID in \
                                                zip(sampleCNAsGermline['CHROM'].map(str.upper), \
                                                        sampleCNAsGermline['START_POS'].map(str), sampleCNAsGermline['END_POS'].map(str))] 

    # for the sample identifier column
    sampleCNAsGermline['SAMPLE_ID'] = ''.join([mice[int(numberFiles/4)], groups[int(numberFiles/8)], regions[numberFiles % 4]])
    
    sampleCNAsGermline['SAMPLE_ID'] = sampleCNAsGermline['SAMPLE_ID'].astype('category')
 
 
    # convert from chr to chromosome number to merge the chromosome number with the genesInfo chromosome number column
    sampleCNAsGermline['CHROM'] = sampleCNAsGermline['CHROM'].str[3:]

    sampleCNAsGermline['CHROM'] = sampleCNAsGermline['CHROM'].astype('uint8')


    # rearrange columns
    sampleCNAsGermline = sampleCNAsGermline.iloc[:, [0, 1, 2, 11, 12, 7, 8, 9, 5, 6, 3, 4, 10]]

    # reset indexes to start on 0 and increment by 1
    sampleCNAsGermline.index = range(len(sampleCNAsGermline.index))



    # Process the CNAs intersections of the genes, which genes are intersected by which CNAs
    sampleCNAsAnnotations = pd.merge(genesInfo, sampleCNAsGermline, on = ['CHROM'], how = 'inner')
    
    # CNA_START_POS <= GENE_START_POS & GENE_END_POS <= CNA_END_POS, the gene is totally inside the CNA region
    sampleCNAsAnnotations = sampleCNAsAnnotations[ (sampleCNAsAnnotations['START_POS_y'] <= sampleCNAsAnnotations['START_POS_x']) & \
                                (sampleCNAsAnnotations['END_POS_x'] <= sampleCNAsAnnotations['END_POS_y'])]


    # rearrange columns
    sampleCNAsAnnotations = sampleCNAsAnnotations.iloc[:, [0, 1, 2, 3, 5, 9, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17]]

    sampleCNAsAnnotations.columns = ['CHROM', 'GENE_START_POS', 'GENE_END_POS', 'GENE_POSITION', 'GENE_NAME', 'SAMPLE_ID', \
                                        'CNA_START_POS', 'CNA_END_POS', 'COPY_ALTERATION_ID', 'COPY_NUMBER', 'MAJOR_COPY', 'MINOR_COPY', \
                                            'CI_UPPER', 'CI_LOWER', 'LOG2', 'BAF' , 'DEPTH']    

    # reset indexes to start on 0 and increment by 1
    sampleCNAsAnnotations.index = range(len(sampleCNAsAnnotations.index))


    # CNAs germline information, genes intersected by CNAs and the mouse sample identifier
    return sampleCNAsGermline, sampleCNAsAnnotations, sampleCNAsGermline.loc[0]['SAMPLE_ID']




# Remove genes from cnasAnnotations that do not have CNAs in every sample and reset the index numbers starting from 0 and incrementing by 1
# Get the copy number of each of the genes for each of the mice samples

"""sampleCNAsAnnotations""" # the intersection between the CNAs and the genes positions for each sample

# returns: the genes intersected by CNAs which have CNAs in all samples and the copy number of those genes alphabetically order
def getMiceGenesCopyNumber(sampleCNAsAnnotations):
    sampleCNAsAnnotations = sampleCNAsAnnotations[~(sampleCNAsAnnotations.GENE_NAME.isin(gWithoutCNAsMice))]


    # Sort the genes alphabetically to use the same order for all mice samples
    sampleCNAsAnnotations = sampleCNAsAnnotations.sort_values(['GENE_NAME'])

    # reset indexes to start on 0 and increment by 1
    sampleCNAsAnnotations.index = range(len(sampleCNAsAnnotations.index)) 


    # Get the copy number of each of the genes 
    sampleCopyNumber = sampleCNAsAnnotations['COPY_NUMBER']
    
    # Get the representative values of the copy number values of each of the genes 
    sampleCopyNumber = np.where(sampleCopyNumber == 2, 0, np.where(sampleCopyNumber > 2, 1, -1)).astype('int8').tolist()   


    # create a list for the copy number values and somatic mutation occurences
    miceGenesInfo = [''] * len(sampleCopyNumber) * 2
    
    
    # All the even positions contain the sample copy number of a gene
    miceGenesInfo[0::2] = sampleCopyNumber
    
    
    return sampleCNAsAnnotations, miceGenesInfo




# Change the value of the occurences of SNVs in a gene in a sample from 0 to 1 

"""numberFiles, mutatedGenes, snvsGenesOccurences""" # each of the mice samples
                                                     # the genes intersected by SNVs in each mice sample
                                                     # if a gene is mutated or not by a SNV
                                
# returns: no return keyword since it updates the miceSamplesGenesInfo list of samples copy number and somatic mutation genes values
def genesWithSNVs(numberFiles, mutatedGenes, snvsGenesOccurences):
    
    # Get the indexes of the genes which are mutated
    mutatedIndexes = np.array(cnasAnnotations[numberFiles].loc[ \
            ((cnasAnnotations[numberFiles]['SAMPLE_ID'].astype(str) + cnasAnnotations[numberFiles]['GENE_NAME']).isin(mutatedGenes))].index).astype('uint32')
    
    # change the occurence of a SNV in a gene to 1 in the list of somatic mutations occurences
    snvsGenesOccurences[mutatedIndexes] = 1
    
    # update the odd positions with the SNVs occurences information
    miceSamplesGenesInfo[numberFiles][1::2] = snvsGenesOccurences




# Load the genes intersected by SNVs and process which genes intersected by CNAs are intersected by which SNVs

"""filename, numberFiles, snvsGenesOccurences""" # The name of the mouse file (sample) accessed in the current iteration
                                                 # The current file number to access the list positions corresponding to each mice, group, and region, and the cnasAnnotations sample 
                                                 # The initial list of mutated genes(at the beginning none of the genes is mutated)
    
# returns: the genes intersected by SNVs and genes intersected by CNAs that are intersected by SNVs
def loadMiceSNVs(filename, numberFiles, snvsGenesOccurences):
    sampleSNVsGenes = pd.read_csv(filename, delimiter=';', \
        dtype = {'CHROM': 'category', 'POS': 'uint32', 'REF': 'category', 'ALT': 'category', 'PROTEIN_EFFECT': 'category', 'FILTER': 'category', 'TUMOR': str},                          
        names = ['GENE_NAME', 'PROTEIN_EFFECT', 'CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'TUMOR'], header = 0,
            usecols = [13, 14, 16, 17, 6, 8, 19, 22]) [['CHROM', 'POS', 'REF', 'ALT', 'GENE_NAME', 'PROTEIN_EFFECT', 'FILTER', 'TUMOR']]


    # consider only the SNVs which filter condition is PASS
    sampleSNVsGenes = sampleSNVsGenes[sampleSNVsGenes['FILTER'] == 'PASS']

    # remove the filter column
    sampleSNVsGenes.drop(columns = ['FILTER'], axis=1, inplace=True)


    # we will not look at X and Y chromosomes (SNVs files only have the X chromosome)
    sampleSNVsGenes = sampleSNVsGenes[sampleSNVsGenes['CHROM'] != "chrX"]
    
    
    # Remove the genes that have SNVs in a specific sample but do not have CNAs in every sample
    # They may or may not have SNVs, but they have to have CNAs in every sample
    
    # On the other hand, there may also be genes with SNVs that did not had any CNA in any sample so they were not detected, and we need to remove these
    
    # This applies to both situations
    sampleSNVsGenes = sampleSNVsGenes[sampleSNVsGenes.GENE_NAME.isin(gCNAsAll)]



    # Process the reference and alternative alleles read counts columns, and the mutation and sample identifiers columns
    # Process the SNVs that intersect the genes which are intersected by CNAs 


    # Reference and alterantive alleles read counts
    # e.g. "0/1:15,10:0.350:6:4:0.600:546,373:3:12", in this case refCounts is 15 and mutCounts is 10

    sampleSNVsGenes['REF_COUNTS'] = sampleSNVsGenes['TUMOR'].str.extract(r':(.*?),', expand = False).astype('uint16')
    sampleSNVsGenes['ALT_COUNTS'] = sampleSNVsGenes['TUMOR'].str.extract(r',(.*?):', expand = False).astype('uint16')
    
    # remove the tumor column
    sampleSNVsGenes.drop(columns = ['TUMOR'], axis=1, inplace=True)
    

    # for the mutation identifier column (to have a unique identifier/key of the SNV)
    sampleSNVsGenes['MUTATION_ID']  = [':'.join(snvID) for snvID in \
                                                zip(sampleSNVsGenes['CHROM'].map(str.upper), sampleSNVsGenes['POS'].map(str), \
                                                         sampleSNVsGenes['REF'], sampleSNVsGenes['ALT'])] 

    # for the sample identifier column
    sampleSNVsGenes['SAMPLE_ID'] = ''.join([mice[int(numberFiles/4)], groups[int(numberFiles/8)], regions[numberFiles % 4]])
 
    sampleSNVsGenes['SAMPLE_ID'] = sampleSNVsGenes['SAMPLE_ID'].astype('category')
    

    # convert from chr to chromosome number to merge the chromosome number with the cnasAnnotations chromosome number column
    sampleSNVsGenes['CHROM'] = sampleSNVsGenes['CHROM'].str[3:]

    sampleSNVsGenes['CHROM'] = sampleSNVsGenes['CHROM'].astype('uint8')


    # rearrange columns
    sampleSNVsGenes = sampleSNVsGenes.iloc[:, [0, 1, 2, 3, 8, 9, 6, 7, 4, 5]]

    # reset indexes to start on 0 and increment by 1
    sampleSNVsGenes.index = range(len(sampleSNVsGenes.index))



    # Intersect the positions of the CNAs that intersect genes with the positions of the SNVs to dicover which genes intersected by CNAs are intersected by SNVs
    sampleGenesAnnotations = pd.merge(cnasAnnotations[numberFiles], sampleSNVsGenes, on = ['CHROM', 'GENE_NAME'], how = 'inner')

    # CNA_START_POS <= SNV_POS & SNV_POS <= CNA_END_POS, the SNV is totally inside the CNA region
    sampleGenesAnnotations = sampleGenesAnnotations[ \
        (sampleGenesAnnotations['CNA_START_POS'] <= sampleGenesAnnotations['POS']) & (sampleGenesAnnotations['POS'] <= sampleGenesAnnotations['CNA_END_POS'])]


    # remove the second sample column resulting from the merge
    sampleGenesAnnotations.drop(columns = ['SAMPLE_ID_y'], axis=1, inplace=True)

    sampleGenesAnnotations.columns = ['CHROM', 'GENE_START_POS', 'GENE_END_POS', 'GENE_POSITION', 'GENE_NAME', 'SAMPLE_ID', 'CNA_START_POS', 'CNA_END_POS', 'COPY_ALTERATION_ID', \
        'COPY_NUMBER', 'MAJOR_COPY', 'MINOR_COPY', 'CI_UPPER', 'CI_LOWER', 'LOG2', 'BAF' , 'DEPTH', 'SNV_POS', 'REF_ALLELE', 'ALT_ALLELE', 'MUTATION_ID', \
            'REF_COUNTS', 'ALT_COUNTS', 'PROTEIN_EFFECT']    
            
    # reset indexes to start on 0 and increment by 1
    sampleGenesAnnotations.index = range(len(sampleGenesAnnotations.index))



    # The genes intersected by SNVs
    mutatedGenes = (sampleGenesAnnotations['SAMPLE_ID'].astype(str) + sampleGenesAnnotations['GENE_NAME']).unique().tolist()
    

    # update the miceSamplesGenesInfo sample list with the indication of the mutated genes
    genesWithSNVs(numberFiles, mutatedGenes, snvsGenesOccurences)


    # Genes intersected by SNVs and genes intersected by CNAs that are intersected by SNVs
    return sampleSNVsGenes, sampleGenesAnnotations




# Load the files with the representative copy number values for the CNAs that intersect the genes in each tumor

"""filename""" # The name of the tumor file accessed in the current iteration

# returns: each TCGA tumor CNAs samples copy number values for each of the genes
def focalScoresCNAsTCGA(filename):
    with open(os.path.join(os.getcwd(), filename), 'r') as f:
        # get the number of columns in the file (to drop the columns Gene ID and Cytoband)
        cNames = csv.DictReader(f, delimiter = ';').fieldnames
        nCols = len(cNames)

        
    # dtype int8 because of the rows in which genes CNAs samples have copy number < 2, with representative value of -1
    dtypesDict = dict.fromkeys(cNames[3: ], 'int8')
       
    # dtype of the Gene Symbol column
    dtypesDict[cNames[0]] = str
    
    
    sampleFocalScores = pd.read_csv(filename, sep = ';', dtype = dtypesDict, usecols = [0, *range(3, nCols)])


    # Trim the gene symbol last two numbers '.xx' of the CNA focal scores (e.g. ENSG00000008128.21 to ENSG00000008128)
    sampleFocalScores['Gene Symbol'] = sampleFocalScores['Gene Symbol'].str.extract(r'(.*)\.', expand = False)



    # merge both the gene symbol and gene identifier columns to associate the gene name to each gene symbol
    sampleFocalScores = pd.merge(sampleFocalScores, genesInfo, left_on = 'Gene Symbol', right_on = 'GENE_ID', how = 'inner')


    # rearrange columns     
    focalScoresCols = len(sampleFocalScores.columns.tolist())

    # Gene name, chromosome, start and end position, and identifier, with the rest of the columns, excluding gene symbol and gene identifier
    sampleFocalScores = sampleFocalScores.iloc[:, [ \
                focalScoresCols - 1, *range(focalScoresCols - 6, focalScoresCols - 2), *range(1, focalScoresCols - 6) ]]


    # remove genes of the chromosomes X and Y
    sampleFocalScores = sampleFocalScores[~(sampleFocalScores.GENE_NAME.isin(genesChrXY))]

 
    # eliminate genes with more than one entry, when different gene ids/symbols correspond to the same gene name
    multipleIdsGene = list(map(lambda geneOccurences: \
                            geneOccurences[0] if geneOccurences[1] > 1 else True, Counter(sampleFocalScores['GENE_NAME'].tolist()).items()))

    sampleFocalScores = sampleFocalScores[~(sampleFocalScores.GENE_NAME.isin(multipleIdsGene))]
    
    
    # Sort the genes alphabetically
    sampleFocalScores = sampleFocalScores.sort_values(['GENE_NAME'])
    

    # reset indexes to start on 0 and increment by 1
    sampleFocalScores.index = range(len(sampleFocalScores.index)) 


    # Each CNA sample copy number representative values for each of the genes
    return sampleFocalScores 




# Get the copy number representative values for each of the patient genes

"""sampleCopySegments, numberFiles""" # The mask copy number segments of each patient 
                                      # The sample number of each tumor                                                        

# returns: the representative copy number values for each patient genes based on the corresponding focal sample scores
def genesCNAsTCGA(sampleCopySegments, numberFiles):

    # All the unique combinations of patient identifier and focal sample idenfier
    patientSamplesCNAs = sampleCopySegments.groupby(['FOCAL_SAMPLE']).first()[['PATIENT_ID']]

    patientSamplesCNAs = patientSamplesCNAs.sort_values(['PATIENT_ID'])

   
    # Get only the columns which focal samples correspond to the patients focal samples identifers
    focalScores = scoresCNAsTCGA[numberFiles][patientSamplesCNAs.index.to_list()]
   
    # The rows become the focal samples identifiers, with the columns being the genes
    focalScores = focalScores.T



    # Join the patients with their corresponding focal samples
    patientCNAsScores = focalScores.merge(patientSamplesCNAs, left_index = True, right_index = True)

    patientCNAsScores = patientCNAsScores.sort_values(['PATIENT_ID'])


    cols = list(patientCNAsScores.columns)
    
    # put patient identifier in the beginning of the columns
    cols = [cols[-1]] + cols [:-1]
    
    patientCNAsScores = patientCNAsScores[cols]


    # if a patient has more than a focal sample with each focal sample having different values for the genes copy number
    if(not patientCNAsScores['PATIENT_ID'].is_unique):
    
        # the decision of the final value that will represent the CNA region copy number of a gene for a specific patient (summarizing all CNA samples)
        # using a sum of different samples values based on a bold approach (which benefits copy numbers values corresponding to abnormal copy number regions)
                
        # group by the patient focal scores samples and sum their respective copy numbers values for each of the genes
        cnasValuesDecisions = patientCNAsScores.groupby(['PATIENT_ID'])[cols].sum()
        
        # CNAs representative values decision, if the sum of all the samples of a patient is bigger than 1, it is presented by 1
        # If it is lower -1, it is represented by -1. If is equal to -1, 1, or 0 it does not change
        patientCNAsScores = list(np.where(cnasValuesDecisions > 1, 1, np.where(cnasValuesDecisions < -1, -1, cnasValuesDecisions)))
        
    else:
        
        cnasValuesDecisions = patientCNAsScores.groupby(['PATIENT_ID'])[cols].first()
        
        # do not include the patient column, only the CNAs decision values
        cnasValuesDecisions = cnasValuesDecisions[cnasValuesDecisions.columns[1:]]
        
        
        # convert each of the dataframe rows in a list, with the representative copy number values for each of the patients
        patientCNAsScores = list(cnasValuesDecisions.values)
   

    # identify each of the lists with patients copy numbers representative values for each of the genes
    patients = cnasValuesDecisions.index.values

    patientCNAsScores = dict(zip(patients, patientCNAsScores))
    
    
    return patientCNAsScores
    
    


# Load the TCGA CNAS segments information, including the CNAs focal scores samples, the CNAs chromosome, start and end position, log2 and patients CNAs samples

"""filename, numberFiles""" # The name of the TCGA file (tumor) accessed in the current iteration
                            # The current file number to access the list positions of the focal scores corresponding to each tumor

# returns: the TCGA CNAs copy number segments of each tumor and each patient CNAs samples copy number segments representative values for each gene
def copySegmentsCNAsTCGA(filename, numberFiles):
    sampleCopySegments = pd.read_csv(filename, delimiter='\t', 
        dtype = {'FOCAL_SAMPLE' : str, 'CHROM': 'category', 'CNA_START_POS': 'uint32', 'CNA_END_POS': 'uint32', \
                 'LOG2': 'float16', 'PATIENT_CNA' : 'category'}, \
            names =  ['FOCAL_SAMPLE', 'CHROM', 'CNA_START_POS', 'CNA_END_POS', 'LOG2', 'PATIENT_CNA'], header = 0,
                usecols = [0, 1, 2, 3, 5, 6]) # not using the number probes column
    

    # for dataframes with many rows, + is faster than ''.join
    sampleCopySegments['COPY_ALTERATION_ID'] = 'CHR' + sampleCopySegments['CHROM'].astype(str) + ':' + \
        sampleCopySegments['CNA_START_POS'].astype(str) + ':' + sampleCopySegments['CNA_END_POS'].astype(str)  
    

    # the type of each tumor based on the file name   
    sampleCopySegments['TUMOR_TYPE'] = \
        [filename[filename.find('-') + 1 : filename.find('_', filename.find('-'))]] * len(sampleCopySegments)
    
    # change data type to reduce memory usage
    sampleCopySegments['TUMOR_TYPE'] = sampleCopySegments['TUMOR_TYPE'].astype('category')


    # To see the CNAs samples corresponding to each of the patients
    # Also used when loading the SNVs to identify the SNV samples of each patient
    sampleCopySegments['PATIENT_ID'] = sampleCopySegments['PATIENT_CNA'].str[:12]

 
    # rearrange columns
    sampleCopySegments = sampleCopySegments.iloc[:, [8, 5, 0, 6, 1, 2, 3, 4, 7]]

    
    # we will not look at the X chromosome (TCGA CNAs segments files only have the X chromosome)
    sampleCopySegments = sampleCopySegments[sampleCopySegments['CHROM'] != "X"]



    # Find which CNAs segments samples identifiers of a patient are part of the focal scores table CNAs identifiers samples columns 
    # The focal scores table has the corresponding copy number representative values of the genes for the CNAs segments

    # All the genes in focal scores are in all the copy segments samples identifiers of each tumor


    # To obtain which mask copy number segment samples are present in the CNAs focal scores dataframe
    # The mask copy number segment samples have some focal scores samples identifiers that are not in the focal scores dataframe
    sampleCopySegments = \
         sampleCopySegments[sampleCopySegments.FOCAL_SAMPLE.isin(scoresCNAsTCGA[numberFiles].columns[5:])]


    # reset indexes to start on 0 and increment by 1
    sampleCopySegments.index = range(len(sampleCopySegments.index))



    # get the representative copy number values for each of the patients
    sampleCNAsGenesTCGA = genesCNAsTCGA(sampleCopySegments, numberFiles)


    # change data type to reduce memory usage
    sampleCopySegments['PATIENT_ID'] = sampleCopySegments['PATIENT_ID'].astype('category')
    sampleCopySegments['FOCAL_SAMPLE'] = sampleCopySegments['FOCAL_SAMPLE'].astype('category')

    
    # The TCGA CNAs copy number segments of each tumor 
    # Each patient CNAs samples copy number segments representative values for each gene
    return sampleCopySegments, sampleCNAsGenesTCGA  
    



# Update the list of patients with the indication of their mutated genes with the patients that do not have any SNVs

"""patientMutatedGenes, patientID""" # list of patients with the indication of their mutated genes
                                     # the identifier of a patient without mutated genes            

# returns: no return keyword, it updates the patientMutatedGenes dataframe with the patients that do not have mutated genes
def patientsWithoutSNVs(patientMutatedGenes, patientID):
       # the patient does not have any mutated genes, so the list of mutated genes indexes is empty
       patientMutatedGenes.loc[len(patientMutatedGenes) + 1] = [patientID, [], np.array([0] * len(gCNAsTumor), dtype = 'i1')]


    

# Update the patientsMutatedGenes dataframe with the mutated genes indication (equal to 1, a gene is mutated)

"""patientSNVsOccurences, genesIndexes""" # the mutations occurences representative values array, initialized with 0s at the beginning
                                          # indexes of the mutated genes of a patient

# returns: no return keyword, it updates the patientsMutatedGenes dataframe mutated genes with the value of 1
def identifySNVsOccurences(patientSNVsOccurences, genesIndexes):
   patientSNVsOccurences[genesIndexes] = 1    




# Identify with 1 the corresponding mutated genes 

"""numberFiles, sampleGenesSNVsPatients, genesNamesIndexes""" # The index of the patient genes table corresponding tumor sample
                                                              # The mutations of the patients genes 
                                                              # The genes names and their corresponding indexes
                                           
# returns: a dictionary with the patients identifiers as keys and the corresponding genes SNVs values, 0 if the gene is not mutated, 1 if it is
def genesSNVsTCGA(numberFiles, sampleGenesSNVsPatients, genesNamesIndexes): 
    
    snvsGenesNamesIndexes = pd.merge(sampleGenesSNVsPatients, genesNamesIndexes, on = ['GENE_NAME'], how = 'inner') \
                                .sort_values('PATIENT_ID')[['PATIENT_ID', 'GENE_INDEX']].drop_duplicates()
    
    
    # list of the mutated genes of each patient
    # observed = True since patient identifier is a category variable
    patientsMutatedGenes = snvsGenesNamesIndexes.groupby('PATIENT_ID', sort = False, observed = True)['GENE_INDEX'].agg(list).reset_index()

    # dtype is integer to be the same type of patientCNAsTCGA arrays
    patientsMutatedGenes['SNVS_OCCURENCES'] = \
        [np.array([0] * len(gCNAsTumor), dtype = 'i1') for patient in range(0, len(patientsMutatedGenes))]    
    
   
    # based on the genes mutated indexes, change the SNV value occurence from 0 to 1
    [identifySNVsOccurences(patientSNVsOccurences, genesIndexes) for patientSNVsOccurences, genesIndexes in \
             zip(patientsMutatedGenes['SNVS_OCCURENCES'], patientsMutatedGenes['GENE_INDEX'])]


    # Some patients that have CNAs in all samples do not have any mutated gene
    pWithoutSNVs = list(set(patientCNAsTCGA[numberFiles].keys()) - set(patientsMutatedGenes['PATIENT_ID']))

    if(len(pWithoutSNVs) > 0):        
        [patientsWithoutSNVs(patientsMutatedGenes, patientID) for patientID in pWithoutSNVs]


    # Sort the patient identifiers alphabetically
    patientsMutatedGenes = patientsMutatedGenes.sort_values(['PATIENT_ID'])


    return dict(zip(patientsMutatedGenes['PATIENT_ID'].values, patientsMutatedGenes['SNVS_OCCURENCES']))


    

# Load the each gene somatic mutations for each of the TCGA tumors, including the gene name, the patient sample identifier, 
# the SNV chromosome, position, reference and alternative allele, the protein effect, among others
# Identify whether a gene that has a CNA has or not any SNVs intersecting it

"""filename, numberFiles, genesNamesIndexes""" # The name of the TCGA file (tumor) accessed in the current iteration
                                               # The current file number to access the patient genes list for a specific tumor
                                               # The genes names and their corresponding indexes

# returns: the TCGA genes SNV for each of the tumors 
# the information of whether a certain gene has or not any SNVs
def loadSNVsTCGA(filename, numberFiles, genesNamesIndexes):
    sampleGenesSNVsPatients = pd.read_csv(filename, delimiter=';', \
        dtype = {'GENE_NAME' : str, 'CHROM': 'category', 'POS': 'uint32', 'REF_ALLELE': 'category', 'ALT_ALLELE': 'category', 'PATIENT_SNV': 'category', \
                 'REF_COUNTS': 'uint16', 'ALT_COUNTS' : 'uint16', 'COUNTS_DEPTH' : 'uint16', 'PROTEIN_EFFECT': 'category', 'FILTER': 'category'}, \
            names = ['GENE_NAME', 'CHROM', 'POS', 'REF_ALLELE', 'ALT_ALLELE', 'PATIENT_SNV', 'REF_COUNTS', \
                  'ALT_COUNTS', 'COUNTS_DEPTH', 'PROTEIN_EFFECT', 'FILTER'], header = 0) # using all columns

    
    # consider only the TCGA SNVs which filter condition is PASS
    sampleGenesSNVsPatients = sampleGenesSNVsPatients[sampleGenesSNVsPatients['FILTER'] == 'PASS']

    # remove the filter column
    sampleGenesSNVsPatients.drop(columns = ['FILTER'], axis=1, inplace=True)


    # we will not look at X and Y chromosomes 
    sampleGenesSNVsPatients = sampleGenesSNVsPatients[(sampleGenesSNVsPatients['CHROM'] != "chrX") & (sampleGenesSNVsPatients['CHROM'] != "chrY")]

    # drop any rows with nan values (just for verification, there seems to be no empty/nan cells)
    sampleGenesSNVsPatients.dropna(axis=0, how='any', inplace=True)


    # the type of each tumor based on the file name
    sampleGenesSNVsPatients['TUMOR_TYPE'] = [filename[filename.find('.') + 1: filename.find(".", filename.find(".") + 1)]] * len(sampleGenesSNVsPatients)

    # change data type to reduce memory usage
    sampleGenesSNVsPatients['TUMOR_TYPE'] = sampleGenesSNVsPatients['TUMOR_TYPE'].astype('category')


    # convert the row values from, e.g., chr1 to 1 and convert from string to uint8
    sampleGenesSNVsPatients['CHROM'] = sampleGenesSNVsPatients['CHROM'].str[3:].astype('uint8')

    # For the the mutation identifier column (to have a unique identifier/key of the SNV)
    sampleGenesSNVsPatients['MUTATION_ID'] = \
        'CHR' + sampleGenesSNVsPatients['CHROM'].astype(str) + ':' + sampleGenesSNVsPatients['POS'].astype(str) + ':' + \
            sampleGenesSNVsPatients['REF_ALLELE'].astype(str) + ':' + sampleGenesSNVsPatients['ALT_ALLELE'].astype(str)
        
        
    # To see the SNVs samples corresponding to each of the patients
    sampleGenesSNVsPatients['PATIENT_ID'] = sampleGenesSNVsPatients['PATIENT_SNV'].str[:12]     
        
    # change data type to reduce memory usage
    sampleGenesSNVsPatients['PATIENT_ID'] = sampleGenesSNVsPatients['PATIENT_ID'].astype('category')
    
    # rearrange columns
    sampleGenesSNVsPatients = sampleGenesSNVsPatients.iloc[:, [12, 5, 0, 1, 2, 3, 4, 11, 6, 7, 8, 9, 10]]



    # The SNVs sample identifiers do not need to match the CNAs sample identifiers

    # remove patients that have SNVs in a specific sample but do not have CNAs in that specific sample
    pWithoutCNAs = list(set(sampleGenesSNVsPatients['PATIENT_ID']) - set(cnasPositionsMaskTCGA[numberFiles]['PATIENT_ID']))

    sampleGenesSNVsPatients = sampleGenesSNVsPatients[~(sampleGenesSNVsPatients.PATIENT_ID.isin(pWithoutCNAs))]


    # removes genes that have SNVs but do no have CNAs in a specific sample
    # since cnasScoresTCGA has the same genes in any tumor, we can use any of the tumor samples
    gWithoutCNAs = list(set(sampleGenesSNVsPatients['GENE_NAME']) - set(gCNAsTumor))

    sampleGenesSNVsPatients = sampleGenesSNVsPatients[~(sampleGenesSNVsPatients.GENE_NAME.isin(gWithoutCNAs))]
    
        
    # Sort the genes alphabetically
    sampleGenesSNVsPatients = sampleGenesSNVsPatients.sort_values(['GENE_NAME'])
    
    # reset indexes to start on 0 and increment by 1
    sampleGenesSNVsPatients.index = range(len(sampleGenesSNVsPatients.index))
    

    # Identify which genes are intersected by SNVs and change the gene SNVs occurence value from 0 to 1 
    patientSNVsTCGA = genesSNVsTCGA(numberFiles, sampleGenesSNVsPatients, genesNamesIndexes)
    

    return sampleGenesSNVsPatients, patientSNVsTCGA




# Updates the patientsGenesInfoTCGA list with the CNAs and SNVs informations 

"""patientGenesInfoTCGA, patientCNAsValues, patientSNVsValues""" # each patient genes copy number and somatic mutations information list
                                                                 # each patient genes copy number representative values
                                                                 # each patient genes somatic mutations occurences informations

# returns: no return keyword, it updates the patientsGenesInfoTCGA list with the CNAs and SNVs informations 
def assignPatientsCNAsSNVsTCGA(patientGenesInfoTCGA, patientCNAsValues, patientSNVsValues):
    # All the even positions contain the sample copy number of a gene
    patientGenesInfoTCGA[0::2] = patientCNAsValues

    # All the odd positions contain the somatic mutations information of a gene
    patientGenesInfoTCGA[1::2] = patientSNVsValues


   

# Build the patients genes information table of the TCGA tumors with the CNAs and SNVs informations together

"""tumorSample""" # the number of the tumor sample

# returns: the patients genes information table of the TCGA tumors with the CNAs and SNVs informations
def buildPatientsCNAsSNVsTCGA(tumorSample):
    # create a list for the copy number values and somatic mutation occurences
    # create using a list comprehension to avoid lists being created as equal copies
    # 'i1' dtype because of the copy number representative values having -1 representing that the gene copy number < 2
    patientsGenesInfoTCGA = [np.array([0.5] * len(gCNAsTumor) * 2, dtype = 'i1') for patient in range(0, len(patientCNAsTCGA[tumorSample]))]

    # insert the patients CNAs and SNVs information for each tumor sample in the patientsGenesInfoTCGA list
    list(map(assignPatientsCNAsSNVsTCGA, patientsGenesInfoTCGA, list(patientCNAsTCGA[tumorSample].values()), list(patientSNVsTCGA[tumorSample].values())))


    return patientsGenesInfoTCGA    
    
    


# Build the matrix that will be used for the manifold with the mice samples and patients TCGA genes 

"""no parameters""" # uses global variables

# returns: the matrix which contains all mice samples and patients TCGA as rows and with the genes that both have in common as columns 
def commonGenesMatrix():
    # Genes that only appear in the miceSamplesGenesInfo matrix
    gMice = list(set(gCNAsAll) - set(gCNAsTumor))
    
    # Genes that only appear in the patientsGenesInfoTCGA matrix
    gPatients = list(set(gCNAsTumor) - set(gCNAsAll))
    
    
    cNamesMiceBoth = pd.DataFrame(cNamesMice, columns = ['GENE_MATRIX'])
    
    # Gets the genes names by removing the _cn and _sm of the copy number and somatic mutation occurrences identifiers
    cNamesMiceBoth = list(cNamesMiceBoth[ \
                            ~(cNamesMiceBoth['GENE_MATRIX'].str.extract(r'(.*?)_', expand = False).isin(gMice))]['GENE_MATRIX'])    
         
        
    # Eliminate genes that only appear in the miceSamplesGenesInfo matrix
    miceSamplesInfoConcat = miceSamplesGenesInfo[cNamesMiceBoth]

    # Insert tumor type column in mice samples genes information    
    miceSamplesInfoConcat.loc[:, 'TUMOR_TYPE'] = len(miceSamplesInfoConcat) * ['MICE_BRCA']
    
    
    cNamesPatientsBothTCGA = pd.DataFrame(cNamesPatientsTCGA, columns = ['GENE_MATRIX'])
    
    # Gets the genes names by removing the _cn and _sm of the copy number and somatic mutation occurrences identifiers
    cNamesPatientsBothTCGA = list(cNamesPatientsBothTCGA[ \
                                    ~(cNamesPatientsBothTCGA['GENE_MATRIX'].str.extract(r'(.*?)_', expand = False).isin(gPatients))]['GENE_MATRIX'])
       
    # include tumor type column in the matrix with the patients TCGA genes information
    cNamesPatientsBothTCGA = cNamesPatientsBothTCGA + ['TUMOR_TYPE']     
        
    
    # Eliminate genes that only appear in the patientsGenesInfoTCGA matrix
    patientsGenesInfoTCGAConcat = patientsGenesInfoTCGA[cNamesPatientsBothTCGA]
    
    
    # Concatenate both dataframes with the same genes as columns and rows with the mice samples and patients TCGA
    miceSamplesPatientsTCGA = pd.concat([miceSamplesInfoConcat, patientsGenesInfoTCGAConcat], axis = 0).sort_index()


    return miceSamplesPatientsTCGA




#--------------------------------------------------------------------------------------------------




#####  Load file with the chromosome, start and end positions, symbol and name of the genes to use in both the mice data and the TCGA data  #####


# Current directory: Dados 2 º Parte/


# To load the genes information
os.chdir("./TCGA")


# Get the genes that are in chromosomes X and Y, and the genes information table
genesChrXY, genesInfo = loadGenesInfo(list(glob.iglob("*.tab", recursive=False))[0])




#--------------------------------------------------------------------------------------------------




#####  Load the CNAs germline files and the SNVs annotations files (mice files)  #####


# To load the CNAs germline files
os.chdir("../CNVS/excel_format_germline")


# Mouse 49, Mouse 55, Mouse 61, Mouse 62
mice = ['M49', 'M55', 'M61', 'M62']

# Control Group, Sunitinib Group
groups = ['Ctrl', 'Sunit']

# Region 1 , Region 2, Region 3, Region 4
regions = ['R1', 'R2', 'R3', 'R4']


# Load the CNAs germline files for each mice samples and obtain which genes are intersected by which CNAs inside a chromosome of a same sample
cnasGenesInfo = list(map(loadMiceCNAs, sorted(glob.iglob("*.csv", recursive=False)), range(0, len(list(glob.iglob("*.csv", recursive=False))))))


# Unzip the returned pair of lists from the loadMiceCNAs method, defining the variables cnasGermline, cnasAnnotations 
# And miceSamplesID, the list of identifiers of the mice samples
cnasGermline, cnasAnnotations, miceSamplesID = map(list, zip(*cnasGenesInfo))


# remove variables to save space
del cnasGenesInfo




# distinct genes (contains all the genes, even the ones which do not have CNAs in every mice sample)
distinctGenes = set(genesInfo['GENE_NAME'])


# Remove the genes that do not have CNAs in every sample 
# Because, even if they have SNVs in every sample, we are considering that they need to have both CNAs and SNVs, with CNAs in every sample
gWithoutCNAsMice = set(itertools.chain(*list(map(lambda sampleCNAs: distinctGenes - set(sampleCNAs['GENE_NAME']), cnasAnnotations))))
                       

# Remove from each sample the CNAs that intersect genes that are not intersected by CNAs in every sample
# cnasGermline contains all the initial CNAs, including those of genes not in every sample

# genesInfo contains all the genesInfo initial information, except the duplicated entries and X and Y chromosomes
# Each mice and TCGA structures contain the genes that have CNAs in every sample

# Get the copy number for each mice sample of the genes alphabetically ordered
cnasMiceSamples = list(map(getMiceGenesCopyNumber, cnasAnnotations))


# Unzip the returned pair of lists from the getMiceGenesCopyNumber method, updating the variable cnasAnnotations and defining the variable miceSamplesGenesInfo
cnasAnnotations, miceSamplesGenesInfo = map(list, zip(*cnasMiceSamples))


# remove variables to save space
del cnasMiceSamples







# To load the SNVs annotations files
os.chdir("../../variantes_somaticas/anotacoes/excel_format")


# Since every sample will have only the genes that have CNAs in every sample, we will have the same genes in all samples (so we can use the genes from any mice sample)
gCNAsAll = cnasAnnotations[0]['GENE_NAME']


# Load which SNVs intersect which genes, and which genes intersected by CNAs are intersected by which SNVs
# With one fix argument initializing the genes somatic mutations occurences for each sample to 0
snvsGenesInfo = list(map(lambda filename, numberFiles: loadMiceSNVs(filename, numberFiles, np.array([0] * len(gCNAsAll)).astype('uint8')), \
                         sorted(glob.iglob("*.csv", recursive=False)), range(0, len(list(glob.iglob("*.csv", recursive=False))))))
   
    
# Unzip the returned pair of lists from the loadMiceSNVs method, defining the variables snvsAnnotations and genesAnnotations
snvsAnnotations, genesAnnotations = map(list, zip(*snvsGenesInfo))


# remove variables to save space
del snvsGenesInfo







# Build the input dataframe for the copy number and somatic mutation mice sample genes information


# for the copy number and somatic mutation info for each gene
cNamesMice = [''] * len(gCNAsAll) * 2

# All the even positions of the array contain the information about the copy number of each gene
# All the odd positions of the array contain the information about the occurence of a somatic mutation on each gene
cNamesMice[0::2] = (gCNAsAll + "_cn").to_list()
cNamesMice[1::2] = (gCNAsAll + "_sm").to_list()


# Dataframe with the copy number and somatic mutation information for all the genes alphabetically ordered in each of the samples
# Since a unique data type can be assigned, we will use int8 also for the somatic mutations information
miceSamplesGenesInfo = pd.DataFrame(np.array(miceSamplesGenesInfo), index = miceSamplesID, columns = cNamesMice, dtype = 'int8')




#--------------------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------------------







#####  Load the TCGA CNAs files and TCGA SNVs files (TCGA files)  #####


# To load the TCGA CNAs focal score files
os.chdir("../../../TCGA/CNV/excel_format")


# Load the focal scores of the CNAs that intersect the genes for each of the tumors 
# After processing it has the same genes in every CNAs focal scores of every tumor
scoresCNAsTCGA = list(map(focalScoresCNAsTCGA, sorted(glob.iglob("*.csv", recursive=False))))







# To load the TCGA CNAs masked copy number segments files
os.chdir("../")


# Load the TCGA CNAs segments information, including the CNAs sample, chromosome, start and end position, log2 and patient samples identifier
cnasInfoTCGA = list(map(lambda filename, numberFiles: copySegmentsCNAsTCGA(filename, numberFiles), \
                  sorted(glob.iglob("*.tab", recursive=False)), range(0, len(list(glob.iglob("*.tab", recursive=False))))))

    
# Unzip the returned pair of lists from the copySegmentsCNAsTCGA method, defining the variables cnasPositionsMaskTCGA and patientCNAsTCGA
cnasPositionsMaskTCGA, patientCNAsTCGA = map(list, zip(*cnasInfoTCGA))


# remove variables to save space
del cnasInfoTCGA







os.chdir("../SNV/excel_format")


# Since every sample has the same genes, we can use the scoresCNAsTCGA genes of the first tumor
gCNAsTumor = scoresCNAsTCGA[0]['GENE_NAME'].tolist()

 
# Load the TCGA somatic mutations information, containing the gene name, chromosome, SNV position, reference and alternative allele, protein effect, among others
# Find the patients with genes mutated in each of the TCGA tumors
# With one fix argument containing a dataframe with the genes names and their corresponding indexes to be used to identify each of the mutated genes names indexes
snvsInfoTCGA = list(map(lambda filename, numberFiles: loadSNVsTCGA(filename, numberFiles, 
                    pd.DataFrame(zip(gCNAsTumor, list(scoresCNAsTCGA[0].index)), columns = ['GENE_NAME', 'GENE_INDEX'])),
                        sorted(glob.iglob("*.csv", recursive=False)), range(0, len(list(glob.iglob("*.csv", recursive=False))))))


# Unzip the returned pair of lists from the loadSNVsTCGA method, defining the variables snvsGenesInfoTCGA and patientSNVsTCGA
snvsGenesInfoTCGA, patientSNVsTCGA = map(list, zip(*snvsInfoTCGA))


# remove variables to save space
del snvsInfoTCGA







# Build the final matrix with the TCGA patients and the copy numbers and somatic mutations genes information


# Get the copy number and somatic mutations information of the patients genes alphabetically ordered for each tumor sample
patientsGenesInfoTCGA = list(map(buildPatientsCNAsSNVsTCGA, range(0, len(patientCNAsTCGA))))


# for the copy number and somatic mutation info for each gene
cNamesPatientsTCGA = [''] * len(gCNAsTumor) * 2

# To build the columns names efficiently
gCNAsTumor = pd.Series(gCNAsTumor)

# All the even positions of the array contain the information about the copy number of each gene
# All the odd positions of the array contain the information about the occurence of a somatic mutation on each gene
cNamesPatientsTCGA[0::2] = (gCNAsTumor + "_cn").to_list()
cNamesPatientsTCGA[1::2] = (gCNAsTumor + "_sm").to_list()


# Dataframe with the copy number and somatic mutation information for all the genes alphabetically ordered in each of the samples
# Since a unique data type can be assigned, we will use int8 also for the somatic mutations information
# Row indexes correspond to the patients names
patientsGenesInfoTCGA = pd.concat(list(map(lambda tumorSample, patientsIDs: pd.DataFrame(np.array(tumorSample), index = list(patientsIDs), \
                            columns = cNamesPatientsTCGA, dtype = 'int8'), patientsGenesInfoTCGA, patientCNAsTCGA)), axis = 0)
   
    
# Get the tumor type of each patient to be used as target variable when applying the machine learning algorithms   
patientsGenesInfoTCGA['TUMOR_TYPE'] = list(itertools.chain(*list(map(lambda tumorType, patientsWithTumor: \
                                           [tumorType.loc[0]['TUMOR_TYPE']] * len(patientsWithTumor), cnasPositionsMaskTCGA, patientCNAsTCGA))))
  
# sort the patients alphabetically
patientsGenesInfoTCGA = patientsGenesInfoTCGA.sort_index()    



    
#--------------------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------------------







#####  Save the miceSamplesGenesInfo and patientsGenesInfoTCGA dataframe objects  #####


# To use each dataframe separately to apply the machine learning algorithms   
os.chdir("../../../InputTablesFilesML")


# Mice samples genes matrix with information about their copy number and somatic mutations occurences
miceSamplesGenesInfo.to_pickle("miceSamplesGenesInfo.pkl")

# Patients TCGA genes matrix with information about their copy number and somatic mutations occurences
patientsGenesInfoTCGA.to_pickle("patientsGenesInfoTCGA.pkl") 




# Merge the mice samples and patients TCGA genes matrix with information about their copy number and somatic mutations occurences in one


# To build a manifold of the mice samples and patients TCGA genes they have to have the same set of genes for analysis
miceSamplesPatientsTCGA = commonGenesMatrix()

# Mice samples and patientsTCGA genes matrix with information about their copy number and somatic mutations occurences
miceSamplesPatientsTCGA.to_pickle("miceSamplesPatientsTCGA.pkl")






        