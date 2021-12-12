#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 18:07:23 2021

@author: david
"""



#####  Libraries  #####


# To open system files
import os

# Pandas also for data management
import pandas as pd

# For transformations between lists and dataframes
import itertools

# import time
# start_time = time.time()




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





#--------------------------------------------------------------------------------------------------




#####  Load SNVs and CNAs input data  #####


# when it equals one, we do not use the CNAs with major copy number = 2 and minor copy number = 0
notAllCNAs = 0


# print(os.getcwd())
# can also change the python directory above in spyder accordingly (Dados/CNVS)


os.chdir("./excel_format_germline")

cnasInfo = []

# sort filenames
for filename in sorted(os.listdir(os.getcwd())):
     with open(os.path.join(os.getcwd(), filename), 'r'):
         # delimiter = ';' on windows csv and '\t' on linux csv
         data = pd.read_csv(filename, delimiter='\t')
         
         # remove chrY since is breast cancer (women only have X chromosome)
         data = data[:-1]
        
         if(notAllCNAs == 0):
             # for the chromosome X CNAs method errors
             data = data.drop( data[((data['cn1'] == 0) & (data['cn2'] == 0))].index )
         else:
             data = data.drop( data[((data['cn1'] == 0) & (data['cn2'] == 0)) | ((data['cn1'] == 2) & (data['cn2'] == 0)) ].index )

  
         # drop rows with nan values
         data = data.dropna() 
    
         # reset indexes to start on 0 and increment by 1
         data = data.reset_index(drop = True)
    
    
         data = data.rename(columns = {data.columns[0]: 'CHROM', data.columns[1]: 'START POS', data.columns[2]: 'END POS',
                                data.columns[8]: 'CN1', data.columns[9]: 'CN2'})
             
         
         data = data[['CHROM', 'START POS', 'END POS', 'CN1', 'CN2']]
        
         cnasInfo.append(data)






os.chdir("../../variantes_somaticas/excel_format")


genomeInfo = []

for filename in sorted(os.listdir(os.getcwd())):
    with open(os.path.join(os.getcwd(), filename), 'r'):
        data = pd.read_csv(filename, delimiter=';')
        
        data = data[data['FILTER'] == 'PASS']
        
        # will not be included because of copy number regions identification method's error
        data = data[data['CHROM'] != "chrX"]
        
        # reset indexes to start on 0 and increment by 1
        data = data.reset_index()
        
        data = data[['CHROM', 'POS', 'TUMOR', 'ALT']]
        
   
        # Add major_cn, minor_cn, normal_cn and tumour_content columns
        data.loc[:, 'major_cn'] = -1
        data.loc[:, 'minor_cn'] = -1
        data.loc[:, 'normal_cn'] = 2
        data.loc[:, 'tumour_content'] = 1
        
        
        genomeInfo.append(data)




#--------------------------------------------------------------------------------------------------




#####  Intersect SNVs and CNAs genome regions by sample  #####


i = 0
    
while (i < len(genomeInfo)):
    j = 0
    
    while j < len(genomeInfo[i]):
        h = 0
        
        compareWith = cnasInfo[i][cnasInfo[i]['CHROM'] == genomeInfo[i].iloc[j]['CHROM']]


        while h < len(compareWith):        
           

            if ((compareWith.iloc[h]['START POS'] <= genomeInfo[i].iloc[j]['POS']) & 
                    (genomeInfo[i].iloc[j]['POS'] <= compareWith.iloc[h]['END POS'])):
                        
                        
                        genomeInfo[i].at[genomeInfo[i].iloc[[j]].index.values[0],'major_cn'] =  compareWith.iloc[h]['CN1']
                        genomeInfo[i].at[genomeInfo[i].iloc[[j]].index.values[0],'minor_cn'] =  compareWith.iloc[h]['CN2']
                            
                        
                        break
     
            h+=1
        
        j+=1
    
        
    # if the row value is -1, that means that it was not possible to estimate cn1 and cn2  
    genomeInfo[i] = genomeInfo[i].drop( genomeInfo[i][ (genomeInfo[i]['minor_cn'] == -1) |
                      (genomeInfo[i]['major_cn'] == -1)  ].index )
    
        
    i+=1




#--------------------------------------------------------------------------------------------------




#####  Create the input files for pyClone  #####


os.chdir("../../CNVS/pyCloneInputFiles")

# tsv file order:
    
# 1 - mutation_id
# 2 - sample_id
# 3 - ref_counts
# 4 - alt_counts 
# 5 - major_cn
# 6 - minor_cn
# 7 - normal_cn
# 8 - tumour_content
    

# Mouse 49, Mouse 55, Mouse 61, Mouse 62
mice = ['M49', 'M55', 'M61', 'M62']

groups = ['Ctrl', 'Sunit']

# Region 1 , Region 2, Region 3, Region 4 
regions = ['R1', 'R2', 'R3', 'R4']


# For the file with all the mice info together
to_output = []


# For each one of the mice files
mice_files_output = [[] for i in itertools.repeat(None, 4)]
     


for i in range(0, len(genomeInfo)):

    
    # For ref_counts and alt_counts columns
    tumorInfo = genomeInfo[i]['TUMOR'].to_list()
    reads = [readCounts(mut) for mut in tumorInfo]
    
    
    for j in range(0, len (reads)):
    
        genomeInfo[i].at[genomeInfo[i].iloc[[j]].index.values[0],'ref_counts'] = reads[j][0]     
        genomeInfo[i].at[genomeInfo[i].iloc[[j]].index.values[0],'alt_counts'] = reads[j][1] 
    
    
        # For mutation_id column
        genomeInfo[i].at[genomeInfo[i].iloc[[j]].index.values[0],'mutation_id'] = \
        ''.join(( str(genomeInfo[i].iloc[j]['CHROM']).upper(), ':', str(genomeInfo[i].iloc[j]['POS']), ':', str(genomeInfo[i].iloc[j]['ALT'])) )
    
    
        # For sample_id column    
        genomeInfo[i].at[genomeInfo[i].iloc[[j]].index.values[0],'sample_id'] = mice[int(i/4)] + groups[int(i/8)] + regions[i % 4]  
            

    genomeInfo[i] = genomeInfo[i][['mutation_id', 'sample_id', 'ref_counts', 'alt_counts', 'major_cn', 'minor_cn', 'normal_cn', 'tumour_content']]
    
    
    # to convert the read counts from float to int on the tsv output
    # https://stackoverflow.com/questions/17092671/python-pandas-output-dataframe-to-csv-with-integers
    
    genomeInfo[i]['ref_counts'] = genomeInfo[i]['ref_counts'].astype(int)
    genomeInfo[i]['alt_counts'] = genomeInfo[i]['alt_counts'].astype(int)



    # For the file with all mice together
    to_output.append(genomeInfo[i].values)


    # For the 4 files    
    mice_files_output[int(i/4)].append(genomeInfo[i].values)




#--------------------------------------------------------------------------------------------------




# For all mice together include mutations in all samples they are not in


to_output = pd.DataFrame (list(itertools.chain.from_iterable((to_output))), \
        columns = ['mutation_id', 'sample_id', 'ref_counts', 'alt_counts', 'major_cn', 'minor_cn', 'normal_cn', 'tumour_content'])
       
    
    
    
# include mutations info in all the samples they are not in   
samples_count = pd.DataFrame(to_output['mutation_id'].value_counts())

# the ones that are not in all samples
to_introduce = samples_count.loc[samples_count['mutation_id'] != 16]


samples_id = pd.DataFrame(to_output.sample_id.unique())



d = {}

h = 0


# each mutation_id not present in all samples
for i in to_introduce.index:
    
    # for each mutation see in which samples it is
    in_samples = pd.DataFrame(to_output[to_output['mutation_id'] == i].sample_id.unique())

    # gives the samples the mutation is not present in
    include_samples = pd.concat([samples_id, in_samples]).drop_duplicates(keep = False)
    
    
    for j in range(0, len(include_samples)):
    
        d[h] = pd.Series(data = {'mutation_id': i, 'sample_id': include_samples.values[j][0], 'ref_counts': 0, 'alt_counts': 0, 'major_cn': 1, 'minor_cn': 1, 'normal_cn': 2, 'tumour_content': 1})
        
        h += 1


df = pd.DataFrame.from_dict(d, "index")

to_output = pd.concat([to_output, df], ignore_index = True)

  
  
  
if(notAllCNAs == 0):
    with open('allMice_genomeInfo.tsv', 'w') as f:    
        f.write(to_output.to_csv(sep='\t', header = True, index=False))
    
else: # with = 1, change the name of the tsv file to open
    with open('allMice_genomeInfoNotAllCNAs.tsv', 'w') as f:    
        f.write(to_output.to_csv(sep='\t', header = True, index=False))





#--------------------------------------------------------------------------------------------------




# For each mouse include mutations in all mouse samples they are not in



miceInfo = [pd.DataFrame(list(itertools.chain.from_iterable(mice_file)), \
        columns = ['mutation_id', 'sample_id', 'ref_counts', 'alt_counts', 'major_cn', 'minor_cn', 'normal_cn', 'tumour_content']) \
                for mice_file in mice_files_output]    
    

    
    
# For each of the 4 mice files include info about mutations in all the samples they are not in   

to_introduce = []

mice_samples = {}

h = 0

for mouse in miceInfo:

    samples_count = pd.DataFrame(mouse['mutation_id'].value_counts())

    # 4 regions samples by each mice file
    to_introduce.append(samples_count.loc[samples_count['mutation_id'] != 4])
    
    mice_samples[h] = pd.DataFrame(mouse.sample_id.unique())

    h += 1




mice_dict = {}

t = 0


# indexes 0 to 3 corresponding to mouse 49, 55, 61 and 62 respectively  
for mouse in to_introduce:
    
    d = {} # dictionary for each mouse mutations' samples to be introduced in
    
    h = 0
    
    # each mutation_id that is not in all 4 samples of a mouse 
    for mut in mouse.index:
    
        in_samples = pd.DataFrame(miceInfo[t][ miceInfo[t]['mutation_id'] == mut ].sample_id.unique())
        
        include_samples = pd.concat([mice_samples[t], in_samples]).drop_duplicates(keep = False)

        d[''] = {} # dictionary of each mutations' samples to be introduced in for a specific mouse
        

        for s in range(0, len(include_samples)):
            
              d[h] = pd.Series(data = {'mutation_id': mut, 'sample_id': include_samples.values[s][0], 'ref_counts': 0, 'alt_counts': 0, 'major_cn': 1, 'minor_cn': 1, 'normal_cn': 2, 'tumour_content': 1})

              h += 1
        

    df = pd.DataFrame.from_dict(d, "index")

    mouse_samples = pd.concat([miceInfo[t], df], ignore_index = True)


    mice_dict[t] = mouse_samples

    t += 1    
    
    

    


# For multiple files:
        
for i in range(0, len(mice_dict)):
    
    if(notAllCNAs == 0):
        with open(mice_dict[i]['sample_id'].values[0][0:3] + "_genomeInfo.tsv", 'w') as f:
            f.write(mice_dict[i].to_csv(sep='\t', header = True, index = False)) 

    else: # with = 1, change the name of the tsv file to open
        with open(mice_dict[i]['sample_id'].values[0][0:3] + "_genomeInfoNotAllCNAs.tsv", 'w') as f:
            f.write(mice_dict[i].to_csv(sep='\t', header = True, index = False)) 
    
    
    
    
    
# print("--- %s seconds ---" % (time.time() - start_time))




