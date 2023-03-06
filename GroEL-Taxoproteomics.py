#!/usr/bin/env python
# coding: utf-8

# In[1]:


#This script was developed by Simon Klaes and Lorenz Adrian for GroEL-Taxoproteomics. 
# If using this script, please cite:
#-------------------------------------------------------------------------
# Klaes, S., Madan, S., Deobald, D., Cooper, M., and Adrian, L., 2023. 
# Population Analysis of Mixed Bacterial Communities by GroEL-Taxoproteomics, 
# Frontiers in Microbiology
#-------------------------------------------------------------------------
# Thank you; it will help us fund more script development.

If you use this script, please cite: 

import itertools
from itertools import product
import pandas as pd
from pandas.core.common import SettingWithCopyWarning
import warnings
import tkinter as tk
from tkinter import filedialog
import os

def make_patterns(s):                                     #This  function generates all possible peptide sequences due to the indistinguishability of leucine and isoleucine. The input is a string of an amino acid sequence, and the output is a generator that yields all possible sequences resulting from replacing the "IL" amino acids in the input string with either "I" or "L".
    IndistinguishableAminoAcids = 'IL'                                                            
    seq = list(s)                         
    indices = [i for i, c in enumerate(seq) if c in IndistinguishableAminoAcids]                                                           
    for t in product(IndistinguishableAminoAcids, repeat=len(indices)):
        for i, c in zip(indices, t):
            seq[i] = c
        yield ''.join(seq)
        
def most_frequent(List):                                  # This function takes a list as an input and returns the most frequent value in the list.
    return max(set(List), key = List.count)              

def most_frequent_frequency_proportion(list_data):       # This function takes a list as an input and returns the proportion of the most frequent value in the list.
    stat_frequency = {}                                  
    stat_proportion = {}                                 
    total = len(list_data)                               
    for e in list_data:                                   
        if str(e) in stat_frequency:                     
            stat_frequency[str(e)] += 1
        else:                                            
            stat_frequency[str(e)] = 1
    for key, value in stat_frequency.items():            
        stat_proportion[key] = value / total             
    return max(stat_proportion.values())                 

def remove_duplicates_from_list(x):                       #This function removes duplicates from a list
    return list(dict.fromkeys(x))

def get_duplicates_from_list(L):                          #This function returns duplicate values from a list
    seen = set()
    seen2 = set()
    seen_add = seen.add
    seen2_add = seen2.add
    for item in L:
        if item in seen:
            seen2_add(item)
        else:
            seen_add(item)
    return list(seen2)

def get_intersection_from_two_lists(x,y):               #This function returns the intersecition of two lists
    return list(set(x).intersection(y))

warnings.filterwarnings("ignore", category=SettingWithCopyWarning)

#ask for all paths needed for the script
root = tk.Tk()
root.withdraw()
collapsed_MetaProSIP_Output_filename = filedialog.askopenfilename(filetypes=[("collapsed MetaproSIP Output", ".tsv .csv")], title='Select the collapsed MetaProSIP Output')
target_database_filename = filedialog.askopenfilename(title='Select the target tsv-database')
cRAP_database_filename = filedialog.askopenfilename(title='Select the contaminants tsv-database')
directory_for_the_output_of_this_script = filedialog.askdirectory(title='Select the folder to save the results in')

#load collapsed_MetaProSIP_Output
print("Loading collapsed MetaProSIP Output from " + collapsed_MetaProSIP_Output_filename)
collapsed_MetaProSIP_Output = pd.read_csv(collapsed_MetaProSIP_Output_filename, sep = '\t')
print("Collapsed MetaProSIP Output loaded successfully")

# load the trypsin-digest tsv of target database and contaminant database in the respective dataframe
print("Loading target database from " + target_database_filename)
df_target_database = pd.read_csv(target_database_filename, sep = '\t')
print("Target database loaded successfully")
print("Loading cRAP database from " + cRAP_database_filename)
df_contaminants = pd.read_csv(cRAP_database_filename, sep = '\t')
print("cRAP database loaded successfully")

#create a set (non redundant) of all tryptics peptides in the target database and contaminant database, respectively
all_tryp_peps_set = set(itertools.chain.from_iterable(df_target_database.TrypPeps.apply(lambda x: x[2:-2].split("', '"))))
all_contaminant_peps_set = set(itertools.chain.from_iterable(df_contaminants.TrypPeps.apply(lambda x: x[2:-2].split("', '"))))

#convert the tryptic peptides column to lists to facilitate counting of hits
df_target_database['TrypPeps'] = df_target_database.TrypPeps.apply(lambda x: x[2:-2].split("', '"))
df_contaminants['TrypPeps'] = df_contaminants.TrypPeps.apply(lambda x: x[2:-2].split("', '"))

#MetaProSIP Output files are loaded subsequently and peptide sequences are extracted to a list "found_seqs". Modifications (Carbamidomythylation and Oxidation) are being removed
sample_names = collapsed_MetaProSIP_Output['Sample'].unique().tolist()
for sample_name in sample_names:
    df_database = df_target_database
    found_seqs = []
    print("Loading sample ", sample_names.index(sample_name)+1, "/", len(sample_names), sample_name)
    MetaProSIP_Output = collapsed_MetaProSIP_Output[collapsed_MetaProSIP_Output['Sample'] == sample_name]
    python_output_filename = os.path.join(directory_for_the_output_of_this_script,sample_name.replace('.mzML', '')+".tsv")
    python_intermediate_filename = os.path.join(directory_for_the_output_of_this_script,sample_name.replace('.mzML', '')+"_intermediate.tsv")
    python_contaminants_output_filename = os.path.join(directory_for_the_output_of_this_script,sample_name.replace('.mzML', '')+"_contaminants.tsv")
    found_seqs = MetaProSIP_Output['Peptide Sequence'].tolist()
    found_seqs = [str(x).replace('(Carbamidomethyl)', '').replace('(Oxidation)', '') for x in found_seqs]
    
    #The found_seqs list is extended to found_seqs_wobbled by adding the isobaric equivalent of the detected peptides if it is also in the database. Peptides that are also in the contaminant database are excluded from found_seqs_wobbled    
    found_seqs_wobbled = []
    found_contaminant_peps_wobbled = []
    
    for i in found_seqs:
        for s in make_patterns(i):
            if s in all_contaminant_peps_set:
                if s not in found_contaminant_peps_wobbled:
                    found_contaminant_peps_wobbled.append(s)
            if s in all_tryp_peps_set:
                if s not in all_contaminant_peps_set:
                    if s not in found_seqs_wobbled:
                        found_seqs_wobbled.append(s)
            
    # set the counts to 0 before start counting and reset the following dataframes
    df_database.loc[:,"HitPepSeq"] = str("")
    df_database.loc[:,"Count"] = 0 
    df_HPSM = pd.DataFrame()
    df_results = pd.DataFrame()
    df_results_min5_family = pd.DataFrame()
    df_results_min5 = pd.DataFrame()
    total_hit_count = 0
    df_contaminants.loc[:,"HitPepSeq"] = str("")
    df_contaminants.loc[:,"Count"] = 0
    total_hit_count_contaminants = 0
    taxonomic_columns = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']

    # loop through contaminant peptides
    for contaminant_peptide in found_contaminant_peps_wobbled:
        # List with True/False - true if peptide is in contaminant database
        mask=[contaminant_peptide in x for x in df_contaminants['TrypPeps']] 
        # apply mask: for all proteins in the contaminant database that contain this tryptic peptide (all rows in which mask is true), "count" column is increased by 1 and the peptide is added to "HitPepSeq" column
        df_contaminants.loc[mask,"Count"] += 1
        df_contaminants.loc[mask,"HitPepSeq"] = df_contaminants.loc[mask,"HitPepSeq"] + contaminant_peptide + str(",")
        # print how often each peptide was found in the contaminant database
        print(contaminant_peptide, found_contaminant_peps_wobbled.index(contaminant_peptide)+1, "/", len(found_contaminant_peps_wobbled), "No. of hits in the contaminant database:", len(df_contaminants[mask]))
        total_hit_count_contaminants += len(df_contaminants[mask])
    # print how often contaminant peptides were found in total in the contaminant database
    print("Total found hits over all contaminant peptides:", total_hit_count_contaminants)
    df_contaminants['HitPepSeq'] = df_contaminants.HitPepSeq.apply(lambda x: x.split(","))
    #save contaminants dataframe as tsv-file
    df_contaminants.to_csv(python_contaminants_output_filename, sep = '\t')

    # loop through the detected GroEL peptides
    for peptide in found_seqs_wobbled:
        # List with True/False - true if peptide is in GroEL database
        mask = [peptide in x for x in df_database['TrypPeps']]
        # apply mask: for all proteins in the GroEL database that contain this tryptic peptide (all rows in which mask is true), "count" column is increased by 1 and the peptide is added to "HitPepSeq" column
        df_database.loc[mask,"Count"] += 1
        df_database.loc[mask,"HitPepSeq"] = df_database.loc[mask,"HitPepSeq"] + peptide + str(",")
        # print how often each peptide was found in the GroEL database
        print(peptide, found_seqs_wobbled.index(peptide)+1, "/", len(found_seqs_wobbled), "No. of hits in the database:", len(df_database[mask]))
        total_hit_count += len(df_database[mask])
    # print how often the peptides were found in total in the GroEL database    
    print("Total found hits over all peptides:", total_hit_count)

    #Create a new dataframe (df_HPSM) based on df_database that merges redundant sets of HitPepSeqs, counts its peptides and gives respective Acc NOs (and count of Acc NOs) of matching organisms + their taxonomic information
    df_database = df_database.drop(columns = ['ProtSeq', 'TrypPeps', 'Unnamed: 0']) #removes unnecessary columns
    df_database = df_database.astype(str) #changes all columns to str, because this is necessary for the following functions
    df_database = df_database.groupby(df_database.HitPepSeq)
    df_HPSM = df_database.agg("first")
    df_HPSM.update(df_database.agg({"AccNo": ",".join, "Kingdom": ",".join, "Phylum": ",".join, "Class": ",".join, "Order": ",".join, "Family": ",".join, "Genus": ",".join})) #merges all rows with the same HitPepSeq into one row
    df_HPSM['AccNo_Count'] = df_HPSM['AccNo'].apply(lambda x: x.count(',') + 1) #give the number of AccNos (number of GroEL proteins with that HitPepSeq)
    df_HPSM['Count'] = df_HPSM['Count'].astype(int) #transforms 'Count' from str to int to sort in the next step
    df_HPSM = df_HPSM.sort_values('Count', ascending = False) #sorts the dataframe by 'Count' in descending order
    df_HPSM.reset_index(level = 0, inplace=True) #sets new index
    df_HPSM['HitPepSeq'] = df_HPSM['HitPepSeq'].str[:-1] #removes last ',' in HitPepSeq
    df_HPSM.to_csv(python_intermediate_filename, sep = '\t') #saves df_HPMS as tsv-file

    #Create a new dataframe (df_results) based on df_HPSM
    df_results = df_HPSM.drop(df_HPSM.index[-1]) #removes last row (no peptide hits)
    df_results['HitPepSeq'] = df_results['HitPepSeq'].str.split(',') #transforms HitPepSeqs from str to lists
    df_results[taxonomic_columns] = df_results[taxonomic_columns].apply(lambda col: col.str.split(',')) #convert the taxonomic information about the organisms of the respective GroEL account numbers to lists to faciliate analysis
    for col in taxonomic_columns: #create new columns with the most frequent classification on the respective taxonomic level and its frequency in decimal numbers
        df_results['Top_' + col] = df_results[col].apply(lambda row: most_frequent(row))
        df_results['Top_' + col + '_fq'] = df_results[col].apply(lambda row: most_frequent_frequency_proportion(row))
    
    #the following is an intermediate step to calculate the TopRank_Count, shared_TopRank_Count, and strict_TopRank_Count:
    unique_count_list = df_results['Count'].unique().tolist()
    df_results['HitPepSeq_of_same_count'] = ""
    df_results['new_HitPepSeq'] = ""
    HitPepSeq_list = []
    #create a new column that contains all peptides that have been detected for each Count ('HitPepSeq_of_same_count') and a new column that contains all peptides that have been detected for each Count but not in for a higher Count ('new_HitPepSeq')
    for unique_count in unique_count_list:
        mask=(df_results['Count'] == unique_count)
        for each_count in df_results.loc[mask,'HitPepSeq']:
            for each_peptide in each_count:
                df_results.loc[mask,'HitPepSeq_of_same_count'] = df_results.loc[mask,'HitPepSeq_of_same_count']+each_peptide+","
                if each_peptide not in HitPepSeq_list:
                    HitPepSeq_list.append(each_peptide)
                    df_results.loc[mask,'new_HitPepSeq'] = df_results.loc[mask,'new_HitPepSeq']+each_peptide+","
    
    df_results['HitPepSeq_of_same_count'] = df_results['HitPepSeq_of_same_count'].str[:-1].str.split(',').apply(lambda row: remove_duplicates_from_list(row))
    df_results['new_HitPepSeq'] = df_results['new_HitPepSeq'].str[:-1].str.split(',').apply(lambda row: remove_duplicates_from_list(row))

    df_results['intersection_HitPepSeqs'] = [list(set(a).intersection(b)) for a, b in zip(df_results.HitPepSeq, df_results.new_HitPepSeq)]
    if total_hit_count > 0:
        df_results['len_intersection_HitPepSeqs'] = df_results['intersection_HitPepSeqs'].str.len()

    list_of_all_intersection_HitPepSeqs = list(itertools.chain.from_iterable(df_results['intersection_HitPepSeqs'].tolist()))
    shared_Top_HitPepSeqs = get_duplicates_from_list(list_of_all_intersection_HitPepSeqs)

    df_results['shared_Top_HitPepSeqs'] = df_results['intersection_HitPepSeqs'].apply (lambda row: get_intersection_from_two_lists(row, shared_Top_HitPepSeqs))
    df_results['single_Top_HitPepSeqs'] = (df_results['intersection_HitPepSeqs'].map(set) - df_results['shared_Top_HitPepSeqs'].map(set)).map(list)
    #adds new columns:  TopRank_Count, shared_TopRank_Count, and strict_TopRank_Count.
    #strict_TopRank_Count counts the number of peptides that are not present in another HitPepSeq with the same amount or more peptides
    #shared_TopRank_Count counts the peptides that are present in another HitPepSeq with the same amount of peptides. 
    #The sum of strict and shared TopRank_Count is called TopRank_Count. We use this value for filtering to achieve maximum parsimony.
    if total_hit_count > 0:
        df_results['shared_TopRank_Count'] = df_results['shared_Top_HitPepSeqs'].str.len()
        df_results['strict_TopRank_Count'] = df_results['len_intersection_HitPepSeqs'] - df_results['shared_TopRank_Count']
        df_results['TopRank_Count'] = df_results['len_intersection_HitPepSeqs']
    
    #if no peptides were found, all Counts are = 0
    if total_hit_count == 0:
        df_results['shared_TopRank_Count'] = 0
        df_results['strict_TopRank_Count'] = 0
        df_results['TopRank_Count'] = 0

    #rearrange columns to improve readability
    df_results = df_results[['HitPepSeq', 'Count', 'TopRank_Count', 'shared_Top_HitPepSeqs', 'shared_TopRank_Count', 'single_Top_HitPepSeqs', 'strict_TopRank_Count', 'AccNo', 'AccNo_Count', 'Kingdom', 'Top_Kingdom', 'Top_Kingdom_fq', 'Phylum', 'Top_Phylum', 'Top_Phylum_fq', 'Class', 'Top_Class', 'Top_Class_fq', 'Order', 'Top_Order', 'Top_Order_fq', 'Family', 'Top_Family', 'Top_Family_fq', 'Genus', 'Top_Genus', 'Top_Genus_fq']]

    #save df_results to tsv-file
    df_results.to_csv(python_output_filename, sep = '\t')

    #Filter df_results to create new dataframe based on taxonomy
    print("Creating family-based output file...")
    family_based_analysis_output_filename = os.path.join(directory_for_the_output_of_this_script,sample_name.replace('.mzML', '')+"family_based_min5.tsv")
    df_results_min5 = df_results[df_results['TopRank_Count'] > 4] #by default GroEL protein groups are filtered by a TopRank_Count >4. If you want to adjust this, change here.
    
    df_results_min5_family = df_results_min5.groupby('Top_Family').agg({'HitPepSeq': 'sum'}) #by default GroEl protein groups are merged at the family level. If you want to adjust this, change here to the level of interest.
    df_results_min5_family.reset_index(level = 0, inplace = True)

    df_results_min5_family['HitPepSeq'] = df_results_min5_family['HitPepSeq'].apply (lambda row: remove_duplicates_from_list(row))

    df_results_min5_family['non_redundant_peptide_count'] = df_results_min5_family['HitPepSeq'].apply(lambda row: len(row))

    df_results_min5_family['Int'] = 0
    
    MetaProSIP_Output.loc[:, 'Peptide Sequence'] = MetaProSIP_Output['Peptide Sequence'].astype(str)
    
    for ind in df_results_min5_family.index:
        peplist = df_results_min5_family['HitPepSeq'][ind]
        for peptide in peplist:
            for IL_peptide in make_patterns(peptide):
                if IL_peptide in MetaProSIP_Output['Peptide Sequence'].apply(lambda x: x.replace('(Carbamidomethyl)', '').replace("(Oxidation)", "")):
                    if IL_peptide not in peplist:
                        peplist.append(IL_peptide)
        mo = MetaProSIP_Output[MetaProSIP_Output['Peptide Sequence'].apply(lambda x: x.replace('(Carbamidomethyl)', '').replace("(Oxidation)", "")).isin(peplist)]
        df_results_min5_family.loc[ind, 'Int'] = df_results_min5_family.loc[ind, 'Int'] + sum(mo['INT 1'])
    
    df_results_min5_family
    df_results_min5_family.to_csv(family_based_analysis_output_filename, sep = '\t')

