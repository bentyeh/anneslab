
# coding: utf-8

# In[1]:


# Import all packages
import os
from multiprocessing.pool import ThreadPool
import numpy as np
import pandas as pd
from Bio import Entrez
from IPython.display import display


# In[2]:


# Global parameters / values
data_dir = "../data/"
dataAux_dir = "../data_aux/"
Klaeger_filename = "Klaeger.csv"
Huang_filename = "Huang.csv"
Annes100_filename = "Annes100.csv"
Annes500_filename = "Annes500.csv"

EntrezEmail = "example@gmail.com"

# Name of column of official gene symbols to be added to data frames
geneSymbolColumn = "GeneSymbol"

# Number of threads to use for looking up official gene symbols via Entrez
# - Possible values
#     None: uses a ThreadPool of a default number of threads as returned by os.cpu_count()
#     1: uses main thread for lookup, no ThreadPool
#     > 1: Uses a ThreadPool of nThreads
# - If too large (say, > 50), may run into server errors such as "HTTP Error 429: Too Many Requests"
nThreads = 50


# In[3]:


# Gene symbol lookup functions

def searchGeneNames(term, email, useSingleIndirectMatch = True):
    '''
    Search official human gene names and aliases in NCBI Gene database for a match to term, returning offical names and IDs.
    Dependencies: Biopython
    
    Args:
    - term: str
        gene name / alias
    - email: str
        email registered with NCBI
    - useSingleIndirectMatch
        If the Entrez Gene Database query returns only 1 NCBI Gene ID, even if the term does not exactly match the gene symbol
        or an alias, use the match.
    
    Returns: dict: str -> list
        "names": list of matched official gene name(s)
        "ids": list of NCBI Gene IDs corresponding to matched official gene name(s)
    '''

    names, ids = [], []
    Entrez.email = email
    handle = Entrez.esearch(db="gene", term='(' + term + '[gene]) AND (Homo sapiens[orgn]) AND alive[prop] NOT newentry[gene]')
    idList = Entrez.read(handle)['IdList']
    for id in idList:
        handle = Entrez.esummary(db='gene', id=id)
        record = Entrez.read(handle)
        name = record['DocumentSummarySet']['DocumentSummary'][0]['Name']
        aliases = record['DocumentSummarySet']['DocumentSummary'][0]['OtherAliases'].split(', ')
        if (term in [name] + aliases):
            names.append(name)
            ids.append(id)
    if useSingleIndirectMatch:
        if len(names) == 0 and len(idList) == 1:
            ids.append(idList[0])
            handle = Entrez.esummary(db='gene', id=id)
            record = Entrez.read(handle)
            names.append(record['DocumentSummarySet']['DocumentSummary'][0]['Name'])
    return({"names": names, "ids": ids})

def geneSymbolLookupFromSeries(series, email):
    '''
    Search official human gene names and aliases in NCBI Gene database for the value in 'Name' index of given pandas Series.
    
    Args
    - series: pandas.Series
        must have 'Name' index
    - email: str
        email registered with NCBI
    
    Returns: str
      If no match found, returns the empty string. Otherwise, returns the first matched official gene symbol.
    '''
    
    term = series['Name']
    match = searchGeneNames(term, email)
    if len(match["names"]) == 0:
        return("")
    if (term in match["names"]):
        return(term)
    return(match["names"][0])

def multiThreadedSearchGeneNames(terms, email, nThreads = None):
    '''
    Search official human gene names and aliases in NCBI Gene database for terms.
    
    Args
    - terms: list of str
        terms to lookup
    - email: str
        email registered with NCBI
    - nThreads: int
        None: Uses a ThreadPool of a default number of threads as returned by os.cpu_count()
        1+: Uses a ThreadPool of nThreads
    
    Returns: list of str
      Where no matches found, returns the empty string. Otherwise, returns the first matched official gene symbol.
    '''
    
    dict_results = []
    gene_symbols = []
    pool = ThreadPool(nThreads)
    print("Using {:d} threads...".format(pool._processes))
    for i in range(len(terms)):
        dict_results.append(pool.apply_async(searchGeneNames, (terms[i], email)))
    pool.close()
    pool.join()
    for i in range(len(terms)):
        match = dict_results[i].get()
        if len(match["names"]) == 0:
            gene_symbols.append("")
        elif terms[i] in match["names"]:
            gene_symbols.append(terms[i])
        else:
            gene_symbols.append(match["names"][0])
    return(gene_symbols)


# ## Process data from Klaeger et al.

# In[4]:


# Parameters
# - highConfidenceOnly: logical
#     Keep only "high confidence" protein-drug interactions.
#       "A protein was considered a high-confidence target if the binding curve showed a sigmoidal shape
#        with a dose-dependent decrease in binding to the Kinobeads." (Klaeger et al., Supplementary Materials)
# - na_values_klaeger: list of str
#     Values from CSV file to be read in as np.nan
# - dtype_klaeger: dict: str -> dtype
#     dtypes of columns of data

highConfidenceOnly = True
na_values_klaeger = ["n.d."]
dtype_klaeger = {'Name': 'str', 'STF1081': np.float64, 'CC401': np.float64}


# In[5]:


# Read in data
# - cannot set dtype of columns yet because of string "n.i." (not inhibited) values in STF1081 and CC401 columns,
#   which should be of dtype np.float
df1 = pd.read_csv(os.path.join(data_dir, Klaeger_filename), na_values=na_values_klaeger)


# In[6]:


# Remove rows (targets) where values are NA
df1.dropna(axis=0, how='any', inplace=True)


# In[7]:


# Remove low-confidence values if specified
if highConfidenceOnly:
    df1.drop(df1.index[df1['STF1081'].str.contains('\(') == True], axis=0, inplace=True)
    df1.drop(df1.index[df1['CC401'].str.contains('\(') == True], axis=0, inplace=True)
else:
    df1['STF1081'] = df1['STF1081'].str.strip('()')
    df1['CC401'] = df1['CC401'].str.strip('()')


# In[8]:


# Convert "n.i." (not inhibited) values to np.inf
df1.replace('n.i.', np.inf, inplace=True)
df1 = df1.astype(dtype_klaeger)


# In[9]:


# Split combined target names into individual target names with their own row
# - Ex: 
#             Name  CC401  STF1081          Name  CC401  STF1081
#  CSNK2A1;CSNK2A3  747.0    11.0  -->  CSNK2A1  747.0      11.0
#                                       CSNK2A3  747.0      11.0
for name in df1.loc[:,'Name']:
    names = name.split(';')
    if (len(names) > 1):
        series = df1.loc[df1.index[df1['Name'] == name]].squeeze()
        for sub_name in names:
            series['Name'] = sub_name
            df1 = df1.append(series, ignore_index=True)
        df1.drop(df1.index[df1['Name'] == name], axis=0, inplace=True)
df1.reset_index(drop=True, inplace=True)


# In[10]:


# Add official gene names as a new column
if (nThreads is None) or nThreads > 1:
    df1[geneSymbolColumn] = multiThreadedSearchGeneNames(df1['Name'].tolist(), EntrezEmail, nThreads)
else:
    df1[geneSymbolColumn] = df1.apply(geneSymbolLookupFromSeries, 1, email = EntrezEmail)

# Show rows where no official gene symbol was found
display(df1.loc[df1[geneSymbolColumn] == ""])


# In[12]:


# For genes with missing official gene symbols, manually add official gene symbols
df1.loc[df1['Name'] == 'Q6ZSR9', geneSymbolColumn] = 'Q6ZSR9'


# In[13]:


# Sort by official gene names
df1.sort_values(by=geneSymbolColumn, inplace=True)
df1.reset_index(drop=True, inplace=True)

# Reorder columns
df1 = df1[[geneSymbolColumn, "Name", "STF1081", "CC401"]]


# In[14]:


# Quick verification of gene symbol lookups

# Display rows where official gene symbol differed from original gene name
display(df1[df1['Name'] != df1[geneSymbolColumn]])

# Confirm there are no duplicate genes (rows)
print("Number of duplicated rows: " + str(sum(df1.duplicated() == True)))


# In[15]:


display(df1)


# In[16]:


df1.to_csv(os.path.join(dataAux_dir, Klaeger_filename), index=False)


# ## Process data from Huang et al.

# In[17]:


df2 = pd.read_csv(os.path.join(data_dir, Huang_filename))


# In[18]:


# Remove rows (targets) where values are NA
df2.dropna(axis=0, how='any', inplace=True)


# In[19]:


# Add official gene names as a new column
if (nThreads is None) or nThreads > 1:
    df2[geneSymbolColumn] = multiThreadedSearchGeneNames(df2['Name'].tolist(), EntrezEmail, nThreads)
else:
    df2[geneSymbolColumn] = df2.apply(geneSymbolLookupFromSeries, 1, email = EntrezEmail)

# Show rows where no official gene symbol was found
display(df2.loc[df2[geneSymbolColumn] == ""])


# In[20]:


# For genes with missing official gene symbols, manually add official gene symbols
df2.loc[df2['Name'] == 'p38 alpha', geneSymbolColumn] = 'MAPK14'
df2.loc[df2['Name'] == 'p38 beta', geneSymbolColumn] = 'MAPK11'
df2.loc[df2['Name'] == 'p38 gamma', geneSymbolColumn] = 'MAPK12'
df2.loc[df2['Name'] == 'p38 delta', geneSymbolColumn] = 'MAPK13'
df2.loc[df2['Name'] == 'PKB beta', geneSymbolColumn] = 'AKT2'
df2.loc[df2['Name'] == 'PKA', geneSymbolColumn] = 'PRKACA'
df2.loc[df2['Name'] == 'CAMKK beta', geneSymbolColumn] = 'CAMKK2'
df2.loc[df2['Name'] == 'GSK3 beta', geneSymbolColumn] = 'GSK3B'
df2.loc[df2['Name'] == 'CDK2-Cyclin A', geneSymbolColumn] = 'CDK2'
df2.loc[df2['Name'] == 'CDK9-Cyclin T1', geneSymbolColumn] = 'CDK9'
df2.loc[df2['Name'] == 'Aurora A', geneSymbolColumn] = 'AURKA'
df2.loc[df2['Name'] == 'Aurora B', geneSymbolColumn] = 'AURKB'
df2.loc[df2['Name'] == 'AMPK (hum)', geneSymbolColumn] = 'PRKAA1'
df2.loc[df2['Name'] == 'CK1 gamma 2', geneSymbolColumn] = 'CSNK1G2'
df2.loc[df2['Name'] == 'CK1 delta', geneSymbolColumn] = 'CSNK1D'
df2.loc[df2['Name'] == 'CK2', geneSymbolColumn] = 'CSNK2A1'
df2.loc[df2['Name'] == 'IKK epsilon', geneSymbolColumn] = 'IKBKE'
df2.loc[df2['Name'] == 'EF2K', geneSymbolColumn] = 'EEF2K' # based on UniProt
df2.loc[df2['Name'] == 'MPSK1', geneSymbolColumn] = 'STK16'
df2.loc[df2['Name'] == 'EPH-A2', geneSymbolColumn] = 'EPHA2'
df2.loc[df2['Name'] == 'EPH-A4', geneSymbolColumn] = 'EPHA4'
df2.loc[df2['Name'] == 'EPH-B1', geneSymbolColumn] = 'EPHB1'
df2.loc[df2['Name'] == 'EPH-B2', geneSymbolColumn] = 'EPHB2'
df2.loc[df2['Name'] == 'EPH-B3', geneSymbolColumn] = 'EPHB3'
df2.loc[df2['Name'] == 'EPH-B4', geneSymbolColumn] = 'EPHB4'
df2.loc[df2['Name'] == 'FGF-R1', geneSymbolColumn] = 'FGFR1'
df2.loc[df2['Name'] == 'IGF-1R', geneSymbolColumn] = 'IGF1R'
df2.loc[df2['Name'] == 'IR', geneSymbolColumn] = '' # unknown
df2.loc[df2['Name'] == 'PINK', geneSymbolColumn] = 'PINK1'


# In[21]:


# Sort by official gene names
df2.sort_values(by=geneSymbolColumn, inplace=True)
df2.reset_index(drop=True, inplace=True)

# Reorder columns
df2 = df2[[geneSymbolColumn, "Name", "STF1081", "HTH01091"]]


# In[22]:


# Quick verification of gene symbol lookups

# Display rows where official gene symbol differed from original gene name
with pd.option_context('display.max_rows', None):
    display(df2[df2['Name'] != df2[geneSymbolColumn]])

# Confirm there are no duplicate genes (rows)
print("Number of duplicated rows: " + str(sum(df2.duplicated() == True)))


# In[23]:


display(df2)


# In[24]:


df2.to_csv(os.path.join(dataAux_dir, Huang_filename), index=False)


# ## Process data from Annes et al.

# In[25]:


# Parameters
# - condense_func: function
#     function condense values of multiple variants (different phosphorylation states, mutants) of the same kinase
condense_func = np.mean


# In[26]:


df3 = pd.read_csv(os.path.join(data_dir, Annes100_filename))
df4 = pd.read_csv(os.path.join(data_dir, Annes500_filename))


# In[27]:


# Remove rows (targets) where values are NA
df3.dropna(axis=0, how='any', inplace=True)
df4.dropna(axis=0, how='any', inplace=True)


# In[28]:


def condenseDuplicatesByKey(df, keyCol, valueCol, func):
    '''
    Condense duplicate rows (based on the keyCol column) by applying a specified function to elements in the valueCol column.
    
    Args
    - df: pandas.DataFrame
    - keyCol: str
        column in which to look for duplicates
    - valueCol: str
        name of column of values to condense
    - func: function
        function to apply to duplicate values. Must take a list and return a single element
    
    Return: pandas.DataFrame
    '''
    for key in df[keyCol]:
        data = df[df[keyCol] == key]
        if data.shape[0] > 1:
            series = data.iloc[0,:].copy() # deep copy to avoid error of assigning value (next line) to a view of slice from df
            series[valueCol] = func(data[valueCol])
            df.drop(df.index[df[keyCol] == key], axis=0, inplace=True)
            df = df.append(series, ignore_index=True)
    return(df)


# In[29]:


df3 = condenseDuplicatesByKey(df3, "Name", "STF1081", condense_func)
df4 = condenseDuplicatesByKey(df4, "Name", "STF1081", condense_func)


# In[30]:


# Add official gene names as a new column
if (nThreads is None) or nThreads > 1:
    df3[geneSymbolColumn] = multiThreadedSearchGeneNames(df3['Name'].tolist(), EntrezEmail, nThreads)
    df4[geneSymbolColumn] = multiThreadedSearchGeneNames(df4['Name'].tolist(), EntrezEmail, nThreads)
else:
    df3[geneSymbolColumn] = df3.apply(geneSymbolLookupFromSeries, 1, email = EntrezEmail)
    df4[geneSymbolColumn] = df4.apply(geneSymbolLookupFromSeries, 1, email = EntrezEmail)

# Show rows where no official gene symbol was found
display(df3.loc[df3[geneSymbolColumn] == ""])
display(df4.loc[df4[geneSymbolColumn] == ""])


# In[31]:


# For genes with missing official gene symbols, manually add official gene symbols
geneSymbolsMap = {
    "MGC42105": "NIM1K",
    "CDPK1": "PF3D7_0217500", # Genus/species: Plasmodium falciparum 3D7 - https://www.ncbi.nlm.nih.gov/gene/812762
    "MAL13P1.279": "PF3D7_1356900", # Genus/species: Plasmodium falciparum 3D7 - https://www.ncbi.nlm.nih.gov/gene/813841
    "pknB": "pknB", # Genus/species: Mycobacterium tuberculosis H37Rv - https://www.ncbi.nlm.nih.gov/gene/887072
    "KIAA0999": "SIK3"}

df3.loc[df3["Name"].isin(geneSymbolsMap), geneSymbolColumn] = list(geneSymbolsMap.values())
df4.loc[df3["Name"].isin(geneSymbolsMap), geneSymbolColumn] = list(geneSymbolsMap.values())


# In[32]:


# Sort by official gene names
df3.sort_values(by=geneSymbolColumn, inplace=True)
df3.reset_index(drop=True, inplace=True)
df4.sort_values(by=geneSymbolColumn, inplace=True)
df4.reset_index(drop=True, inplace=True)

# Reorder columns
df3 = df3[[geneSymbolColumn, "Name", "STF1081", "STF1285"]]
df4 = df4[[geneSymbolColumn, "Name", "STF1081", "STF1285"]]


# In[33]:


# Quick verification of gene symbol lookups

# Display rows where official gene symbol differed from original gene name
display(df3[df3['Name'] != df3[geneSymbolColumn]])
display(df4[df4['Name'] != df4[geneSymbolColumn]])

# Confirm there are no duplicate genes (rows)
print("Number of duplicated rows: " + str(sum(df3.duplicated() == True)))
print("Number of duplicated rows: " + str(sum(df4.duplicated() == True)))


# In[34]:


display(df3)


# In[35]:


display(df4)


# In[36]:


df3.to_csv(os.path.join(dataAux_dir, Annes100_filename), index=False)
df4.to_csv(os.path.join(dataAux_dir, Annes500_filename), index=False)

