
# coding: utf-8

# In[1]:


# Import all packages
import os
import numpy as np
import pandas as pd
from IPython.display import display


# In[2]:


# Global parameters / values
dataAux_dir = "../data_aux/"
results_dir = "../results/"
Klaeger_filename = "Klaeger.csv"
Huang_filename = "Huang.csv"
Annes100_filename = "Annes100.csv"
Annes500_filename = "Annes500.csv"

# Name of column of official gene symbols
geneSymbolColumn = "GeneSymbol"

# Thresholds for converting kinase activity metrics to boolean values
# - use_diff: bool
#     Use difference between STF1081 and another compound as value to threshold
#       For Annes and Huang datasets: use the absolute difference and compare to threshold diff_percent_thresh
#       For Klaeger dataset: use fold difference and compare to threshold diff_fold_thresh
# - diff_percent_thresh: int or float
#     Threshold absolute difference (% control (Annes) or % activity remaining (Huang)) between a less toxic compound
#     and STF1081 at which a target (kinase) is considered a potential target for the toxicity of STF1081
# - diff_fold_thresh: int or float
#     Threshold fold difference (Kd_app) between a less toxic compound and STF1081 at which a target (kinase)
#     is considered a potential target for the toxicity of STF1081
# - min_percent_thresh: int or float
#     Threshold % control (Annes) or % activity remaining (Huang) at which a target is considered to be inhibited by STF1081
# - max_percent_thresh: int or float
#     Threshold % control (Annes) or % activity remaining (Huang) at which a target is considered to be not inhibited by
#     a less toxic compound
use_diff = True
diff_percent_thresh = 20
diff_fold_thresh = 20
min_percent_thresh = 25
max_percent_thresh = 75

# Compare 100nM STF1081 to 100nM STF1285 (Annes100_filename) or 500nM STF1285 (Annes500_filename)
Annes_filename = Annes100_filename


# In[3]:


## Read in data
df1 = pd.read_csv(os.path.join(dataAux_dir, Klaeger_filename))
df2 = pd.read_csv(os.path.join(dataAux_dir, Huang_filename))
df3 = pd.read_csv(os.path.join(dataAux_dir, Annes_filename))


# ## Boolean filtering
# 
# Add new column `bool` to each DataFrame that indicates whether target is both inhibited by STF1081 and *not* inhibited by a less toxic compound. Find the intersection across all datasets.

# In[4]:


df1['bool'] = np.nan
df2['bool'] = np.nan
df3['bool'] = np.nan


# In[5]:


if use_diff:
    df1['bool'] = df1['CC401'] / df1['STF1081'] >= diff_fold_thresh
    df2['bool'] = df2['HTH01091'] - df2['STF1081'] >= diff_percent_thresh
    df3['bool'] = df3['STF1285'] - df3['STF1081'] >= diff_percent_thresh
else:
    df1.loc[np.isinf(df1['CC401']) & (df1['STF1081'] < np.inf), 'bool'] = True
    df2.loc[(df2['HTH01091'] >= max_percent_thresh) & (df2['STF1081'] <= min_percent_thresh), 'bool'] = True
    df3.loc[(df3['STF1285'] >= max_percent_thresh) & (df3['STF1081'] <= min_percent_thresh), 'bool'] = True


# In[6]:


intersect = pd.merge(df1, df2, how='inner', on=[geneSymbolColumn, 'bool'])
intersect = pd.merge(intersect, df3, how='inner', on=[geneSymbolColumn, 'bool'])
intersect = intersect.loc[intersect["bool"] == True, geneSymbolColumn]
intersect.sort_values(inplace=True)
intersect.reset_index(drop=True, inplace=True)
display(intersect)


# In[7]:


with open(os.path.join(results_dir, "intersect.txt"), "w") as f:
    f.write("\n".join(list(intersect)))
    f.write("\n")


# ## Rank ordering

# In[8]:


df1['diff'] = df1['CC401'] - df1['STF1081']
df1.sort_values(by=['diff','STF1081'], ascending=[False,True], inplace=True)
df1.reset_index(drop=True, inplace=True)
df1['rank'] = df1.index / len(df1)

df2['diff'] = df2['HTH01091'] - df2['STF1081']
df2.sort_values(by=['diff','STF1081'], ascending=[False,True], inplace=True)
df2.reset_index(drop=True, inplace=True)
df2['rank'] = df2.index / len(df2)

df3['diff'] = df3['STF1285'] - df3['STF1081']
df3.sort_values(by=['diff','STF1081'], ascending=[False,True], inplace=True)
df3.reset_index(drop=True, inplace=True)
df3['rank'] = df3.index / len(df3)


# In[9]:


targets = list(set(df1[geneSymbolColumn]) & set(df2[geneSymbolColumn]) & set(df3[geneSymbolColumn]))


# In[10]:


rank = pd.Series({target: (df1.loc[df1[geneSymbolColumn] == target, 'rank'].values[0] +
                           df2.loc[df2[geneSymbolColumn] == target, 'rank'].values[0] +
                           df3.loc[df3[geneSymbolColumn] == target, 'rank'].values[0]) for target in targets})
rank.sort_values(ascending=True, inplace=True)
with pd.option_context('display.max_rows', None):
    display(rank)


# In[11]:


rank.to_csv(os.path.join(results_dir, "rank.tsv"), index=True, sep="\t")


# ## Boolean filtering based on increasing STF-1285 concentration
# 
# This is not meant to identify toxic targets of STF1081 but instead to get an idea of targets that may be responsible for toxicity of STF-1285 at higher concentrations.

# In[12]:


# Parameters
diff_percent_thresh = 50
min_percent_thresh = 25
max_percent_thresh = 75


# In[13]:


df4 = pd.read_csv(os.path.join(dataAux_dir, Annes100_filename))
df5 = pd.read_csv(os.path.join(dataAux_dir, Annes500_filename))


# In[14]:


# merge 100 nM and 500 nM datasets
stf1285 = pd.merge(df4, df5, how='inner', on=[geneSymbolColumn], suffixes=(" (100 nM)", " (500 nM)"))
stf1285 = stf1285[[geneSymbolColumn, "STF1285 (100 nM)", "STF1285 (500 nM)"]]
stf1285['diff'] = stf1285["STF1285 (100 nM)"] - stf1285["STF1285 (500 nM)"]
stf1285.sort_values(by='diff', ascending=False, inplace=True)
stf1285.reset_index(drop=True, inplace=True)


# In[15]:


stf1285['bool'] = (stf1285['diff'] >= diff_percent_thresh) #& \
                  #(stf1285['STF1285 (500 nM)'] <= min_percent_thresh) & \
                  #(stf1285['STF1285 (100 nM)'] >= max_percent_thresh)
print(stf1285.loc[stf1285['bool'] == True, geneSymbolColumn].sort_values().reset_index(drop=True))

