#%%
import os
import numpy as np
import pandas as pd

#%%
data_direc = "Y:/lab/current/Users/Brent/ChIP-seq/PurR_ChIP_data"
analysis_direc = "C:/Users/civwa/Documents/GitHub/PurR_ChIP_analysis"

#%%
os.chdir(data_direc)
info_fname = 'PurR_chip_and_input_pairs.csv'

info_df = pd.read_csv(info_fname)
print(info_df)
#%%
chip_name_wt_00 = 'Run_2691/jawang/Sample_124312/Sample_124312_sampled.npy'
chip_wt_00 = np.load(chip_name_wt_00)
#%%
input_name_wt_00 = 'Run_2848/jawang/Sample_129017/Sample_129017_sampled.npy'
input_wt_00 = np.load(input_name_wt_00)

#%%
chip_name_wt_10 = 'Run_2691/jawang/Sample_124312/Sample_124312_sampled.npy'
chip_wt_10 = np.load(chip_name_wt_10)

#%%
input_name_wt_10 = 'Run_2848/jawang/Sample_129018/Sample_129018_sampled.npy'
input_wt_10 = np.load(input_name_wt_10)

#%%
chip_name_p0_10 = 'Run_2691/jawang/Sample_124317/Sample_124317_sampled.npy'
chip_p0_10 = np.load(chip_name_p0_10)

#%%
input_name_p0_10 = 'Run_2848/jawang/Sample_129020/Sample_129020_sampled.npy'
input_p0_10 = np.load(input_name_p0_10)

#%%
chip_name_p0_00 = 'Run_2691/jawang/Sample_124315/Sample_124315_sampled.npy'
chip_p0_00 = np.load(chip_name_p0_00)

#%%
input_name_p0_00 = 'Run_2848/jawang/Sample_129019/Sample_129019_sampled.npy'
input_p0_00 = np.load(input_name_p0_00)

#%%
enrich_wt_10 = np.log2(chip_wt_10/input_wt_10)
enrich_wt_00 = np.log2(chip_wt_00/input_wt_00)

enrich_p0_10 = np.log2(chip_p0_10/input_p0_10)
enrich_p0_00 = np.log2(chip_p0_00/input_p0_00)

#%%
# 
effect_of_p = (enrich_wt_10 - enrich_wt_00) - (enrich_p0_10 - enrich_p0_00)

#%%
np.save('effect_of_ppGPp', effect_of_p)

#%%
# Take median row-wise??? Might need to be axis=1. Check by doing med.shape to ensure you have 400,000 x 1 array
med = np.median(effect_of_p, axis=0)
med.shape

#%% once you have all three replicate's medians
# play with axis to ensure you get correct shape
mean_vals = np.mean([rep1, rep2, rep3], axis=0)
mean_vals.shape
np.save('mean_effect', mean_vals)

sd_vals = np.sd([rep1, rep2, rep3], axis=0)
