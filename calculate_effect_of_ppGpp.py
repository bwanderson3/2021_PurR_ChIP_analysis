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
chip_name_wt_00_1 = 'Run_2691/jawang/Sample_124312/Sample_124312_sampled.npy'
chip_wt_00_1 = np.load(chip_name_wt_00_1)
#%%
input_name_wt_00_1 = 'Run_2848/jawang/Sample_129017/Sample_129017_sampled.npy'
input_wt_00_1 = np.load(input_name_wt_00_1)

#%%
chip_name_wt_10_1 = 'Run_2691/jawang/Sample_124314/Sample_124314_sampled.npy'
chip_wt_10_1 = np.load(chip_name_wt_10_1)

#%%
input_name_wt_10_1 = 'Run_2848/jawang/Sample_129018/Sample_129018_sampled.npy'
input_wt_10_1 = np.load(input_name_wt_10_1)

#%%
chip_name_p0_00_1 = 'Run_2691/jawang/Sample_124315/Sample_124315_sampled.npy'
chip_p0_00_1 = np.load(chip_name_p0_00_1)

#%%
input_name_p0_00_1 = 'Run_2848/jawang/Sample_129021/Sample_129021_sampled.npy'
input_p0_00_1 = np.load(input_name_p0_00_1)

#%%
chip_name_p0_10_1 = 'Run_2691/jawang/Sample_124317/Sample_124317_sampled.npy'
chip_p0_10_1 = np.load(chip_name_p0_10_1)

#%%
input_name_p0_10_1 = 'Run_2848/jawang/Sample_129022/Sample_129022_sampled.npy'
input_p0_10_1 = np.load(input_name_p0_10_1)

#%%
chip_name_wt_00_2 = 'Run_2848/jawang/Sample_129005/Sample_129005_sampled.npy'
chip_wt_00_2 = np.load(chip_name_wt_00_2)

#%%
input_name_wt_00_2 = 'Run_2848/jawang/Sample_129015/Sample_129015_sampled.npy'
input_wt_00_2 = np.load(input_name_wt_00_2)

#%%
chip_name_wt_10_2 = 'Run_2848/jawang/Sample_129006/Sample_129006_sampled.npy'
chip_wt_10_2 = np.load(chip_name_wt_10_2)

#%%
input_name_wt_10_2 = 'Run_2848/jawang/Sample_129016/Sample_129016_sampled.npy'
input_wt_10_2 = np.load(input_name_wt_10_2)

#%%
chip_name_p0_00_2 = 'Run_2848/jawang/Sample_129009/Sample_129009_sampled.npy'
chip_p0_00_2 = np.load(chip_name_p0_00_2)

#%%
input_name_p0_00_2 = 'Run_2848/jawang/Sample_129019/Sample_129019_sampled.npy'
input_p0_00_2 = np.load(input_name_p0_00_2)

#%%
chip_name_p0_10_2  = 'Run_2848/jawang/Sample_129010/Sample_129010_sampled.npy'
chip_p0_10_2 = np.load(chip_name_p0_10_2)

#%%
input_name_p0_10_2 = 'Run_2848/jawang/Sample_129020/Sample_129020_sampled.npy'
input_p0_10_2 = np.load(input_name_p0_10_2)

#%%
chip_name_wt_00_3 = 'Run_2848/jawang/Sample_129007/Sample_129007_sampled.npy'
chip_wt_00_3 = np.load(chip_name_wt_00_3)

#%%
input_name_wt_00_3 = 'Run_2848/jawang/Sample_129017/Sample_129017_sampled.npy'
input_wt_00_3 = np.load(input_name_wt_00_3)

#%%
chip_name_wt_10_3 = 'Run_2848/jawang/Sample_129008/Sample_129008_sampled.npy'
chip_wt_10_3 = np.load(chip_name_wt_10_3)

#%%
input_name_wt_10_3 = 'Run_2848/jawang/Sample_129018/Sample_129018_sampled.npy'
input_wt_10_3 = np.load(input_name_wt_10_3)

#%%
chip_name_p0_00_3 = 'Run_2848/jawang/Sample_129011/Sample_129011_sampled.npy'
chip_p0_00_3 = np.load(chip_name_p0_00_3)

#%%
input_name_p0_00_3 = 'Run_2848/jawang/Sample_129021/Sample_129021_sampled.npy'
input_p0_00_3 = np.load(input_name_p0_00_3)

#%%
chip_name_p0_10_3  = 'Run_2848/jawang/Sample_129012/Sample_129012_sampled.npy'
chip_p0_10_3 = np.load(chip_name_p0_10_3)

#%%
input_name_p0_10_3 = 'Run_2848/jawang/Sample_129022/Sample_129022_sampled.npy'
input_p0_10_3 = np.load(input_name_p0_10_3)

#%%
enrich_wt_10_1 = np.log2(chip_wt_10_1/input_wt_10_1)
enrich_wt_00_1 = np.log2(chip_wt_00_1/input_wt_00_1)

enrich_p0_10_1 = np.log2(chip_p0_10_1/input_p0_10_1)
enrich_p0_00_1 = np.log2(chip_p0_00_1/input_p0_00_1)

#%%
enrich_wt_10_2 = np.log2(chip_wt_10_2/input_wt_10_2)
enrich_wt_00_2 = np.log2(chip_wt_00_2/input_wt_00_2)

enrich_p0_10_2 = np.log2(chip_p0_10_2/input_p0_10_2)
enrich_p0_00_2 = np.log2(chip_p0_00_2/input_p0_00_2)

#%%
enrich_wt_10_3 = np.log2(chip_wt_10_3/input_wt_10_3)
enrich_wt_00_3 = np.log2(chip_wt_00_3/input_wt_00_3)

enrich_p0_10_3 = np.log2(chip_p0_10_3/input_p0_10_3)
enrich_p0_00_3 = np.log2(chip_p0_00_3/input_p0_00_3)

#%%
# 
effect_of_p_1 = (enrich_wt_10_1 - enrich_wt_00_1) - (enrich_p0_10_1 - enrich_p0_00_1)

effect_of_p_2 = (enrich_wt_10_2 - enrich_wt_00_2) - (enrich_p0_10_2 - enrich_p0_00_2)

effect_of_p_3 = (enrich_wt_10_3 - enrich_wt_00_3) - (enrich_p0_10_3 - enrich_p0_00_3)

#%%
np.save('effect_of_ppGpp', effect_of_p)

#%%
# Take median row-wise. Check by doing med.shape to ensure you have 400,000 x 1 array
med_1 = np.median(effect_of_p_1, axis=1)
med_1.shape

med_2 = np.median(effect_of_p_2, axis=1)
med_2.shape

med_3 = np.median(effect_of_p_3, axis=1)
med_3.shape

#%% once you have all three replicates' medians
# play with axis to ensure you get correct shape
mean_vals = np.mean([med_1, med_2, med_3], axis=0)
mean_vals.shape

#%%
sd_vals = np.std([med_1, med_2, med_3], axis=0)
sd_vals.shape

#%%
np.save('mean_ppGpp_effect', mean_vals)
np.save('sd_ppGpp_effect', sd_vals)


# %%
