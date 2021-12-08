#%%
import os
import numpy as np
import pandas as pd
import glob

from matplotlib import pyplot as plt

#%%
analysis_direc = '/home/wanglab/src/PurR_ChIP_analysis'

#pre_inv_sample_ids = [line.strip() for line in open('pre_inversion_samples.txt')]
#post_inv_sample_ids = [line.strip() for line in open('post_inversion_samples.txt')]

data_direc = '/mnt/jadelab/lab/current/Users/Brent/ChIP-seq/PurR_ChIP_data/'
Run_2691_samples = [os.path.join(data_direc, 'Run_2691/jawang', direc) for direc in os.listdir('/mnt/jadelab/lab/current/Users/Brent/ChIP-seq/PurR_ChIP_data/Run_2691/jawang') if direc.startswith('Sample_')]
Run_2848_samples = [os.path.join(data_direc, 'Run_2848/jawang', direc) for direc in os.listdir('/mnt/jadelab/lab/current/Users/Brent/ChIP-seq/PurR_ChIP_data/Run_2848/jawang') if direc.startswith('Sample_')]

data_direcs = Run_2691_samples + Run_2848_samples

#%%
coverage_dict = {}

for direc in data_direcs:

    os.chdir(direc)
    sample_id = direc.split('/')[-1]

    inFileName = "{}_summary_median.npy".format(sample_id)

    if not os.path.isfile(inFileName):
        continue

    coverage = np.load(inFileName)

    coverage_dict[sample_id] = coverage

#%%
os.chdir(analysis_direc)

coverage_df = pd.DataFrame(coverage_dict)

coverage_data = coverage_df.values.astype('float32')

#%%
windowLen = 250
win = np.hanning(windowLen)

plt.plot(win)

#%%
def smooth(x,window_len=11,window='hanning'):
    
    if window_len<3:
        return(x)

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='same')[window_len-1:-1*(window_len-1)]
    return(y)

#%%
smoothed_data = np.apply_along_axis(smooth, 
                                    0, 
                                    coverage_data, 
                                    window_len=10)

#%%
smoothed_df = pd.DataFrame(smoothed_data, columns=coverage_df.columns)




#%%
def plotLocus(smoothed, original, start, end, sample_index=0):

    plt.plot(smoothed[start:end, sample_index], 'r-', linewidth=4, alpha=0.5)
    plt.plot(original[start:end, sample_index], linewidth=4, alpha=0.5)
    plt.show();


#%%
plotLocus(smoothed_data, coverage_data, start=5300, end=5600, sample_index=0)


#%%
print(smoothed_data.shape)
print(coverage_data.shape)


#%%
print(smoothed_data)


#%%
smoothed_df.to_csv(os.path.join(data_direc, "smoothed_data_100bp.csv"))
