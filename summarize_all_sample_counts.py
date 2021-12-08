import numpy as np 
import pandas as pd 
import os 
import subprocess
from dfply import *

def extract_dir_str(series_element):
    return series_element

data_direc = "/mnt/jadelab/lab/current/Users/Brent/ChIP-seq/PurR_ChIP_data"
info_file_name = "PurR_chip_and_input_pairs.csv"

info_file = os.path.join(data_direc, info_file_name)

info_df = pd.read_csv(info_file)

path_list = []
for index, row in info_df.iterrows():
    # print(data_direc)
    # print(row['Dir'])
    # print(row['Sample'])
    path_list.append(os.path.join(
        data_direc, 
        str(row['Dir']), 
        'jawang', 
        'Sample_{}'.format(str(row['Sample'])), 
        'Sample_{}_sampled.npy'.format(str(row['Sample']))
        ))

info_df = (info_df >> mutate(path=path_list))

for index, row in info_df.iterrows():
    elements = str(row['path']).split('/')[:-1]
    direc_path = '/'.join(elements)
    os.chdir(direc_path)
    # comm = '''python3 /home/wanglab/src/2018_Lrp_ChIP/ChIP_analysis/bootstrap_sam_file.py summarize\
    #             --samples_path {} --out_prefix Sample_{}_summary
    #         '''.format(str(row['path']), str(row['Sample']))
    p = subprocess.Popen(["python3", 
                        "/home/wanglab/src/2018_Lrp_ChIP/ChIP_analysis/bootstrap_sam_file.py",
                        "summarize",
                        "--samples_path",
                        str(row['path']),
                        "--out_prefix",
                        "Sample_{}_summary".format(str(row['Sample']))])
    p.communicate()
