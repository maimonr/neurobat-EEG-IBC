#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pyxdf
import scipy.io
import numpy as np
import glob
import pathlib 
import os
import itertools


# In[17]:


wd = 'C:\\Recordings\\emotiv_testing\\20211210\\'
fnames = glob.glob(os.path.join(wd,'*.xdf'))


# In[18]:


emotiv_ids = ['E2020789','E20207BE']
exp_names = ['mutual_gaze','eyes_closed','visual_flicker_20hz','finger_tapping','metronome_180bpm','finger_tapping_metronome_180bpm']
sync_flag = True
fs = 128
mat_data ={}

for exp_name in exp_names:
    mat_data[exp_name] = []

for file_name in fnames:
    
    exp_name = file_name.split('\\')[-1].split('.')[0][:-2]
    exp_num = file_name.split('\\')[-1].split('_')[-1][0]
    
    streams,headers = pyxdf.load_xdf(file_name,synchronize_clocks=sync_flag,dejitter_timestamps=sync_flag)
    used_streams = []
    for e_id in emotiv_ids:
        emotiv_idx = [(e_id in s['info']['source_id'][0] and s['info']['type'][0] == 'EEG' and int(float(s['info']['nominal_srate'][0])) == fs) for s in streams]
        used_streams.append(streams[int(np.argwhere(emotiv_idx)[0])])
    
    exp_mat_data = []
    for s in used_streams:
        chNames = [ch['label'] for ch in s['info']['desc'][0]['channels'][0]['channel']]
        exp_mat_data.append({'data':s['time_series'],'timestamps':s['time_stamps'],'chNames':chNames,'source_id':s['info']['source_id'],'fs':fs})
    
    mat_data[exp_name].append(exp_mat_data)

mat_data['ids'] = emotiv_ids
scipy.io.savemat('C:/Users/BatLab/eeg-notebooks/emotiv_data_20211210.mat',mat_data)


# In[ ]:




