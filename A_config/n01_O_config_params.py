


import scipy
import numpy as np





################################
######## GENERAL PARAMS ######## 
################################

srate = 500

sujet_pilot = ['NS217']

sujet_list = ['LH018']


conditions = ['RB', 'HV', 'AoRoCo', 'AoR-Co', 'A+RoCc', 'A+R-Cc', 'A+RoC-', 'A+R-C-']
conditions_VCCP = ['AoRoCo', 'AoR-Co', 'A+RoCc', 'A+R-Cc', 'A+RoC-', 'A+R-C-']

aux_chanlist = ['resp', 'co2', 'gsm']

correspondance_cond = {'baseline' : 'RB', 'HV' : 'HV', '(1/2)A0(6)R-C0' : 'AoR-Co', '(1)A+(6)R-Cc' : 'A+R-Cc', '(2)A+(1)R0C-' : 'A+RoC-', '(2)A+(1)R0Cc' : 'A+RoCc', 
                       '(2)A+(1)R0Cc' : 'A+RoCc', '(1)A0(1)R0C0' : 'AoRoCo', '(1)A+(6)R-C-' : 'A+R-C-'}

freq_band_dict = {'theta' : [4,8], 'alpha' : [8,12], 'beta' : [12,50], 'gamma' : [60, 150]}


ROI_short_list = ['entorhinal', 'Amygdala', 'medialorbitofrontal', 'Hippocampus', 
       'lateralorbitofrontal', 'Thalamus', 'posteriorcingulate',
       'isthmuscingulate']





########################################
######## PATH DEFINITION ########
########################################

import socket
import os
import platform
 
PC_OS = platform.system()
PC_ID = socket.gethostname()
init_workdir = os.getcwd()

if PC_ID == 'jules-ubuntu1':

    path_main_workdir = '/home/jules/Documents/VCCP_JULES/Scripts'
    path_general = '/home/jules/Documents/VCCP_JULES'
    n_core = 25

    
path_data = os.path.join(path_general, 'Data')
path_prep = os.path.join(path_general, 'Analyses', 'preprocessing')
path_anatomy = os.path.join(path_general, 'Data', 'anatomy')
path_precompute = os.path.join(path_general, 'Analyses', 'precompute') 
path_results = os.path.join(path_general, 'Analyses', 'results') 
path_memmap = os.path.join(path_general, 'memmap') 

os.chdir(init_workdir)




################################
######## PRECOMPUTE TF ########
################################

#### chunk data
chunk_time = 2 #sec
exclude_cycle_duration_jules = 2 #sec
exclude_inspi_duration_jules = 1 #sec

#### stretch
stretch_point_TF = 250
stretch_TF_auto = False
ratio_stretch_TF = 0.50

#### TF & ITPC
nfrex = 150
ncycle_list = [7, 41]
freq_list = [2, 150]
srate_dw = 10
wavetime = np.arange(-3,3,1/srate)
frex = np.logspace(np.log10(freq_list[0]), np.log10(freq_list[1]), nfrex) 
cycles = np.logspace(np.log10(ncycle_list[0]), np.log10(ncycle_list[1]), nfrex).astype('int')
Pxx_wavelet_norm = 1000

#### STATS
n_surrogates_tf = 1000
tf_percentile_sel_stats = 2 # for both side
tf_stats_percentile_cluster = 95
norm_method = 'rscore'# 'zscore', 'dB'
tf_stats_percentile_cluster_allplot = 99
df_extraction_Cxy = 0.02 #Hz around respi median


#### plot
tf_plot_percentile_scale = 1 #for one side




################################
######## PRECOMPUTE Cxy ########
################################


zero_pad_coeff = 15

def get_params_spectral_analysis(srate):

    nwind = int( 20*srate ) # window length in seconds*srate
    nfft = nwind*zero_pad_coeff # if no zero padding nfft = nwind
    noverlap = np.round(nwind/2) # number of points of overlap here 50%
    hannw = scipy.signal.windows.hann(nwind) # hann window

    return nwind, nfft, noverlap, hannw

#### plot Pxx Cxy  
if zero_pad_coeff - 5 <= 0:
    remove_zero_pad = 0
remove_zero_pad = zero_pad_coeff - 5

#### coh
freq_surrogates = [0, 2]
Cxy_extraction_range = 0.025











