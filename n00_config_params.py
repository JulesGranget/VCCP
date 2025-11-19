
import numpy as np


################################
######## GENERAL PARAMS ######## 
################################

srate = 500


sujet_list = ['NS217',
              ]


conditions = ['AoRoCo', 'AoR-Co', 'A+RoCc', 'A+R-Cc', 'A+RoC-', 'A+R-C-']





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

    path_main_workdir = '/home/jules/Documents/VCCL_JULES/Scripts'
    path_general = '/home/jules/Documents/VCCL_JULES'
    n_core = 25

    
path_data = os.path.join(path_general, 'Data')
path_prep = os.path.join(path_general, 'Analyses', 'preprocessing')
path_precompute = os.path.join(path_general, 'Analyses', 'precompute') 
path_results = os.path.join(path_general, 'Analyses', 'results') 

os.chdir(init_workdir)




################################
######## STRETCH ########
################################

stretch_TF_auto = False
ratio_stretch_TF = 0.5




