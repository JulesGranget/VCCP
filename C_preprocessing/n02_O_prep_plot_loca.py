



from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *

import os
import numpy as np
import matplotlib.pyplot as plt
import glob
import pandas as pd


debug = False

def generate_plot_loca(sujet):

    path_open_chanlist_raw = os.path.join(path_prep, sujet, f"{sujet}_chanlist.xlsx")
    path_open_correspondance = os.path.join(path_data, sujet, 'anatomy', f"{sujet}_Electrodes_Natus_TDT_correspondence.xlsx")
    
    df_chanlist_raw = pd.read_excel(path_open_chanlist_raw)
    df_correspondance = pd.read_excel(path_open_correspondance)

    auxchan = auxchan_allsujet[sujet]
    chanlist_raw_noaux = df_chanlist_raw.query(f"chan not in {auxchan}").values.reshape(-1)

    #### identify missplots
    chanlist_not_in_correspondance = []
    chanlist_in_correspondance = []
    for chan_i, chan in enumerate(chanlist_raw_noaux):
        if chan in df_correspondance['Label'].values:
            chanlist_in_correspondance.append(chan)
        else:
            chanlist_not_in_correspondance.append(chan)  

    #### export chanlist
    path_export_chanlists = os.path.join(path_prep, sujet, 'anatomy')

    chanlist_in_correspondance_textfile = open(os.path.join(path_export_chanlists, f"{sujet}_chanlist_in_correspondance.txt"), "w")
    for chan in chanlist_in_correspondance:
        chanlist_in_correspondance_textfile.write(chan + "\n")
    chanlist_in_correspondance_textfile.close()

    chanlist_not_in_correspondance_textfile = open(os.path.join(path_export_chanlists, f"{sujet}_chanlist_not_in_correspondance.txt"), "w")
    for chan in chanlist_not_in_correspondance:
        chanlist_not_in_correspondance_textfile.write(chan + "\n")
    chanlist_not_in_correspondance_textfile.close()

    #### construct plotloca
    df_plot_loca = []			

    for chan in chanlist_in_correspondance:

        _df_sel = df_correspondance.query(f"Label == '{chan}'")

        _df_plot_loca = pd.DataFrame({'chan' : [chan], 'FS_volumetric' : [_df_sel['FS_vol'].values[0]], 'side' : [0], 'type' : [0], 'SOZ' : [0], 'Spikey_or_bad' : [0], 'Out' : [0], 
                                      'sig_inspected' : [0], 'loca_inspected' : [0], 'double_inspection' : [0], 'select' : [1], 'FS_volumetric_corrected' : [_df_sel['FS_vol'].values[0]], 
                                      'LEPTO_coords_1' : [_df_sel['LEPTO_coords_1'].values[0]], 'LEPTO_coords_2' : [_df_sel['LEPTO_coords_2'].values[0]], 
                                      'LEPTO_coords_3' : [_df_sel['LEPTO_coords_3'].values[0]], 'fsaverage_coords_1' : [_df_sel['fsaverage_coords_1'].values[0]], 
                                        'fsaverage_coords_2' : [_df_sel['fsaverage_coords_2'].values[0]], 'fsaverage_coords_3' : [_df_sel['fsaverage_coords_3'].values[0]]})
        
        df_plot_loca.append(_df_plot_loca)
        
    df_plot_loca = pd.concat(df_plot_loca)
    
    #### correct localist
    loca_list = df_plot_loca['FS_volumetric'].values
    loca_list_modified = modify_loca_name(loca_list)
    df_plot_loca['FS_volumetric'] = loca_list_modified[:,0]
    df_plot_loca['FS_volumetric_corrected'] = loca_list_modified[:,0]
    df_plot_loca['side'] = loca_list_modified[:,2]
    df_plot_loca['type'] = loca_list_modified[:,1]

    #### export
    path_export_plotloca = os.path.join(path_prep, sujet, 'anatomy')
    df_plot_loca.to_excel(os.path.join(path_export_plotloca, f"{sujet}_plot_loca.xlsx"))






################################
######## EXECUTE ########
################################


if __name__== '__main__':

    #sujet = sujet_list[0]
    for sujet in sujet_list:

        print(f'#### {sujet}')

        #### execute
        if os.path.exists(os.path.join(path_prep, sujet, 'anatomy', f"{sujet}_plot_loca.xlsx")):
            print('#### MONOPOL ALREADY COMPUTED ####')
        else:
            generate_plot_loca(sujet)








