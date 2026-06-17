
from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *

import plotly
import plotly.graph_objects as go
import h5py



debug = False









########################
######## EXPORT ########
########################

#sujet = sujet_list[0]
def export_data(sujet):

    #### get data
    df_loca = get_df_loca(sujet)
    path_chanlist_raw = os.path.join(path_prep, sujet, f"{sujet}_chanlist.xlsx")
    chanlist_raw = pd.read_excel(path_chanlist_raw)['chan'].values

    path_data_raw = os.path.join(path_prep, sujet, f"{sujet}_data_prep.mat")
    f = h5py.File(path_data_raw, 'r')
    #print(f.keys())
    data_prep = f['eeg_export'][:]

    ########################
    dset = f['eeg_export']

    bad_chunks = []

    for start in range(0, dset.shape[0], dset.chunks[0]):
        stop = min(start + dset.chunks[0], dset.shape[0])
        try:
            _ = dset[start:stop, :]
        except Exception as e:
            print("Bad chunk:", start, stop, e)
            bad_chunks.append((start, stop))

    print("Number bad chunks:", len(bad_chunks))
    ########################

    #### filter good chan
    select_chan = df_loca.query(f"select == 1")['chan'].values
    select_chan_i = np.array([np.where(chanlist_raw == chan)[0][0] for chan in select_chan])
    select_chan_aux_i = np.array([np.where(chanlist_raw == chan)[0][0] for chan in auxchan_allsujet[sujet]])
    select_allchan_i = np.concat([select_chan_i, select_chan_aux_i])
    select_chanlist_tot = np.concat([select_chan, np.array(auxchan_allsujet[sujet])])

    data_filtered = data_prep[select_allchan_i]

    #### chunk conditions & export
    path_trigger = os.path.join(path_prep, sujet, f"{sujet}_trigger.xlsx")
    trial_info = pd.read_excel(path_trigger)
    trial_info['cond'] = [correspondance_cond[_cond] for _cond in trial_info['cond'].values]

    count_cond = {cond : 0 for cond in conditions}
    count_trial_list = []
    for trial in trial_info['cond']:

        count_cond[trial] += 1
        count_trial_list.append(f"{trial}_{count_cond[trial]}")
    
    trial_info['cond'] = count_trial_list

    path_export_trial = os.path.join(path_prep, sujet, 'trial_exports')

    #trial = 'RB_1'
    for trial in trial_info['cond']:
        
        _trial_i = trial_info.query(f"cond == '{trial}'")
        start, stop = _trial_i['start'][0], _trial_i['stop'][0]

        data_trial = data_filtered[:,start:stop]

        coords_xr_data_trial = {'chan' : select_chanlist_tot, 'time' : np.arange(data_trial.size)/srate}
        xr_data_trial = xr.DataArray(data_trial, dims=coords_xr_data_trial.keys(), coords=coords_xr_data_trial)

        xr_data_trial.to_netcdf(os.path.join(path_export_trial, f"{sujet}_{trial}.nc"))


    

            







################################
######## EXECUTE ########
################################

if __name__ == '__main__':

    sujet = 'NS217'
    sujet = 'LH018'

    export_data(sujet)



