



from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *


debug = False




########################
######## Cxy ########
########################



#sujet = sujet_list[0]
def precompute_Cxy_sujet(sujet):

    path_export_Cxy = os.path.join(path_precompute, 'Cxy')
    if os.path.exists(os.path.join(path_export_Cxy, f"{sujet}_surr_Cxy.xr")):
        print(f'Cxy ALREADY COMPUTED {sujet}')
        return
    else:
        print(f"COMPUTE Cxy {sujet}")
    
    #### config
    nwind, nfft, noverlap, hannw = get_params_spectral_analysis(srate)
    chanlist, chanlist_aux, chanlist_all, localist = get_chanlist_localist(sujet)

    #### precompute frequency bins and mask
    hzCxy = np.linspace(0, srate / 2, int(nfft / 2 + 1))
    mask_hzCxy = (hzCxy >= freq_surrogates[0]) & (hzCxy < freq_surrogates[1])
    hzCxy_respi = hzCxy[mask_hzCxy]

    #### extract Cxy
    trial_list = get_alltrials_sujet(sujet)

    Cxy_sujet_allcond = []

    #cond = conditions[-1]
    for cond in conditions:

        trial_list_cond = [trial for trial in trial_list if trial.find(cond) != -1]

        Pxx_trial_list = []

        for trial in trial_list_cond:

            print(f"{trial}")

            #### load
            data_raw = load_trial_data_sujet(sujet, trial)
            data = data_raw.sel(chan=chanlist)
            respi = data_raw.sel(chan='respi').values

            #### compute Cxy
            Cxy_chan = []

            for chan in chanlist:

                hzCxy_tmp, Cxy = scipy.signal.coherence(data.sel(chan=chan).values, respi, fs=srate, window=hannw, nperseg=nwind, noverlap=noverlap, nfft=nfft)
                Cxy_chan.append(Cxy[mask_hzCxy])

            Cxy_chan = np.vstack(Cxy_chan)
            Pxx_trial_list.append(Cxy_chan)

        if len(trial_list_cond) != 1:

            Pxx_cond = np.median(np.stack(Pxx_trial_list), axis=0)

        else:

            Pxx_cond = Pxx_trial_list[0]

        Cxy_sujet_allcond.append(Pxx_cond)

    #### get xr
    xr_coords = {'cond' : conditions, 'chan' : chanlist, 'freq' : hzCxy_respi}
    Cxy_sujet_allcond = np.stack(Cxy_sujet_allcond)
    xr_Cxy_sujet = xr.DataArray(data=Cxy_sujet_allcond, dims=xr_coords.keys(), coords=xr_coords)

    #### export
    xr_Cxy_sujet.to_netcdf(os.path.join(path_export_Cxy, f"{sujet}_Cxy.nc"))

    print('Done')









################################
######## EXECUTE ########
################################


if __name__ == '__main__':


    for sujet in sujet_list:

        precompute_Cxy_sujet(sujet)

    






