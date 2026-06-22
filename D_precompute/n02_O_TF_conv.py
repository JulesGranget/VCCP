

import joblib
import plotly

from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *








################################
######## PRECOMPUTE TF ########
################################

#sujet, norm_param = sujet_list[0], 'rscore'
def precompute_tf(sujet, norm_param):

    #### get trial names
    trial_list = get_alltrials_sujet(sujet)

    #### precompute
    #trial = trial_list[5]
    for trial in trial_list:
        
        #### config
        path_export_TF = os.path.join(path_precompute, 'TF', sujet)
        path_load_trial = os.path.join(path_prep, sujet, 'trial_exports')

        if os.path.exists(os.path.join(path_export_TF, f"{sujet}_{trial}_TF_conv.npy")):
            print(f'{sujet} {trial} ALREADY COMPUTED', flush=True)
            return
        else:
            print(f'TF PRECOMPUTE {sujet} {trial}', flush=True)

        trial_name = f"{sujet}_{trial}.nc"

        xr_data = xr.load_dataarray(os.path.join(path_load_trial, trial_name))
        chanlist, chanlist_aux, chanlist_all, localist = get_chanlist_localist(sujet)
        resp = xr_data.sel(chan='respi').data
        CO2 = xr_data.sel(chan='CO2').data
        
        #### select wavelet parameters
        wavelets = get_wavelets()

        #### compute
        print('CONV', flush=True)

        tf_allconv = np.memmap(os.path.join(path_memmap, f'memmap_{sujet}_{trial}_tf_conv.npy'), mode="w+", shape=(len(chanlist), nfrex, xr_data.shape[-1]), dtype=np.float32)
                
        #chan_i, chan = 0, chanlist[0]
        # for chan_i, _ in enumerate(chanlist):
        def extract_chan_conv(chan_i,chan):

            print_advancement(chan_i, len(chanlist), steps=[25, 50, 75])

            for fi in range(nfrex):
                
                tf_allconv[chan_i, fi] = abs(scipy.signal.fftconvolve(xr_data.sel(chan=chan).values, wavelets[fi], 'same'))**2

        joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(extract_chan_conv)(chan_i,chan) for chan_i, chan in enumerate(chanlist))

        #### extract and norm
        print(f'NORM PARAMS', flush=True)

        if norm_param == 'rscore':

            rscore_param = np.memmap(os.path.join(path_memmap, f'memmap_{sujet}_{trial}_rscore_param.npy'), mode="w+", shape=(2, len(chanlist), nfrex), dtype=np.float32)

            # for chan_i, chan in enumerate(chanlist):
            def extract_chan_rscore_params(chan_i): 

                print_advancement(chan_i, chanlist.shape[0], steps=[25, 50, 75])

                for fi in range(nfrex):
                    
                    _med = np.median(tf_allconv[chan_i, fi,:])
                    _mad = scipy.stats.median_abs_deviation(tf_allconv[chan_i, fi,:].reshape(-1), axis=0)

                    if debug:

                        _cycles = tf_allconv[chan_i, fi,:]

                        for cycle_i in range(_cycles.shape[0]):
                            plt.plot(_cycles[cycle_i])
                        plt.show()

                        min, max = np.percentile(_cycles, 1), np.percentile(_cycles, 99)
                        plt.pcolormesh(_cycles, vmin=min, vmax=max)
                        plt.colorbar()
                        plt.show()

                    rscore_param[0,chan_i,fi] = _mad
                    rscore_param[1,chan_i,fi] = _med

            joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(extract_chan_rscore_params)(chan_i) for chan_i, _ in enumerate(chanlist))

            print(f'NORM', flush=True)

            tf_allconv_norm = np.memmap(os.path.join(path_memmap, f'memmap_{sujet}_{trial}_tf_conv_norm.npy'), mode="w+", shape=(len(chanlist), nfrex, xr_data.shape[-1]), dtype=np.float32)

            def norm_chan_conv(chan_i):

                print_advancement(chan_i, len(chanlist), steps=[25, 50, 75])

                mat = tf_allconv[chan_i]

                mad = rscore_param[0, chan_i, :]
                med = rscore_param[1, chan_i, :]

                mat_norm = (mat - med[:, None]) * 0.6745 / mad[:, None]

                if debug:

                    time_plot = 120*srate

                    plt.pcolormesh(mat[:,:time_plot])
                    plt.colorbar()
                    plt.show()

                    plt.plot(np.median(mat[:10], axis=0))
                    plt.show()

                    vmin, vmax = np.percentile(mat_norm[:,:time_plot].reshape(-1), 1), np.percentile(mat_norm[:,:time_plot].reshape(-1), 99)
                    plt.pcolormesh(mat_norm[:,:time_plot], vmin=vmin, vmax=vmax)
                    plt.colorbar()
                    plt.show()

                    plt.hist(mat_norm.ravel(), bins=200)
                    plt.show()

                tf_allconv_norm[chan_i] = mat_norm.astype(np.float32, copy=False)

            joblib.Parallel(n_jobs=n_core, prefer="threads")(joblib.delayed(norm_chan_conv)(chan_i) for chan_i in range(len(chanlist)))

        if norm_param == None:

            print(f'NO NORM', flush=True)

            tf_allconv_norm[:] = tf_allconv

        if debug:

            for chan_i in range(10):
                _tf = tf_allconv_norm[chan_i]
                min, max = np.percentile(_tf, 1), np.percentile(_tf, 99)
                plt.pcolormesh(_tf, vmin=min, vmax=max)
                plt.colorbar()
                plt.show()

            plt.hist(_tf.reshape(-1))
            plt.show()

            _tf.reshape(-1).max()
          
        #### stretch
        print('STRETCH')

        cond_name, trial_num = trial.split('_')[0], trial[-1]
        respfeature_trial = get_respfeatures_relabel(sujet).query(f"cond == '{cond_name}' and trial == {trial_num}")

        if debug:

            plt.scatter(respfeature_trial['inspi_index'].values, respfeature_trial['inspi_index'].values, label='inspi')
            plt.scatter(respfeature_trial['expi_index'].values, respfeature_trial['expi_index'].values, label='expi')
            plt.scatter(respfeature_trial['next_inspi_index'].values, respfeature_trial['next_inspi_index'].values, label='next_inspi')
            plt.legend()
            plt.show()

            cycle_times = respfeature_trial[['inspi_time', 'expi_time', 'next_inspi_time']].values
            (cycle_times[1:, 0] == cycle_times[:-1, -1]).all()

        tf_allchan_stretch = []

        for chan_i, _ in enumerate(chanlist):

            print_advancement(chan_i, len(chanlist), [25,50,75])

            tf_allchan_stretch.append(stretch_data_tf(respfeature_trial, stretch_point_TF, tf_allconv_norm[chan_i], srate)[0].astype(np.float32))

        tf_allchan_stretch = np.array(tf_allchan_stretch)

        resp_stretch = stretch_data(respfeature_trial, stretch_point_TF, resp, srate)[0]
        CO2_stretch = stretch_data(respfeature_trial, stretch_point_TF, CO2, srate)[0]

        if debug:
            for cycle_i in range(resp_stretch.shape[0]):
                plt.plot(resp_stretch[cycle_i])
            plt.show()

            plt.pcolormesh(np.median(tf_allchan_stretch[0], axis=[0]))
            plt.show()

        #### extract pre/post/condition
        if debug:
            plt.hist(respfeature_trial['cycle_duration'], bins=50)
            plt.show()

            plt.hist(respfeature_trial['inspi_duration'], bins=50)
            plt.show()
        
        #### save stretch tf
        print(f'SAVE TF STRETCH {sujet} {trial}', flush=True)
        dict_xr_tf_coords = {'chan' : chanlist, 'cycle' : np.arange(tf_allchan_stretch.shape[1]), 'freq' : frex, 'phase' : np.arange(stretch_point_TF)} 
        xr_tf_stretch_export = xr.DataArray(tf_allchan_stretch, dims=dict_xr_tf_coords.keys(), coords=dict_xr_tf_coords)
        xr_tf_stretch_export.to_netcdf(os.path.join(path_export_TF, f"{sujet}_{trial}_tf_stretch.nc"))

        #### remove
        os.remove(os.path.join(path_memmap, f'memmap_{sujet}_{trial}_tf_conv.npy'))
        os.remove(os.path.join(path_memmap, f'memmap_{sujet}_{trial}_rscore_param.npy'))
        os.remove(os.path.join(path_memmap, f'memmap_{sujet}_{trial}_tf_conv_norm.npy'))

        del tf_allconv, rscore_param, tf_allconv_norm









################################
######## EXTRACT POWER ########
################################


def extract_power(sujet):

    print("EXTRACT POWER", flush=True)

    #### verify if already computed
    path_export_Pxx = os.path.join(path_precompute, 'Pxx')

    if os.path.exists(os.path.join(path_export_Pxx, f'{sujet}_xr_Pxx.nc')):
        print(f'{sujet} ALREADY COMPUTED', flush=True)
        return

    #### params
    chanlist, _, _, localist = get_chanlist_localist(sujet)

    idx = np.arange(stretch_point_TF)
    inspi_sel = idx <= stretch_point_TF / 2
    expi_sel  = idx >  stretch_point_TF / 2

    bands = list(freq_band_dict.keys())
    phases = ["inspi", "expi"]

    band_frex_sel = {band: (frex >= freq[0]) & (frex <= freq[-1])
                    for band, freq in freq_band_dict.items()}
    
    #### load
    path_load_precompute_tf = os.path.join(path_precompute, 'TF', sujet)
    trial_list = get_alltrials_sujet(sujet)

    #### extract Pxx
    xr_Pxx = []

    for cond in conditions:

        trial_list_cond = [trial for trial in trial_list if trial.find(cond) != -1]

        Pxx_trial_list = []

        for trial in trial_list_cond:

            print(f"{trial}")

            #### load
            xr_tf = xr.open_dataarray(os.path.join(path_load_precompute_tf, f"{sujet}_{trial}_tf_stretch.nc"))

            #### extract
            Pxx_trial = np.zeros((1, xr_tf['chan'].shape[0], len(bands), len(phases), xr_tf['cycle'].shape[0]), dtype=np.float32)

            for chan_i, chan_name in enumerate(chanlist):

                ncycle = xr_tf['cycle'].shape[0]

                for cycle_i in range(ncycle):

                    for band_i, (band, frex_sel) in enumerate(band_frex_sel.items()):

                        _tf = xr_tf.sel(chan=chan_name, cycle=cycle_i, freq=frex[frex_sel]).values

                        pxx_inspi = np.median(_tf[:, inspi_sel])
                        pxx_expi  = np.median(_tf[:, expi_sel])

                        Pxx_trial[0, chan_i, band_i, 0, cycle_i] = pxx_inspi
                        Pxx_trial[0, chan_i, band_i, 1, cycle_i] = pxx_expi

            Pxx_trial_list.append(Pxx_trial)

        xr_Pxx.append(np.concat(Pxx_trial_list, axis=-1))

    #### construct xr
    max_cycle = np.array([_Pxx.shape[-1] for _Pxx in xr_Pxx]).max()
    xr_Pxx_data = np.full((1, len(conditions), xr_tf['chan'].shape[0], len(bands), len(phases), max_cycle), np.nan, dtype=np.float32)

    for cond_i, cond in enumerate(conditions):
        
        cond_max_cycle = xr_Pxx[cond_i].shape[-1]
        xr_Pxx_data[0,cond_i,:,:,:,:cond_max_cycle] = xr_Pxx[cond_i]

    da_Pxx = xr.DataArray(
        xr_Pxx_data,
        dims=("sujet", "cond", "chan", "band", "phase", "cycle"),
        coords={
            "sujet" : [sujet],
            "cond": conditions,
            "chan": chanlist,
            "band": bands,
            "phase": phases,
            "cycle": np.arange(max_cycle),
            "ROI": ("chan", localist['FS_volumetric_corrected'].values),   
        },
        name="Pxx"
    )


    #### SAVE
    print('SAVE EXTRACT PXX', flush=True)
    da_Pxx.to_netcdf(os.path.join(path_export_Pxx, f'{sujet}_xr_Pxx.nc')) 

    print('done', flush=True)






################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    #### PRECOMPUTE
    norm_param = 'rscore'
    
    #sujet = sujet_list[0]
    for sujet in sujet_list:

        precompute_tf(sujet, norm_param)
        extract_power(sujet)





                        