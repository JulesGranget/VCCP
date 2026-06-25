

import joblib
import plotly

from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *








################################
######## PRECOMPUTE TF ########
################################

#sujet, norm_param = sujet_list[0], 'rscore'
def extract_stretch_TF_raw(sujet, norm_param):

    #### config
    path_export_TF = os.path.join(path_precompute, 'TF', sujet)
    path_load_trial = os.path.join(path_prep, sujet, 'trial_exports')

    #### load
    trial_list = get_alltrials_sujet(sujet)
    respfeature = get_respfeatures_raw(sujet)

    #### precompute
    #cond = conditions[3]
    for cond in conditions:

        if os.path.exists(os.path.join(path_export_TF, f"{sujet}_{cond}_tf_stretch_cycle_cleaned_RAW.nc")):
            print(f'{sujet} {cond} ALREADY COMPUTED', flush=True)
            return
        else:
            print(f'#### TF PRECOMPUTE {sujet} {cond} ####', flush=True)

        trial_list_cond = [trial for trial in trial_list if trial.find(cond) != -1]
        tf_stretch_cond = []

        #trial = trial_list[-1]
        for trial in trial_list_cond:

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
            respfeature_trial_raw = respfeature.query(f"cond == '{cond_name}' and trial == {trial_num}")

            if debug:

                plt.scatter(respfeature_trial_raw['inspi_index'].values, respfeature_trial_raw['inspi_index'].values, label='inspi')
                plt.scatter(respfeature_trial_raw['expi_index'].values, respfeature_trial_raw['expi_index'].values, label='expi')
                plt.scatter(respfeature_trial_raw['next_inspi_index'].values, respfeature_trial_raw['next_inspi_index'].values, label='next_inspi')
                plt.legend()
                plt.show()

                cycle_times = respfeature_trial_raw[['inspi_time', 'expi_time', 'next_inspi_time']].values
                (cycle_times[1:, 0] == cycle_times[:-1, -1]).all()

            tf_allchan_stretch = []

            for chan_i, _ in enumerate(chanlist):

                print_advancement(chan_i, len(chanlist), [25,50,75])

                tf_allchan_stretch.append(stretch_data_tf(respfeature_trial_raw, stretch_point_TF, tf_allconv_norm[chan_i], srate)[0].astype(np.float32))

            tf_allchan_stretch = np.array(tf_allchan_stretch)

            resp_stretch = stretch_data(respfeature_trial_raw, stretch_point_TF, resp, srate)[0]
            CO2_stretch = stretch_data(respfeature_trial_raw, stretch_point_TF, CO2, srate)[0]

            if debug:
                for cycle_i in range(resp_stretch.shape[0]):
                    plt.plot(resp_stretch[cycle_i])
                plt.show()

                plt.pcolormesh(np.median(tf_allchan_stretch[0], axis=[0]))
                plt.show()

            #### sel good cycles
            sel_vec_good_cycles = respfeature_trial_raw.query(f"select == 1")['cycle'].values - respfeature_trial_raw['cycle'].values[0] 
            tf_allchan_stretch_filtered = tf_allchan_stretch[:,sel_vec_good_cycles]
            resp_stretch_filtered = resp_stretch[sel_vec_good_cycles]
            CO2_stretch_filtered = CO2_stretch[sel_vec_good_cycles]

            if debug:
                for cycle_i in range(resp_stretch_filtered.shape[0]):
                    plt.plot(resp_stretch_filtered[cycle_i])
                plt.show()

                for cycle_i in range(CO2_stretch_filtered.shape[0]):
                    plt.plot(CO2_stretch_filtered[cycle_i])
                plt.show()

                plt.pcolormesh(np.median(tf_allchan_stretch_filtered[0], axis=[0]))
                plt.show()

            tf_stretch_cond.append(tf_allchan_stretch_filtered)

        #### concat
        tf_stretch_cond = np.concat(tf_stretch_cond, axis=1)

        #### save stretch tf
        print(f'SAVE TF STRETCH {sujet} {cond}', flush=True)
        dict_xr_tf_coords = {'chan' : chanlist, 'cycle' : np.arange(tf_stretch_cond.shape[1]), 'freq' : frex, 'phase' : np.arange(stretch_point_TF)} 
        xr_tf_stretch_export = xr.DataArray(tf_stretch_cond, dims=dict_xr_tf_coords.keys(), coords=dict_xr_tf_coords)
        xr_tf_stretch_export.to_netcdf(os.path.join(path_export_TF, f"{sujet}_{cond}_tf_stretch_cycle_cleaned_RAW.nc"))

        #### remove
        os.remove(os.path.join(path_memmap, f'memmap_{sujet}_{trial}_tf_conv.npy'))
        os.remove(os.path.join(path_memmap, f'memmap_{sujet}_{trial}_rscore_param.npy'))
        os.remove(os.path.join(path_memmap, f'memmap_{sujet}_{trial}_tf_conv_norm.npy'))

        del tf_allconv, rscore_param, tf_allconv_norm, 




def extract_stretch_TF_relabel(sujet):

    #### config
    path_load_TF = os.path.join(path_precompute, 'TF', sujet)
    chanlist, _, _, localist = get_chanlist_localist(sujet)

    #### load
    respfeatures_relabel_raw = get_respfeatures_relabel_raw(sujet)
    respfeatures_relabel_raw = respfeatures_relabel_raw.query(f"select == 1")
    respfeatures_relabel_raw = respfeatures_relabel_raw.sort_values(['cond', 'cycle']).reset_index(drop=True)

    cycle_i_vec = []
    cycle_count = 0

    for row_i, row_val in respfeatures_relabel_raw.iterrows():

        if row_i != 0:
            cond_row = row_val['cond']
            if cond_row == cond_prev:
                cycle_count += 1
            else:
                cycle_count = 0
        cycle_i_vec.append(cycle_count)
        cond_prev = row_val['cond']

    respfeatures_relabel_raw['cycle'] = cycle_i_vec

    #### extract cond VCCP
    for cond in conditions_VCCP:

        print(f"RELABEL {cond}")

        _rf = respfeatures_relabel_raw.query(f"cond_relabel == '{cond}'")
        _sel_params = _rf[['cond', 'cycle']]

        cond_relabel_cycles = []

        for cond_target in _sel_params['cond'].unique():
            
            _tf_stretch = xr.load_dataarray(os.path.join(path_load_TF, f"{sujet}_{cond_target}_tf_stretch_cycle_cleaned.nc"))
            _sel_cycle_cond_relabel = _sel_params.query(f"cond == '{cond_target}'")['cycle'].values

            __tf_stretch_cond_relabel = _tf_stretch[:,_sel_cycle_cond_relabel].values
            cond_relabel_cycles.append(__tf_stretch_cond_relabel)

        tf_stretch_cond_relabel = np.concat(cond_relabel_cycles, axis=1)

        #### export
        dict_xr_tf_coords = {'chan' : chanlist, 'cycle' : np.arange(tf_stretch_cond_relabel.shape[1]), 'freq' : frex, 'phase' : np.arange(stretch_point_TF)} 
        xr_tf_stretch_export = xr.DataArray(tf_stretch_cond_relabel, dims=dict_xr_tf_coords.keys(), coords=dict_xr_tf_coords)
        xr_tf_stretch_export.to_netcdf(os.path.join(path_load_TF, f"{sujet}_{cond}_tf_stretch_cycle_cleaned_RELABEL.nc"))





################################
######## EXTRACT POWER ########
################################


def extract_power_RAW(sujet):

    print("EXTRACT POWER", flush=True)

    #### verify if already computed
    path_export_Pxx = os.path.join(path_precompute, 'Pxx')

    if os.path.exists(os.path.join(path_export_Pxx, f'{sujet}_xr_Pxx_RAW.nc')):
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

    #### extract Pxx
    xr_Pxx = []

    for cond in conditions:

        print(f"{cond}")

        #### load
        xr_tf = xr.open_dataarray(os.path.join(path_load_precompute_tf, f"{sujet}_{cond}_tf_stretch_cycle_cleaned_RAW.nc"))

        #### extract
        Pxx_cond = np.zeros((1, xr_tf['chan'].shape[0], len(bands), len(phases), xr_tf['cycle'].shape[0]), dtype=np.float32)

        for chan_i, chan_name in enumerate(chanlist):

            ncycle = xr_tf['cycle'].shape[0]

            for cycle_i in range(ncycle):

                for band_i, (band, frex_sel) in enumerate(band_frex_sel.items()):

                    _tf = xr_tf.sel(chan=chan_name, cycle=cycle_i, freq=frex[frex_sel]).values

                    pxx_inspi = np.median(_tf[:, inspi_sel])
                    pxx_expi  = np.median(_tf[:, expi_sel])

                    Pxx_cond[0, chan_i, band_i, 0, cycle_i] = pxx_inspi
                    Pxx_cond[0, chan_i, band_i, 1, cycle_i] = pxx_expi

        xr_Pxx.append(Pxx_cond)

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
            "ROI": ("chan", localist['loca'].values),   
        },
        name="Pxx"
    )

    #### SAVE
    print('SAVE EXTRACT PXX', flush=True)
    da_Pxx.to_netcdf(os.path.join(path_export_Pxx, f'{sujet}_xr_Pxx_RAW.nc')) 

    print('done', flush=True)




def extract_power_RELABEL(sujet):

    print("EXTRACT POWER", flush=True)

    #### verify if already computed
    path_export_Pxx = os.path.join(path_precompute, 'Pxx')

    if os.path.exists(os.path.join(path_export_Pxx, f'{sujet}_xr_Pxx_RELABEL.nc')):
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

    #### extract Pxx
    xr_Pxx = []

    for cond in conditions_VCCP:

        print(f"{cond}")

        #### load
        xr_tf = xr.open_dataarray(os.path.join(path_load_precompute_tf, f"{sujet}_{cond}_tf_stretch_cycle_cleaned_RELABEL.nc"))

        #### extract
        Pxx_cond = np.zeros((1, xr_tf['chan'].shape[0], len(bands), len(phases), xr_tf['cycle'].shape[0]), dtype=np.float32)

        for chan_i, chan_name in enumerate(chanlist):

            ncycle = xr_tf['cycle'].shape[0]

            for cycle_i in range(ncycle):

                for band_i, (band, frex_sel) in enumerate(band_frex_sel.items()):

                    _tf = xr_tf.sel(chan=chan_name, cycle=cycle_i, freq=frex[frex_sel]).values

                    pxx_inspi = np.median(_tf[:, inspi_sel])
                    pxx_expi  = np.median(_tf[:, expi_sel])

                    Pxx_cond[0, chan_i, band_i, 0, cycle_i] = pxx_inspi
                    Pxx_cond[0, chan_i, band_i, 1, cycle_i] = pxx_expi

        xr_Pxx.append(Pxx_cond)

    #### construct xr
    max_cycle = np.array([_Pxx.shape[-1] for _Pxx in xr_Pxx]).max()
    xr_Pxx_data = np.full((1, len(conditions_VCCP), xr_tf['chan'].shape[0], len(bands), len(phases), max_cycle), np.nan, dtype=np.float32)

    for cond_i, cond in enumerate(conditions_VCCP):
        
        cond_max_cycle = xr_Pxx[cond_i].shape[-1]
        xr_Pxx_data[0,cond_i,:,:,:,:cond_max_cycle] = xr_Pxx[cond_i]

    da_Pxx = xr.DataArray(
        xr_Pxx_data,
        dims=("sujet", "cond", "chan", "band", "phase", "cycle"),
        coords={
            "sujet" : [sujet],
            "cond": conditions_VCCP,
            "chan": chanlist,
            "band": bands,
            "phase": phases,
            "cycle": np.arange(max_cycle),
            "ROI": ("chan", localist['loca'].values),   
        },
        name="Pxx"
    )

    #### SAVE
    print('SAVE EXTRACT PXX', flush=True)
    da_Pxx.to_netcdf(os.path.join(path_export_Pxx, f'{sujet}_xr_Pxx_RELABEL.nc')) 

    print('done', flush=True)







################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    #### PRECOMPUTE
    norm_param = 'rscore'
    
    #sujet = sujet_list[0]
    for sujet in sujet_list:

        extract_stretch_TF_raw(sujet, norm_param)
        extract_stretch_TF_relabel(sujet)
        extract_power_RAW(sujet)
        extract_power_RELABEL(sujet)





                        