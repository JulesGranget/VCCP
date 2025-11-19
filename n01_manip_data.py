
from n00_config_params import *
from n00bis_config_analysis_functions import *
from n01bis_patient_info import *

import plotly.graph_objects as go
import seaborn as sns
import neo

debug = False



################################
######## LOAD DATA ########
################################

#y = _CO2
def find_best_lag(_respi, y, srate_sl):

    if debug:

        plt.plot(scipy.stats.zscore(_respi), label='respi')
        plt.plot(scipy.stats.zscore(y), label=f'CO2')
        plt.show()

    t = np.arange(y.size)
    p = np.polyfit(t, y, deg=10)
    trend = np.polyval(p, t)

    y_det = y - trend

    if debug:
        plt.plot(y, label="Original")
        plt.plot(trend, label="Fitted trend (deg=2)")
        plt.legend()
        plt.show()

        plt.plot(y, label="Original")
        plt.plot(y_det, label="detrended")
        plt.legend()
        plt.show()

    _respiz = scipy.stats.zscore(_respi)
    _peaks_respi, _ = scipy.signal.find_peaks(_respiz, height=_respiz.std(), distance=srate_sl*3)
    
    if debug:

        plt.plot(_respi)
        plt.scatter(_peaks_respi, _respi[_peaks_respi], color='r')
        plt.show()

    _peaks_y, _ = scipy.signal.find_peaks(y_det, height=y_det.std(), distance=srate_sl*3)

    if debug:

        plt.plot(y_det)
        plt.scatter(_peaks_y, y_det[_peaks_y], color='r')
        plt.legend()
        plt.show()

        imp_respi, imp_y = np.zeros(y.size), np.zeros(y.size)
        plt.scatter(_peaks_y, imp_respi[_peaks_y], label='y')
        plt.scatter(_peaks_respi, imp_y[_peaks_respi], label='respi')
        plt.legend()
        plt.show()
        
    min_peak_list = []

    for i, _py in enumerate(_peaks_y):
        _peak_diff = [_prespi - _py for _prespi in _peaks_respi] 
        min_peak_list.append(_peak_diff[np.abs(_peak_diff).argmin()])

    _shift_val_y = int(np.median(min_peak_list))

    if debug:

        _shifted_y = np.roll(y, _shift_val_y)
        plt.plot(scipy.stats.zscore(_respi), label='respi')
        plt.plot(scipy.stats.zscore(y), label=f'raw_y')
        plt.plot(scipy.stats.zscore(_shifted_y), label=f'shift_y')
        plt.legend()
        plt.show()

    return _shift_val_y



#sujet = 'NS217'
def load_physio_data(sujet):

    #### params
    # CO2_delay = 2.9 #sec
    # O2_delay = 1.59 #sec

    CO2_delay = 0 #sec
    O2_delay = 0 #sec

    #### Extraction data
    if sujet == 'jules':
        # df = pd.read_csv(os.path.join(path_data, 'hypotask_Jules_clean.csv'))
        # df = pd.read_csv(os.path.join(path_data, 'hypotask_Jules_clean_2.csv'))
        df = pd.read_csv(os.path.join(path_data, 'pilot', 'CPvolitional_Jules_3_clean.csv'))
    elif sujet == 'jose':
        df = pd.read_csv(os.path.join(path_data, 'pilot', 'jose_volitional_capnic_CP_protocol_clean.csv'))
    else:
        os.chdir(os.path.join(path_data, sujet))
        _file_sl = [file for file in os.listdir() if file.find('_cleaned.csv') != -1][0]
        df = pd.read_csv(os.path.join(path_data, sujet, _file_sl))
    
    time_raw = df['Session Time'].values[1:]
    flow_expi = df['Module 1: Raw Flow'].values[1:].astype(float)
    pressure = df['Module 2: Raw Pressure'].values[1:].astype(float)
    O2_raw = df['Internal Sensor: O2 Concentration'].values[1:].astype(float)
    flow_inspi = df['Module 4: Raw Flow'].values[1:].astype(float)
    trig = df['Digital I/O: Input 1'].values[1:].astype(float)
    CO2_raw = df['Serial Port 1: CO2 Concentration'].values[1:].astype(float)

    O2 = np.roll(O2_raw, O2_delay*srate_sl)
    CO2 = np.roll(CO2_raw, CO2_delay*srate_sl)

    if debug:

        select_time = 10*60*srate_sl
        x_scale = 2

        plt.plot(scipy.stats.zscore(flow_expi[:select_time]), label='flow_expi')
        plt.plot(scipy.stats.zscore(pressure[:select_time])+1*x_scale, label='pressure')
        plt.plot(scipy.stats.zscore(O2[:select_time])+2*x_scale, label='O2')
        plt.plot(scipy.stats.zscore(flow_inspi[:select_time])+3*x_scale, label='flow_inspi')
        plt.plot(scipy.stats.zscore(trig[:select_time])+4*x_scale, label='trig')
        plt.plot(scipy.stats.zscore(CO2[:select_time])+5*x_scale, label='CO2')
        plt.legend()
        plt.show()

        plt.plot(flow_expi[:select_time], label='flow_expi')
        plt.plot(flow_inspi[:select_time], label='flow_inspi')
        plt.legend()
        plt.show()

        plt.plot(O2[:select_time], label='O2')
        plt.plot(O2_raw[:select_time], label='O2_raw')
        plt.legend()
        plt.show()

        plt.plot(CO2[:select_time], label='CO2')
        plt.plot(CO2_raw[:select_time], label='CO2_raw')
        plt.legend()
        plt.show()

        plt.plot(scipy.stats.zscore(flow_expi), label='flow_expi')
        # plt.plot(scipy.stats.zscore(pressure), label='pressure')
        plt.plot(scipy.stats.zscore(CO2), label='CO2')
        plt.legend()
        plt.show()


    #### constructing respi

    respi = flow_inspi*-1 + flow_expi

    if debug:

        select_time = 10*60*srate_sl
        plt.plot(respi[:select_time])
        plt.show()

        plt.plot(scipy.stats.zscore(respi))
        plt.plot(scipy.stats.zscore(trig))
        plt.show()

        plt.plot(scipy.stats.zscore(respi))
        # plt.plot(scipy.stats.zscore(CO2))
        plt.plot(scipy.stats.zscore(CO2_raw))
        plt.show()

    #### extract respi/CO2

    if debug:

        plt.plot(scipy.stats.zscore(trig))
        plt.plot(scipy.stats.zscore(respi))
        plt.show()

        np.where(np.diff(trig) != 0)

    if sujet in ['jose', 'jules']:

        trig_dict = trig_dict_pilote[sujet]

    else:

        trig_dict = trig_dict_allpatient[sujet]

    if debug:
        
        plt.plot(trig)
        for cond in trig_dict:
            plt.vlines(trig_dict[cond], ymin=trig.min()+trig.std(), ymax=trig.max()-trig.std(), color='r')
        plt.show()

    data_dict_resp = {}
    data_dict_CO2_raw = {}
    data_dict_O2_raw = {}

    for cond in trig_dict:
        
        data_dict_resp[cond] = respi[trig_dict[cond][0]:trig_dict[cond][1]]
        data_dict_CO2_raw[cond] = CO2[trig_dict[cond][0]:trig_dict[cond][1]]
        data_dict_O2_raw[cond] = O2[trig_dict[cond][0]:trig_dict[cond][1]]

    return data_dict_resp, data_dict_CO2_raw, data_dict_O2_raw








def preproc_physio_sig(data_dict_resp, data_dict_CO2_raw, data_dict_O2_raw, srate_sl):
    
    if sujet in ['jose', 'jules']:

        trig_dict = trig_dict_pilote[sujet]

    else:

        trig_dict = trig_dict_allpatient[sujet]

    #### preproc respi

    data_dict_resp_clean = {}

    for cond in trig_dict:

        resp_clean = physio.preprocess(data_dict_resp[cond], srate_sl, band=25., btype='lowpass', ftype='bessel', order=5, normalize=False)
        resp_clean_smooth = physio.smooth_signal(resp_clean, srate_sl, win_shape='gaussian', sigma_ms=40.0)

        data_dict_resp_clean[cond] = resp_clean_smooth

    #### check instruction following and CO2 desync
    if debug:

        plt.plot(data_dict_resp['AoRoCo'], label='AoRoCo')
        plt.plot(data_dict_resp_clean['AoRoCo'], label='AoRoCo_clean')
        plt.legend()
        plt.show()

        # compare every blocks to baseline
        for cond in trig_dict:
            plt.plot(data_dict_resp_clean['AoRoCo01'], label='baseline')
            plt.plot(data_dict_resp_clean[cond], label='cond')
            plt.legend()
            plt.title(cond)
            plt.show()

        # check CO2 desync
        for cond in trig_dict:
            plt.plot(scipy.stats.zscore(data_dict_resp_clean[cond]), label='resp')
            plt.plot(scipy.stats.zscore(data_dict_CO2_raw[cond]), label='CO2')
            plt.legend()
            plt.title(cond)
            plt.show()

    #### correction respi timing to O2 / CO2
    data_dict_CO2 = {}
    data_dict_O2 = {}

    dict_shift_CO2 = {}
    dict_shift_O2 = {}

    #cond = 'A+R-Cc01'
    for cond in trig_dict:
        
        _respi = data_dict_resp_clean[cond]
        _CO2 = data_dict_CO2_raw[cond]
        _O2 = data_dict_O2_raw[cond]

        _CO2_lag = find_best_lag(_respi, _CO2, srate_sl)
        _O2_lag = find_best_lag(_respi, _O2, srate_sl)

        data_dict_CO2[cond] = np.roll(_CO2, _CO2_lag)
        data_dict_O2[cond] = np.roll(_O2, _O2_lag)

        dict_shift_CO2[cond] = _CO2_lag
        dict_shift_O2[cond] = _O2_lag

    #### extract respfeatures

    respfeatures = {}
    fig_detect_CO2 = {}
    fig_detect_O2 = {}

    for cond in trig_dict:

        _resp = data_dict_resp_clean[cond]
        _CO2 = data_dict_CO2[cond]
        _O2 = data_dict_O2[cond]
        
        cycles = physio.detect_respiration_cycles(_resp, srate_sl, method="crossing_baseline", baseline_mode='zero')

        resp_features_i = physio.compute_respiration_cycle_features(_resp, srate_sl, cycles, baseline=None)

        select_vec = np.ones((resp_features_i.index.shape[0]), dtype='int')
        resp_features_i.insert(resp_features_i.columns.shape[0], 'select', select_vec)

        time_vec = np.arange(_resp.shape[0])/srate_sl

        #### fig final
        fig_final_resp, ax = plt.subplots(figsize=(18, 10))
        ax.plot(time_vec, _resp)
        ax.scatter(cycles[:,0]/srate_sl, _resp[cycles[:,0]], color='g', label='inspi_selected')
        ax.scatter(cycles[:,1]/srate_sl, _resp[cycles[:,1]], color='c', label='expi_selected', marker='s')
        plt.legend()
        plt.title(cond)
        # plt.show()
        plt.close()

        fig_final_CO2, ax = plt.subplots(figsize=(18, 10))
        ax.plot(time_vec, _CO2)
        ax.scatter(cycles[:,0]/srate_sl, _CO2[cycles[:,0]], color='g', label='inspi_selected')
        ax.scatter(cycles[:,1]/srate_sl, _CO2[cycles[:,1]], color='c', label='expi_selected', marker='s')
        plt.legend()
        plt.title(cond)
        # plt.show()
        plt.close()

        fig_final_O2, ax = plt.subplots(figsize=(18, 10))
        ax.plot(time_vec, _O2)
        ax.scatter(cycles[:,0]/srate_sl, _O2[cycles[:,0]], color='g', label='inspi_selected')
        ax.scatter(cycles[:,1]/srate_sl, _O2[cycles[:,1]], color='c', label='expi_selected', marker='s')
        plt.legend()
        plt.title(cond)
        # plt.show()
        plt.close()
        
        respfeatures[cond] = [resp_features_i, fig_final_resp]
        fig_detect_CO2[cond] = [fig_final_CO2]
        fig_detect_O2[cond] = [fig_final_O2]

    #### stretch respi

    data_dict_resp_stretch = {}
    data_dict_CO2_stretch = {}
    data_dict_O2_stretch = {}

    nb_point_by_cycle = 500

    block_num = {'Ro' : 2, 'R-' : 3}

    for cond in conditions:

        _respfeature = respfeatures[f'{cond}01'][0]
        _resp = data_dict_resp_clean[f'{cond}01']
        _CO2 = data_dict_CO2[f'{cond}01']
        _O2 = data_dict_O2[f'{cond}01']

        _resp_stretch, _ = stretch_data(_respfeature, nb_point_by_cycle, _resp, srate_sl)
        _CO2_stretch, _ = stretch_data(_respfeature, nb_point_by_cycle, _CO2, srate_sl)
        _O2_stretch, _ = stretch_data(_respfeature, nb_point_by_cycle, _O2, srate_sl)
        
        if cond.find('Ro') != -1:
            
            for _block_i in range(block_num['Ro'])[1:]:

                if sujet == 'NS217' and cond == 'A+RoC-' and _block_i == 1:
                    continue

                _respfeature = respfeatures[f'{cond}0{_block_i+1}'][0]
                _resp = data_dict_resp_clean[f'{cond}0{_block_i+1}']
                _CO2 = data_dict_CO2[f'{cond}0{_block_i+1}']
                _O2 = data_dict_O2[f'{cond}0{_block_i+1}']

                _resp_stretch = np.concatenate([_resp_stretch, stretch_data(_respfeature, nb_point_by_cycle, _resp, srate_sl)[0]])
                _CO2_stretch = np.concatenate([_CO2_stretch, stretch_data(_respfeature, nb_point_by_cycle, _CO2, srate_sl)[0]])
                _O2_stretch = np.concatenate([_O2_stretch, stretch_data(_respfeature, nb_point_by_cycle, _O2, srate_sl)[0]])

        elif cond.find('R-') != -1:

            for _block_i in range(block_num['R-'])[1:]:

                if sujet == 'NS217' and cond == 'A+R-C-' and _block_i == 2:
                    continue

                _respfeature = respfeatures[f'{cond}0{_block_i+1}'][0]
                _resp = data_dict_resp_clean[f'{cond}0{_block_i+1}']
                _CO2 = data_dict_CO2[f'{cond}0{_block_i+1}']
                _O2 = data_dict_O2[f'{cond}0{_block_i+1}']

                _resp_stretch = np.concatenate([_resp_stretch, stretch_data(_respfeature, nb_point_by_cycle, _resp, srate_sl)[0]])
                _CO2_stretch = np.concatenate([_CO2_stretch, stretch_data(_respfeature, nb_point_by_cycle, _CO2, srate_sl)[0]])
                _O2_stretch = np.concatenate([_O2_stretch, stretch_data(_respfeature, nb_point_by_cycle, _O2, srate_sl)[0]])

        if debug:
            plt.plot(_resp_stretch.mean(axis=0))
            plt.show()

            plt.plot(_CO2_stretch.mean(axis=0))
            plt.show()

            plt.plot(_O2_stretch.mean(axis=0))
            plt.show()

        data_dict_resp_stretch[cond] = _resp_stretch
        data_dict_CO2_stretch[cond] = _CO2_stretch
        data_dict_O2_stretch[cond] = _O2_stretch
    
    #### pilote
    if sujet in ['jules', 'jose']:

        data_dict_resp_stretch = {}
        data_dict_CO2_stretch = {}
        data_dict_O2_stretch = {}

        nb_point_by_cycle = 500

        for cond in conditions:

            if cond == 'AoR-Co':

                _respfeature = respfeatures['AoR-Co01'][0]
                _resp = data_dict_resp_clean['AoR-Co01']
                _CO2 = data_dict_CO2['AoR-Co01']
                _O2 = data_dict_O2['AoR-Co01']

                _resp_stretch, _ = stretch_data(_respfeature, nb_point_by_cycle, _resp, srate_sl)
                _CO2_stretch, _ = stretch_data(_respfeature, nb_point_by_cycle, _CO2, srate_sl)
                _O2_stretch, _ = stretch_data(_respfeature, nb_point_by_cycle, _O2, srate_sl)

                _respfeature = respfeatures['AoR-Co02'][0]
                _resp = data_dict_resp_clean['AoR-Co02']
                _CO2 = data_dict_CO2['AoR-Co02']
                _O2 = data_dict_O2['AoR-Co02']

                _resp_stretch = np.concatenate([_resp_stretch, stretch_data(_respfeature, nb_point_by_cycle, _resp, srate_sl)[0]])
                _CO2_stretch = np.concatenate([_CO2_stretch, stretch_data(_respfeature, nb_point_by_cycle, _CO2, srate_sl)[0]])
                _O2_stretch = np.concatenate([_O2_stretch, stretch_data(_respfeature, nb_point_by_cycle, _O2, srate_sl)[0]])

            elif cond == 'A+R-Cc':

                _respfeature = respfeatures['A+R-Cc01'][0]
                _resp = data_dict_resp_clean['A+R-Cc01']
                _CO2 = data_dict_CO2['A+R-Cc01']
                _O2 = data_dict_O2['A+R-Cc01']

                _resp_stretch, _ = stretch_data(_respfeature, nb_point_by_cycle, _resp, srate_sl)
                _CO2_stretch, _ = stretch_data(_respfeature, nb_point_by_cycle, _CO2, srate_sl)
                _O2_stretch, _ = stretch_data(_respfeature, nb_point_by_cycle, _O2, srate_sl)

                _respfeature = respfeatures['A+R-Cc02'][0]
                _resp = data_dict_resp_clean['A+R-Cc02']
                _CO2 = data_dict_CO2['A+R-Cc02']
                _O2 = data_dict_O2['A+R-Cc02']

                _resp_stretch = np.concatenate([_resp_stretch, stretch_data(_respfeature, nb_point_by_cycle, _resp, srate_sl)[0]])
                _CO2_stretch = np.concatenate([_CO2_stretch, stretch_data(_respfeature, nb_point_by_cycle, _CO2, srate_sl)[0]])
                _O2_stretch = np.concatenate([_O2_stretch, stretch_data(_respfeature, nb_point_by_cycle, _O2, srate_sl)[0]])

            elif cond == 'A+R-C-':

                _respfeature = respfeatures['A+R-C-01'][0]
                _resp = data_dict_resp_clean['A+R-C-01']
                _CO2 = data_dict_CO2['A+R-C-01']
                _O2 = data_dict_O2['A+R-C-01']

                _resp_stretch, _ = stretch_data(_respfeature, nb_point_by_cycle, _resp, srate_sl)
                _CO2_stretch, _ = stretch_data(_respfeature, nb_point_by_cycle, _CO2, srate_sl)
                _O2_stretch, _ = stretch_data(_respfeature, nb_point_by_cycle, _O2, srate_sl)

                _respfeature = respfeatures['A+R-C-02'][0]
                _resp = data_dict_resp_clean['A+R-C-02']
                _CO2 = data_dict_CO2['A+R-C-02']
                _O2 = data_dict_O2['A+R-C-02']

                _resp_stretch = np.concatenate([_resp_stretch, stretch_data(_respfeature, nb_point_by_cycle, _resp, srate_sl)[0]])
                _CO2_stretch = np.concatenate([_CO2_stretch, stretch_data(_respfeature, nb_point_by_cycle, _CO2, srate_sl)[0]])
                _O2_stretch = np.concatenate([_O2_stretch, stretch_data(_respfeature, nb_point_by_cycle, _O2, srate_sl)[0]])

            else:
        
                _respfeature = respfeatures[cond][0]
                _resp = data_dict_resp_clean[cond]
                _CO2 = data_dict_CO2[cond]
                _O2 = data_dict_O2[cond]

                _resp_stretch, _ = stretch_data(_respfeature, nb_point_by_cycle, _resp, srate_sl)
                _CO2_stretch, _ = stretch_data(_respfeature, nb_point_by_cycle, _CO2, srate_sl)
                _O2_stretch, _ = stretch_data(_respfeature, nb_point_by_cycle, _O2, srate_sl)

            if debug:
                plt.plot(_resp_stretch.mean(axis=0))
                plt.show()

                plt.plot(_CO2_stretch.mean(axis=0))
                plt.show()

                plt.plot(_O2_stretch.mean(axis=0))
                plt.show()

        data_dict_resp_stretch[cond] = _resp_stretch
        data_dict_CO2_stretch[cond] = _CO2_stretch
        data_dict_O2_stretch[cond] = _O2_stretch

        if debug:

            for i in range(data_dict_CO2_stretch['A+RoCc'].shape[0]):
                plt.plot(data_dict_CO2_stretch['A+RoCc'][i,:])
            plt.show()

        

    #### extract CO2/O2
    endtidal_CO2 = pd.DataFrame()
    endtidal_O2 = pd.DataFrame()

    for cond in conditions:

        _CO2_stretch = data_dict_CO2_stretch[cond]
        _end_tidal_add = [_CO2_stretch[cycle_i][int(nb_point_by_cycle/2):].max() for cycle_i in range(_CO2_stretch.shape[0])]
        _cond_add = [cond] * _CO2_stretch.shape[0]
        endtidal_CO2 = pd.concat([endtidal_CO2, pd.DataFrame({'cond' : _cond_add, 'end_tidal' : _end_tidal_add})])

        _O2_stretch = data_dict_O2_stretch[cond]
        _end_tidal_add = [_O2_stretch[cycle_i][int(nb_point_by_cycle/2):].max() for cycle_i in range(_O2_stretch.shape[0])]
        _cond_add = [cond] * _O2_stretch.shape[0]
        endtidal_O2 = pd.concat([endtidal_O2, pd.DataFrame({'cond' : _cond_add, 'end_tidal' : _end_tidal_add})])
    
    #### export features

    os.chdir(os.path.join(path_results, 'respi', sujet, 'detection'))

    for cond in trig_dict:

        respfeatures[cond][0].to_excel(f"{cond}_respfeatures.xlsx")
        respfeatures[cond][1].savefig(f"{cond}_respfig0.jpeg")

        fig_detect_CO2[cond][0].savefig(f"{cond}_CO2fig0.jpeg")
        fig_detect_O2[cond][0].savefig(f"{cond}_O2fig0.jpeg")

    #### export mean fig

    os.chdir(os.path.join(path_results, 'respi', sujet, 'mean_figure'))

    fig_mean_resp, ax = plt.subplots()

    for cond in conditions:

        _data_stretch = data_dict_resp_stretch[cond]
        _mean = _data_stretch.mean(axis=0)
        _std = _data_stretch.std(axis=0)
        (line,) = ax.plot(_mean, label=f"{cond} n:{_data_stretch.shape[0]}")
        ax.fill_between(range(_mean.size), _mean - _std, _mean + _std, color=line.get_color(), alpha=0.1)

    plt.legend()
    plt.title('RESP')
    # plt.show()

    fig_mean_resp.savefig(f"MEAN_RESP.jpeg")

    fig_mean_CO2, ax = plt.subplots()

    for cond in conditions:

        _data_stretch = data_dict_CO2_stretch[cond]
        _mean = _data_stretch.mean(axis=0)
        _std = _data_stretch.std(axis=0)
        (line,) = ax.plot(_mean, label=f"{cond} n:{_data_stretch.shape[0]}")
        ax.fill_between(range(_mean.size), _mean - _std, _mean + _std, color=line.get_color(), alpha=0.1)

    plt.legend()
    plt.title('CO2')
    # plt.show()

    fig_mean_CO2.savefig(f"MEAN_CO2.jpeg")

    fig_mean_O2, ax = plt.subplots()

    for cond in conditions:

        _data_stretch = data_dict_O2_stretch[cond]
        _mean = _data_stretch.mean(axis=0)
        _std = _data_stretch.std(axis=0)
        (line,) = ax.plot(_mean, label=f"{cond} n:{_data_stretch.shape[0]}")
        ax.fill_between(range(_mean.size), _mean - _std, _mean + _std, color=line.get_color(), alpha=0.1)

    plt.legend()
    plt.title('O2')
    # plt.show()

    fig_mean_O2.savefig(f"MEAN_O2.jpeg")

    #### extract cond all respfeatures

    respfeatures_allcond = pd.DataFrame()

    block_num = {'Ro' : 2, 'R-' : 3}

    for block in trig_dict:

        _respfeature_add = respfeatures[block][0]
        _respfeature_add['cond'] = [block] * _respfeature_add.shape[0]
        respfeatures_allcond = pd.concat([respfeatures_allcond, _respfeature_add])

    for cond in conditions:

        if cond.find('Ro') != -1:

            for block_i in range(block_num['Ro']):

                respfeatures_allcond["cond"] = respfeatures_allcond["cond"].replace(f"{cond}0{block_i+1}", cond)

        elif cond.find('R-') != -1:

            for block_i in range(block_num['R-']):

                respfeatures_allcond["cond"] = respfeatures_allcond["cond"].replace(f"{cond}0{block_i+1}", cond)

    exclude_hyperventilation_vec = respfeatures_allcond['cond'] != 'A+RoCo00'
    respfeatures_allcond = respfeatures_allcond[exclude_hyperventilation_vec]

    #### reroll cycles & plot respfeatures
    # thresh_C = np.median(endtidal_CO2['end_tidal'].values)
    # thresh_A = np.median(respfeatures_allcond['total_amplitude'].values)
    # thresh_R = np.median(respfeatures_allcond['cycle_freq'].values)

    # sd_factor = 0.5
    # Co_upT = np.median(endtidal_CO2.query(f"cond == 'AoRoCo'")['end_tidal'].values) + sd_factor * scipy.stats.median_abs_deviation(endtidal_CO2.query(f"cond == 'AoRoCo'")['end_tidal'].values)
    # Co_dwT = np.median(endtidal_CO2.query(f"cond == 'AoRoCo'")['end_tidal'].values) - sd_factor * scipy.stats.median_abs_deviation(endtidal_CO2.query(f"cond == 'AoRoCo'")['end_tidal'].values)
    # Ao_upT = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values) + sd_factor * scipy.stats.median_abs_deviation(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values)
    # Ao_dwT = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values) - sd_factor * scipy.stats.median_abs_deviation(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values)
    # Ro_upT = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values) + sd_factor * scipy.stats.median_abs_deviation(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values)
    # Ro_dwT = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values) - sd_factor * scipy.stats.median_abs_deviation(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values)

    Co_upT = np.percentile(endtidal_CO2.query(f"cond == 'AoRoCo'")['end_tidal'].values, 100) 
    Co_dwT = np.percentile(endtidal_CO2.query(f"cond == 'AoRoCo'")['end_tidal'].values, 25)
    Ao_upT = np.percentile(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values, 85)
    Ao_dwT = np.percentile(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values, 0)
    Ro_upT = np.percentile(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values, 100)
    Ro_dwT = np.percentile(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values, 25)

    ARC_code_conditions = {'AoRoCo' : [1, 1, 1], 'AoR-Co' : [1, 0, 1], 'A+RoCc' : [0, 1, 1], 
                           'A+R-Cc' : [0, 0, 1], 'A+RoC-' : [0, 1, 0], 'A+R-C-' : [0, 0, 0]} #code ARC 1 means baseline

    df_resp_cycle_id = pd.DataFrame()

    for cycle_i in range(respfeatures_allcond.shape[0]):

        _A, _R, _C = respfeatures_allcond.iloc[cycle_i]['total_amplitude'], respfeatures_allcond.iloc[cycle_i]['cycle_freq'], endtidal_CO2.iloc[cycle_i]['end_tidal']
        
        if _A < Ao_upT and _A > Ao_dwT:
            _A_id = 1
        elif _A < Ao_dwT:
            _A_id = 'NA'
        elif _A > Ao_upT:
            _A_id = 0

        if _R < Ro_upT and _R > Ro_dwT:
            _R_id = 1
        elif _R < Ro_dwT:
            _R_id = 0
        elif _R > Ro_upT:
            _R_id = 'NA'

        if _C < Co_upT and _C > Co_dwT:
            _C_id = 1
        elif _C < Co_dwT:
            _C_id = 0
        elif _C > Co_upT:
            _C_id = 'NA'

        if 'NA' in [_A_id, _R_id, _C_id]:
            _cond_id = ['NA']
        else:
            _cond_id = [k for k, v in ARC_code_conditions.items() if v == [_A_id, _R_id, _C_id]]
        
        if len(_cond_id) > 1:
            raise ValueError('too much cond_id')
        if len(_cond_id) == 0:
            _cond_id = ['NA']

        _df = pd.DataFrame({'A' : [_A_id], 'R' : [_R_id], 'C' : [_C_id], 'cond' : [_cond_id[0]]})
        df_resp_cycle_id = pd.concat([df_resp_cycle_id, _df])

    respfeatures_allcond_balanced = respfeatures_allcond.copy()
    respfeatures_allcond_balanced['cond'] = df_resp_cycle_id['cond'].values

    endtidal_CO2_balanced = endtidal_CO2.copy()
    endtidal_CO2_balanced['cond'] = df_resp_cycle_id['cond'].values

    if debug:

        counts, _, _ = plt.hist(endtidal_CO2['end_tidal'], bins=50)
        plt.vlines([Co_dwT, Co_upT], ymin=0, ymax=counts.max(), color='r')
        plt.show()

        counts, _, _ = plt.hist(respfeatures_allcond['total_amplitude'], bins=50)
        plt.vlines([Ao_dwT, Ao_upT], ymin=0, ymax=counts.max(), color='r')
        plt.show()

        counts, _, _ = plt.hist(respfeatures_allcond['cycle_freq'], bins=50)
        plt.vlines([Ro_dwT, Ro_upT], ymin=0, ymax=counts.max(), color='r')
        plt.show()



        counts, _, _ = plt.hist(endtidal_CO2.query(f"cond == 'AoRoCo'")['end_tidal'], bins=50)
        plt.vlines([Co_dwT, Co_upT], ymin=0, ymax=counts.max(), color='r')
        plt.title('CO2')
        plt.show()

        counts, _, _ = plt.hist(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'], bins=50)
        plt.vlines([Ao_dwT, Ao_upT], ymin=0, ymax=counts.max(), color='r')
        plt.title('Amplitude')
        plt.show()

        counts, _, _ = plt.hist(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'], bins=50)
        plt.vlines([Ro_dwT, Ro_upT], ymin=0, ymax=counts.max(), color='r')
        plt.title('Freq')
        plt.show()


    #### export respfeatures
    respfeatures_allcond
    respfeatures_allcond_balanced
    respfeatures_allinfo = respfeatures_allcond.copy()
    df_resp_cycle_id = df_resp_cycle_id.rename(columns={'cond' : 'cond_bal'})
    respfeatures_allinfo = respfeatures_allinfo.rename(columns={'cond' : 'cond_raw'})
    respfeatures_allinfo = pd.concat([respfeatures_allinfo.reset_index(drop=True), df_resp_cycle_id.reset_index(drop=True)], axis=1)
    respfeatures_allinfo['etCO2'] = endtidal_CO2['end_tidal'].reset_index(drop=True)
    respfeatures_allinfo['etO2'] = endtidal_O2['end_tidal'].reset_index(drop=True)

    os.chdir(os.path.join(path_results, 'respi', 'respfeatures'))
    respfeatures_allinfo.to_excel(f"{sujet}_respfeatures.xlsx")

    #### baseline fig non balanced
    os.chdir(os.path.join(path_results, 'respi', sujet, 'cycles_identification'))

    color_cond = {
    'AoRoCo': 'red',
    'AoR-Co': 'blue',
    'A+RoCc': 'black',
    'A+R-Cc': 'cyan',
    'A+RoC-': 'green',
    'A+R-C-': 'yellow'
    }

    fig_baseline_non_balanced = go.Figure()

    cond = 'AoRoCo'
    _sel_vec = respfeatures_allcond['cond'] == cond
    _sel_CO2 = endtidal_CO2['end_tidal'].values[_sel_vec]
    _sel_amp = respfeatures_allcond['total_amplitude'].values[_sel_vec]
    _sel_freq = respfeatures_allcond['cycle_freq'].values[_sel_vec]

    fig_baseline_non_balanced.add_trace(go.Scatter3d(
        x=_sel_CO2,
        y=_sel_amp,
        z=_sel_freq,
        mode='markers',
        marker=dict(
            size=4,
            color=color_cond[cond]    # your color dictionary
        ),
        name=f"{cond}"               # label in legend
    ))

    _sel_vec = respfeatures_allcond['cond'] != cond
    _sel_CO2 = endtidal_CO2['end_tidal'].values[_sel_vec]
    _sel_amp = respfeatures_allcond['total_amplitude'].values[_sel_vec]
    _sel_freq = respfeatures_allcond['cycle_freq'].values[_sel_vec]

    fig_baseline_non_balanced.add_trace(go.Scatter3d(
        x=_sel_CO2,
        y=_sel_amp,
        z=_sel_freq,
        mode='markers',
        marker=dict(
            size=4,
            color='black'    # your color dictionary
        ),
        name=f"!= {cond}"               # label in legend
    ))

    fig_baseline_non_balanced.update_layout(
        scene=dict(
            xaxis_title="CO2",
            yaxis_title="AMP",
            zaxis_title="FREQ"
        ),
        legend=dict(title="Condition"),
        width=900,
        height=700,
        title={'text':"RAW BASELINE CYCLE", 'x':0.5, 'xanchor':"center", 'yanchor':"top", 'y':0.95}

    )

    # fig_baseline_non_balanced.show()
    fig_baseline_non_balanced.write_html(f"{sujet}_baseline_non_balanced.html")

    #### baseline fig balanced
    fig_baseline_balanced = go.Figure()

    cond = 'AoRoCo'
    _sel_vec = respfeatures_allcond_balanced['cond'] == cond
    _sel_CO2 = endtidal_CO2_balanced['end_tidal'].values[_sel_vec]
    _sel_amp = respfeatures_allcond_balanced['total_amplitude'].values[_sel_vec]
    _sel_freq = respfeatures_allcond_balanced['cycle_freq'].values[_sel_vec]

    fig_baseline_balanced.add_trace(go.Scatter3d(
        x=_sel_CO2,
        y=_sel_amp,
        z=_sel_freq,
        mode='markers',
        marker=dict(
            size=4,
            color=color_cond[cond]    # your color dictionary
        ),
        name=f"{cond}"               # label in legend
    ))

    _sel_vec = respfeatures_allcond_balanced['cond'] != cond
    _sel_CO2 = endtidal_CO2_balanced['end_tidal'].values[_sel_vec]
    _sel_amp = respfeatures_allcond_balanced['total_amplitude'].values[_sel_vec]
    _sel_freq = respfeatures_allcond_balanced['cycle_freq'].values[_sel_vec]

    fig_baseline_balanced.add_trace(go.Scatter3d(
        x=_sel_CO2,
        y=_sel_amp,
        z=_sel_freq,
        mode='markers',
        marker=dict(
            size=4,
            color='black'    # your color dictionary
        ),
        name=f"!= {cond}"               # label in legend
    ))

    



    fig_baseline_balanced.update_layout(
        scene=dict(
            xaxis_title="CO2",
            yaxis_title="AMP",
            zaxis_title="FREQ"
        ),
        legend=dict(title="Condition"),
        width=900,
        height=700,
        title={'text':"RAW BASELINE CYCLE BALANCED", 'x':0.5, 'xanchor':"center", 'yanchor':"top", 'y':0.95}

    )

    # fig_baseline_balanced.show()
    fig_baseline_balanced.write_html(f"{sujet}_baseline_balanced.html")

    #### raw cond sel
    fig_raw = go.Figure()

    for cond in conditions:
        _sel_vec = respfeatures_allcond['cond'] == cond
        _sel_CO2 = endtidal_CO2['end_tidal'].values[_sel_vec]
        _sel_amp = respfeatures_allcond['total_amplitude'].values[_sel_vec]
        _sel_freq = respfeatures_allcond['cycle_freq'].values[_sel_vec]

        fig_raw.add_trace(go.Scatter3d(
            x=_sel_CO2,
            y=_sel_amp,
            z=_sel_freq,
            mode='markers',
            marker=dict(
                size=4,
                color=color_cond[cond]    # your color dictionary
            ),
            name=f"{cond}"               # label in legend
        ))

    fig_raw.update_layout(
        scene=dict(
            xaxis_title="CO2",
            yaxis_title="AMP",
            zaxis_title="FREQ"
        ),
        legend=dict(title="Condition"),
        width=900,
        height=700,
        title={'text':"RAW", 'x':0.5, 'xanchor':"center", 'yanchor':"top", 'y':0.95}

    )

    # fig_raw.show()
    fig_raw.write_html(f"{sujet}_cond_raw.html")

    #### balanced cond sel
    fig_balanced = go.Figure()

    for cond in conditions:
        _sel_vec = respfeatures_allcond_balanced['cond'] == cond
        _sel_CO2 = endtidal_CO2_balanced['end_tidal'].values[_sel_vec]
        _sel_amp = respfeatures_allcond_balanced['total_amplitude'].values[_sel_vec]
        _sel_freq = respfeatures_allcond_balanced['cycle_freq'].values[_sel_vec]

        fig_balanced.add_trace(go.Scatter3d(
            x=_sel_CO2,
            y=_sel_amp,
            z=_sel_freq,
            mode='markers',
            marker=dict(
                size=4,
                color=color_cond[cond]    # your color dictionary
            ),
            name=f"{cond}"               # label in legend
        ))

    fig_balanced.update_layout(
        scene=dict(
            xaxis_title="CO2",
            yaxis_title="AMP",
            zaxis_title="FREQ"
        ),
        legend=dict(title="Condition"),
        width=900,
        height=700,
        title={'text':"BALANCED", 'x':0.5, 'xanchor':"center", 'yanchor':"top", 'y':0.95}

    )

    # fig_balanced.show()
    fig_balanced.write_html(f"{sujet}_cond_balanced.html")

    #### histplot
    os.chdir(os.path.join(path_results, 'respi', sujet, 'mean_figure'))

    list_to_keep = ["cycle_duration", "inspi_duration", "expi_duration", "cycle_freq", "cycle_ratio"]
    df_plot = respfeatures_allcond.melt(value_vars=list_to_keep, id_vars=[col for col in respfeatures_allcond.columns if col not in list_to_keep],
                  var_name="feature_name", 
                  value_name="value")

    g = sns.catplot(kind='box', data=df_plot, x='cond', y='value', col="feature_name", sharey=False, showfliers=False)
    plt.suptitle('RESP DURATION')
    plt.tight_layout()
    # plt.show()

    g.savefig("RESP_DURATION.png")


    list_to_keep = ["inspi_volume", "expi_volume", "total_volume", "inspi_amplitude", "expi_amplitude", "total_amplitude"]
    df_plot = respfeatures_allcond.melt(value_vars=list_to_keep, id_vars=[col for col in respfeatures_allcond.columns if col not in list_to_keep],
                  var_name="feature_name", 
                  value_name="value")

    g = sns.catplot(kind='box', data=df_plot, x='cond', y='value', col="feature_name", sharey=False, showfliers=False)
    plt.suptitle('RESP VOLUME')
    plt.tight_layout()
    # plt.show()

    g.savefig("RESP_VOLUME.png")



    df_plot = endtidal_CO2

    g = sns.catplot(kind='box', data=df_plot, x='cond', y='end_tidal', sharey=False, showfliers=False)
    plt.suptitle('endtidal CO2')
    plt.tight_layout()
    # plt.show()

    g.savefig("endtidalCO2.png")

    df_plot = endtidal_O2

    g = sns.catplot(kind='box', data=df_plot, x='cond', y='end_tidal', sharey=False, showfliers=False)
    plt.suptitle('endtidal O2')
    plt.tight_layout()
    # plt.show()

    g.savefig("endtidalO2.png")

    plt.close('all')

    #### summary figure
    fig_summary, axs = plt.subplots(ncols=3, figsize=(20,6))

    ax = axs[0]
    df_plot = endtidal_CO2
    sns.boxplot(data=df_plot, x='cond', y='end_tidal', ax=ax, showfliers=False)
    ax.set_title('endtital_CO2')

    ax = axs[1]
    sns.boxplot(data=respfeatures_allcond, x='cond', y='cycle_freq', ax=ax, showfliers=False)
    ax.set_title('cycle_freq')

    ax = axs[2]
    sns.boxplot(data=respfeatures_allcond, x='cond', y='total_amplitude', ax=ax, showfliers=False)
    ax.set_title('total_amplitude')

    plt.suptitle(sujet)
    # plt.show()

    fig_summary.savefig(f"SUMMARY_{sujet}.png")











def load_data(sujet):

    _path_data = os.path.join(path_data, sujet, datafile_name_allsujet[sujet]['B'])
    block = tdt.read_block(_path_data)

    _path_data = os.path.join(path_data, sujet, 'Baseline', 'FirstDeriv')
    reader = neo.io.TdtIO(dirname=_path_data)

    import tdt

    data = extract_tdt(reader)

    # Use read_block() to load the full dataset
    
def extract_tdt(reader):
    
    block = reader.read_block()

    data_allseg = {}
    srate_allseg = {}
    chanlist_allseg = {}
    trig_allseg = {}

    for seg_i, seg in enumerate(block.segments):

        print(f"EXTRACT SEGMENT : {seg_i}")

        chanlist_allseg[seg_i] = {}
        data_allseg[seg_i] = {}  # shape: (n_samples, n_channels)
        srate_allseg[seg_i] = {}

        for asig_i, asig in enumerate(seg.analogsignals):

            print(asig_i)

            aa = getattr(asig, "array_annotations", None)
            if aa and "channel_names" in aa:
                ch_names = [str(x) for x in aa["channel_names"]]
            else:
                if "channel_names" in asig.annotations:
                    ch_names = [str(x) for x in asig.annotations["channel_names"]]
                else:
                    arr = np.array(asig)
                    n_ch = arr.shape[1] if arr.ndim == 2 else 1
                    ch_names = [asig.name if (n_ch == 1 and asig.name) else f"ch{i}" for i in range(n_ch)]

            chanlist_allseg[seg_i][asig_i] = ch_names
            data_allseg[seg_i][asig_i] = np.array(asig)  # shape: (n_samples, n_channels)
            srate_allseg[seg_i][asig_i] = float(asig.sampling_rate.magnitude)

        trig_allseg[seg_i] = {}

        for ev_i, ev in enumerate(seg.events):

            print(f"EXTRACT event : {ev_i}")

            trig = ev.times.magnitude  # quantities array
            trig_allseg[seg_i][ev_i] = trig

    #### aggregates
    seg_num = len(data_allseg)
    if seg_num != 1:
        raise ValueError('WARNING SEG NUM != 1')
    
    sig_block_num = len(data_allseg[0])
    for block_i in range(len(data_allseg[0])):
        data_allseg[0][block_i].shape

    sig_block_num = len(data_allseg[0])
    for block_i in range(len(data_allseg[0])):
        len(chanlist_allseg[0][block_i])

    sig_block_num = len(data_allseg[0])
    for block_i in range(len(data_allseg[0])):
        chanlist_allseg[0][block_i]

    block_id_list = ['analog', 'Resp1', 'dRSP', 'sRSp', 'EEG', 'Fil', 'dsRP', 'EEG', 'rSte', 'fil']
    aux_data = data_allseg[0][block_id_list.index('analog')]
    srate_aux = srate_allseg[0][block_id_list.index('analog')]

    np.concat([data_allseg[0][block_id_list.index('analog')], data_allseg[0][block_id_list.index('analog')], data_allseg[0][block_id_list.index('analog')], data_allseg[0][block_id_list.index('analog')]])


    return data_allseg, srate_allseg, chanlist_allseg, trig_allseg


len(chanlist_allseg[0][10]).keys()



################################
######## LOAD DATA ########
################################

if __name__ == '__main__':

    srate_sl = 100

    sujet = 'jules'
    sujet = 'jose'
    sujet = 'NS217'

    data_dict_resp, data_dict_CO2_raw, data_dict_O2_raw = load_physio_data(sujet, srate_sl)


    preproc_physio_sig(data_dict_resp, data_dict_CO2_raw, data_dict_O2_raw, srate_sl)





  git config --global user.email "you@example.com"
  git config --global user.name "Your Name"
















