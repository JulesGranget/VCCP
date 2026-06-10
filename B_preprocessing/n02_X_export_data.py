
from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *

import plotly
import plotly.graph_objects as go
import seaborn as sns
import neo
import h5py


debug = False








################################
######## LOAD DATA ########
################################

#y = _CO2
def find_best_lag(_respi, y):

    if debug:

        time_vec = np.arange(_respi.shape[0]) / srate

        plt.plot(time_vec, scipy.stats.zscore(_respi), label='respi')
        plt.plot(time_vec, scipy.stats.zscore(y), label=f'CO2')
        plt.legend()
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
    _peaks_respi, _ = scipy.signal.find_peaks(_respiz, height=_respiz.std(), distance=srate*3)
    
    if debug:

        plt.plot(_respi)
        plt.scatter(_peaks_respi, _respi[_peaks_respi], color='r')
        plt.show()

    _peaks_y, _ = scipy.signal.find_peaks(y_det, height=y_det.std(), distance=srate*3)

    if debug:

        time_vec = np.arange(_respi.shape[0]) / srate

        plt.plot(y_det)
        plt.scatter(_peaks_y, y_det[_peaks_y], color='r')
        plt.legend()
        plt.show()

        imp_respi, imp_y = np.zeros(y.size), np.zeros(y.size)
        plt.scatter(time_vec[_peaks_y], imp_respi[_peaks_y], label='y')
        plt.scatter(time_vec[_peaks_respi], imp_y[_peaks_respi], label='respi')
        plt.legend()
        plt.show()
        
    min_peak_list = []

    for i, _py in enumerate(_peaks_y):
        _peak_diff = [_prespi - _py for _prespi in _peaks_respi] 
        min_peak_list.append(_peak_diff[np.abs(_peak_diff).argmin()])

    _shift_val_y = int(np.median(min_peak_list))

    if debug:

        time_vec = np.arange(_respi.shape[0]) / srate

        _shifted_y = np.roll(y, _shift_val_y)
        plt.plot(time_vec, scipy.stats.zscore(_respi), label='respi')
        plt.plot(time_vec, scipy.stats.zscore(y), label=f'raw_y')
        plt.plot(time_vec, scipy.stats.zscore(_shifted_y), label=f'shift_y')
        plt.legend()
        plt.show()

    return _shift_val_y



#sujet = 'NS217'
def export_data_sujet_pilote(sujet):    

    #### Extraction data
    os.chdir(os.path.join(path_prep))

    trial_info = pd.read_excel(f"{sujet}_trialsinfo.xlsx")
    trial_list = [_trial.replace('0', 'o') for _trial in trial_info['Type'].values]

    anatomy = pd.read_excel(f"{sujet}_anatomy.xlsx")
    
    f = h5py.File(f"{sujet}_VCCP_B30_P1_ecog_epoch_smartlab.mat", "r")    
    # print(list(f.keys()))
    # f["ecog"].keys()
    # f["ecog"]['epochs'].keys()
    # f["ecog"]['epochs']['chan_info'][:]

    n_trial = f["ecog"]["epochs"]['eeg'][:].shape[-1]

    data = []
    respi = []
    inspi_starts = []
    expi_starts = []
    gsm_signal = []
    co2 = []
    trial_list_i = []

    for _trial_i, _trial in enumerate(trial_list):

        ref = f["ecog"]["epochs"]['eeg'][0][_trial_i]
        data.append(np.transpose(f[ref][()], (1,0)))

        ref = f["ecog"]["epochs"]['inhonsets'][0][_trial_i]
        inspi_starts.append(f[ref][()].reshape(-1))

        ref = f["ecog"]["epochs"]['inhoffsets'][0][_trial_i]
        expi_starts.append(f[ref][()].reshape(-1))

        ref = f["ecog"]["epochs"]['iflow'][()][0][_trial_i]
        iflow = f[ref][()].reshape(-1)

        ref = f["ecog"]["epochs"]['eflow'][()][0][_trial_i]
        eflow = f[ref][()].reshape(-1)

        # ref = f["ecog"]["epochs"]['pressure'][()][0][_trial_i]
        # pressure = f[ref][()].reshape(-1)

        ref = f["ecog"]["epochs"]['gsm_mixer'][()][0][_trial_i]
        gsm_signal.append(f[ref][()].reshape(-1))

        ref = f["ecog"]["epochs"]['co2'][()][0][_trial_i]
        co2.append(f[ref][()].reshape(-1))

        _respi = iflow * -1 + eflow
        
        respi.append(_respi)

        if _trial_i == 0:
            trial_list_i.append([0, _respi.size])
        else:
            trial_list_i.append([trial_list_i[_trial_i-1][-1], _respi.size + trial_list_i[_trial_i-1][-1]])

        if debug:

            plt.plot(_respi)
            plt.show()

            plt.plot(iflow)
            plt.show()

            plt.plot(eflow)
            plt.show()

            plt.plot(gsm_signal)
            plt.show()

    #### inspect signal
    if debug:

        for trial_i in range(n_trial):

            data[trial_i].shape

        _sig = np.concat(data, axis=1)[0,:]
        plt.plot(_sig)
        plt.vlines(trial_list_i, ymin=_sig.min(), ymax=_sig.max(), colors='r')
        plt.show()

        _sig = np.concat(respi)
        plt.plot(_sig)
        plt.vlines([_cond_i[0] for _cond_i in trial_list_i], ymin=_sig.min(), ymax=_sig.max(), colors='r')
        plt.show()

        _sig = np.concat(gsm_signal)
        plt.plot(_sig)
        plt.vlines(trial_list_i, ymin=_sig.min(), ymax=_sig.max(), colors='r')
        plt.show()

        _sig = np.concat(co2)
        plt.plot(_sig)
        plt.vlines(trial_list_i, ymin=_sig.min(), ymax=_sig.max(), colors='r')
        plt.show()

    #### inspect to fill patient info
    trial_count = {_cond : np.where(np.array(trial_list) == _cond)[0].size for _cond in conditions}

        #### preproc respi
    respi_clean = []

    for _trial_i, _respi in enumerate(respi):

        resp_clean = physio.preprocess(_respi, srate, band=25., btype='lowpass', ftype='bessel', order=5, normalize=False)
        resp_clean_smooth = physio.smooth_signal(resp_clean, srate, win_shape='gaussian', sigma_ms=40.0)

        respi_clean.append(resp_clean_smooth)

        #### check instruction following and CO2 desync
        if debug:

            plt.plot(respi_clean['AoRoCo'], label='AoRoCo_clean')
            plt.legend()
            plt.show()

            # compare every blocks to baseline
            plt.plot(respi_clean['AoRoCo01'], label='baseline')
            plt.plot(respi_clean[cond], label='cond')
            plt.legend()
            plt.title(cond)
            plt.show()

            # check CO2 desync
            plt.plot(scipy.stats.zscore(respi_clean[cond]), label='resp')
            plt.legend()
            plt.title(cond)
            plt.show()

        #### correction respi timing to CO2
    co2_corr = []
    for _trial_i, _trial in enumerate(trial_list):
        
        _respi = respi_clean[_trial_i]
        _CO2 = co2[_trial_i]

        _CO2_lag = find_best_lag(_respi, _CO2)

        _CO2_corr = np.roll(_CO2, _CO2_lag)

        co2_corr.append(_CO2_corr)

    if debug:

        _trial_i = 0
        plt.plot(scipy.stats.zscore(respi_clean[_trial_i]), label='respi')
        plt.plot(scipy.stats.zscore(co2[_trial_i]), label='co2_raw')
        plt.plot(scipy.stats.zscore(co2_corr[_trial_i]), label='co2_corr')
        plt.legend()
        plt.show()

    #### export
    os.chdir(os.path.join(path_prep, f"{sujet}_exports"))

    for cond in conditions:

        if cond in trial_count:

            _trial_i_list = np.where(np.array(trial_list) == cond)[0]

            for _trial_i_count, _trial_i in enumerate(_trial_i_list):

                _aux_export = np.vstack([respi_clean[_trial_i], co2_corr[_trial_i], gsm_signal[_trial_i]])
                np.save(f"{sujet}_{cond}0{_trial_i_count+1}_aux.npy", _aux_export)

                _data_export = data[_trial_i]
                np.save(f"{sujet}_{cond}0{_trial_i_count+1}_data.npy", _data_export)
                
    
#sujet = sujet_list[0]
def get_data_raw(sujet):

    #### params
    path_filetrial = os.path.join(path_data, sujet, f"{sujet}_Experiment_Log.xlsx")
    trial_info = pd.read_excel(path_filetrial, sheet_name='trialsinfo')

    path_filepreproc = os.path.join(path_prep, sujet, f"{sujet}_VCCL_ecog_epoch_smartlab.mat")
    f = h5py.File(path_filepreproc, "r")    
    # print(list(f.keys()))
    # f["ecog"].keys()
    # f["ecog"]["ftrip"].keys()
    # f["ecog"]['epochs'].keys()
    # f["ecog"]['epochs']['eeglabels']
    # f["ecog"]['epochs']['eeglabels_FSurf'][:]

    n_trial = f["ecog"]["epochs"]['eeg'][:].shape[-1]

    if n_trial != trial_info.shape[0]:
        raise ValueError(f"!!! TRIAL NUMBER IN DATA AND PROTOCOL ARE DIFFERENT !!!")
    
    #### chan / loca list
    chanlist = []
    refs = f["ecog"]["epochs"]["eeglabels"]
    for i in range(refs.shape[0]):
        r = refs[i, 0]
        r_data = f[r][()]

        # flatten and convert numeric codes to characters
        s = "".join(chr(x) for x in r_data.flatten())
        chanlist.append(s)

    localist = []
    refs = f["ecog"]["epochs"]["eeglabels_FSurf"]
    for i, r in enumerate(refs):
        # r = refs[i, 0]
        # r_data = f[r][()]
        # # flatten and convert numeric codes to characters
        # s = "".join(chr(x) for x in r_data.flatten())
        localist.append(r)

    #### data
    data = []
    respi = []
    co2 = []
    trial_list_i = []

    for _trial_i in range(n_trial):

        ref = f["ecog"]["epochs"]['eeg'][0][_trial_i]
        data.append(np.transpose(f[ref][()], (1,0)))

        ref = f["ecog"]["epochs"]['iflow'][()][0][_trial_i]
        iflow = f[ref][()].reshape(-1)

        ref = f["ecog"]["epochs"]['eflow'][()][0][_trial_i]
        eflow = f[ref][()].reshape(-1)

        # ref = f["ecog"]["epochs"]['pressure'][()][0][_trial_i]
        # pressure = f[ref][()].reshape(-1)

        ref = f["ecog"]["epochs"]['co2'][()][0][_trial_i]
        co2.append(f[ref][()].reshape(-1))

        _respi = iflow * -1 + eflow
        respi.append(_respi)

        if _trial_i == 0:
            trial_list_i.append([0, _respi.size])
        else:
            trial_list_i.append([trial_list_i[_trial_i-1][-1], _respi.size + trial_list_i[_trial_i-1][-1]])

        if debug:

            plt.plot(_respi)
            plt.show()

            plt.plot(iflow)
            plt.plot(eflow)
            plt.show()

    #### verify
    if np.concat(data, axis=1).shape[0] != len(chanlist) != len(localist):
        raise ValueError(f"!!! DATA AND CHANLIST ARE NOT THE SAME LENGTH !!!")

    #### inspect signal
    if debug:

        for trial_i in range(n_trial):

            data[trial_i].shape

        _sig = np.concat(data, axis=1)[0,:]
        plt.plot(_sig)
        plt.vlines(trial_list_i, ymin=_sig.min(), ymax=_sig.max(), colors='r')
        plt.show()

        _sig = np.concat(respi)
        plt.plot(_sig)
        plt.vlines([_cond_i[0] for _cond_i in trial_list_i], ymin=_sig.min(), ymax=_sig.max(), colors='r')
        plt.show()

        _sig = np.concat(co2)
        plt.plot(_sig)
        plt.vlines(trial_list_i, ymin=_sig.min(), ymax=_sig.max(), colors='r')
        plt.show()

    #### preproc respi
    respi_clean = []

    for _trial_i, _respi in enumerate(respi):

        resp_clean = physio.preprocess(_respi, srate, band=25., btype='lowpass', ftype='bessel', order=5, normalize=False)
        resp_clean_smooth = physio.smooth_signal(resp_clean, srate, win_shape='gaussian', sigma_ms=40.0)

        respi_clean.append(resp_clean_smooth)

    #### correction respi timing to CO2
    co2_corr = []
    for _trial_i in range(n_trial):
        
        _respi = respi_clean[_trial_i]
        _CO2 = co2[_trial_i]

        _CO2_lag = find_best_lag(_respi, _CO2)

        _CO2_corr = np.roll(_CO2, _CO2_lag)

        co2_corr.append(_CO2_corr)

    if debug:

        _trial_i = 0
        plt.plot(scipy.stats.zscore(respi_clean[_trial_i]), label='respi')
        plt.plot(scipy.stats.zscore(co2[_trial_i]), label='co2_raw')
        plt.plot(scipy.stats.zscore(co2_corr[_trial_i]), label='co2_corr')
        plt.legend()
        plt.show()

    return data, respi_clean, co2_corr, chanlist, localist







########################################
######## GENERATE RAW LOCA FILE ########
########################################



def generate_raw_loca_file(sujet):

    #### get data
    data, resp_clean, co2_corr, chanlist, localist = get_data_raw(sujet)

    #### get loca info
    path_filechanloca = os.path.join(path_anatomy, f"{sujet}_Electrodes_Natus_TDT_correspondence.xlsx")
    df_loca_Nat_TDT_corresp = pd.read_excel(path_filechanloca)

    #### compare
    chan_in_data_and_sheet = np.array([1 if row_val['Label'] in chanlist else 0 for row_i, row_val in df_loca_Nat_TDT_corresp.iterrows()])
    chan_in_data_and_not_in_sheet = np.array([0 if chan in df_loca_Nat_TDT_corresp['Label'].values else 1 for chan in chanlist])

    chan_not_in_data = df_loca_Nat_TDT_corresp[~chan_in_data_and_sheet.astype('bool')]['Label'].values # remove from sheet
    chan_not_in_sheet = np.array(chanlist)[chan_in_data_and_not_in_sheet.astype('bool')] # remove from data

    chanlist_mask_final = ~chan_in_data_and_not_in_sheet.astype('bool')
    chanlist_final = np.array(chanlist)[chanlist_mask_final]
    localist_final = np.array(localist)[chanlist_mask_final]

    #### verify
    if (chan_in_data_and_sheet.sum() + chan_not_in_sheet.shape[0]) != data[0].shape[0]:
        raise ValueError('!!! CHAN NUMBER DOES NOT CORRESPOND WITH SHEET !!!')

    #### write
    path_export_chanlist = os.path.join(path_anatomy, f"{sujet}_chanlist_info", f"{sujet}_chan_not_in_data.txt")
    with open(path_export_chanlist, "w") as f:
        for _chan in chan_not_in_data:
            f.write(str(_chan) + "\n")

    path_export_chanlist = os.path.join(path_anatomy, f"{sujet}_chanlist_info", f"{sujet}_chan_not_in_sheet.txt")
    with open(path_export_chanlist, "w") as f:
        for _chan in chan_not_in_sheet:
            f.write(str(_chan) + "\n")

    path_export_chanlist = os.path.join(path_anatomy, f"{sujet}_chanlist_info", f"{sujet}_chan_selected.txt")
    with open(path_export_chanlist, "w") as f:
        for _chan in chanlist_final.tolist():
            f.write(str(_chan) + "\n")

    #### construct df_loca
    df_loca_Nat_TDT_corresp_filtered = df_loca_Nat_TDT_corresp[chan_in_data_and_sheet.astype('bool')][['Label', 'FS_vol', 'Desikan_Killiany', 'SOZ (1, 0/empty)', 'Spikey (1, 0/empty)', 'Out (1, 0/empty)','BAD (1, 0/empty)',
                                                                                                        'LEPTO_coords_1', 'LEPTO_coords_2', 'LEPTO_coords_3',
                                                                                                        'fsaverage_coords_1', 'fsaverage_coords_2', 'fsaverage_coords_3']]
    df_loca_Nat_TDT_corresp_filtered = df_loca_Nat_TDT_corresp_filtered.rename(columns={'FS_vol' : 'FD_vol_raw', 'Desikan_Killiany' : 'Desikan_Killiany_raw',
                                                                                        'SOZ (1, 0/empty)' : 'SOZ_raw', 'Spikey (1, 0/empty)' : 'Spikey_raw', 'Out (1, 0/empty)' : 'Out_raw',
                                                                                        'BAD (1, 0/empty)' : 'BAD_raw'})
    df_loca_Nat_TDT_corresp_filtered['Label_copy'] = df_loca_Nat_TDT_corresp_filtered['Label']
    df_loca_Nat_TDT_corresp_filtered['Spikey_jules'] = np.zeros((df_loca_Nat_TDT_corresp_filtered.shape[0]))
    df_loca_Nat_TDT_corresp_filtered['sig_inspected'] = np.zeros((df_loca_Nat_TDT_corresp_filtered.shape[0]))
    df_loca_Nat_TDT_corresp_filtered['Desikan_Killiany_jules'] = df_loca_Nat_TDT_corresp_filtered['Desikan_Killiany_raw']
    df_loca_Nat_TDT_corresp_filtered['loca_inspected'] = np.zeros((df_loca_Nat_TDT_corresp_filtered.shape[0]))
    df_loca_Nat_TDT_corresp_filtered['select'] = np.ones((df_loca_Nat_TDT_corresp_filtered.shape[0]))

    path_export_locafile = os.path.join(path_anatomy, f"{sujet}_locafile.xlsx")
    df_loca_Nat_TDT_corresp_filtered.to_excel(path_export_locafile)
    

















################################
######## INSPECT DATA ########
################################



#sujet = sujet_list[2]
def visualize_sig(sujet, step_group=20):

    #### get data
    print('IMPORT')
    data, resp, co2, chanlist, localist = get_data_raw(sujet)

    #### extract length
    chunk_time_list = np.cumsum([_chunk_data.shape[-1] for _chunk_data in data])

    #### generate linear signals
    print('CONSTRUCT LINEAR SIG')
    data_linear = np.concat(data, axis=1)
    resp_linear = np.concat(resp)
    co2_linear = np.concat(co2)

    #### downsample
    time_vec = np.arange(data_linear.shape[-1])
    data_down = scipy.signal.resample_poly(data_linear, up=1, down=10, axis=1, padtype='line')
    respi_down = scipy.signal.resample_poly(resp_linear, up=1, down=10, padtype='line')
    co2_down = scipy.signal.resample_poly(co2_linear, up=1, down=10, padtype='line')
    time_down = scipy.signal.resample_poly(time_vec, up=1, down=10, padtype='line')

    if debug:

        for chan_i in range(10):
            plt.plot(data_down[chan_i])
            plt.show()

        plt.plot(time_down)
        plt.show()

    #### plot
    print("PLOT")
    changroup_plot_list = np.arange(len(chanlist), step=step_group)

    #changroup_plot = changroup_plot_list[3]
    for changroup_plot in changroup_plot_list:

        fig = plotly.graph_objects.Figure()

            #### respi
        fig.add_trace(
                go.Scatter(
                    x=time_down,
                    y=scipy.stats.zscore(respi_down),
                    mode="lines",
                    name='respi'
                )
            )
        
            #### CO2
        y = scipy.stats.zscore(co2_down) - (1 * 4)
        
        fig.add_trace(
                go.Scatter(
                    x=time_down,
                    y=y,
                    mode="lines",
                    name='CO2'
                )
            )

            #### DATA
        chan_sel = np.arange(changroup_plot, changroup_plot+step_group)
        for chan_i_i, chan_i in enumerate(chan_sel):

            y = scipy.stats.zscore(data_down[chan_i]) - ((chan_i_i+2) * 4)

            fig.add_trace(
                go.Scatter(
                    x=time_down,
                    y=y,
                    mode="lines",
                    name=f"{localist[chan_i]} {chanlist[chan_i]}"
                )
            )

        for line in chunk_time_list:
            fig.add_vline(
                x=line,
                line_width=2,
                line_dash="dash",
                line_color="red"
            )

        fig.update_layout(
            template="simple_white",
            title=f"{sujet}",
            xaxis_title="Samples",
            yaxis_title="Amplitude (stacked)",
            height=800
        )

        fig.show()











########################
######## EXPORT ########
########################

#sujet = sujet_list[0]
def export_data(sujet):

    #### get data
    data, resp_clean, co2_corr, chanlist, localist = get_data_raw(sujet)

    path_filetrial = os.path.join(path_data, sujet, f"{sujet}_Experiment_Log.xlsx")
    trial_info = pd.read_excel(path_filetrial, sheet_name='trialsinfo')

    import json
    path_locainfo = os.path.join(path_anatomy, f"{sujet}_SEECOG", "data", f"electrodes.json")
    with open(path_locainfo, "r") as f:
        locadata_raw = json.load(f)

    locaraw = []

    for chan_i, _ in enumerate(locadata_raw):

        locaraw.append(pd.DataFrame({'chan' : [locadata_raw[chan_i]['elecid']], 'loca' : [locadata_raw[chan_i]['anat']]}))

    df_locaraw = pd.concat(locaraw)
    df_locaraw.iloc[:16]

    #### filter data

        #### first filter with whats there is in correspondance sheet and data themselves
    path_export_chanlist = os.path.join(path_anatomy, f"{sujet}_chanlist_info", f"{sujet}_chan_selected.txt")

    with open(path_export_chanlist, "r") as f:
        chanlist_sel = [line.strip() for line in f]

    chanlist_mask = np.array([True if chan in chanlist_sel else False for chan in chanlist])

    data_first_filter = [data_trial[chanlist_mask] for _, data_trial in enumerate(data)]
    chanlist_first_filter = np.array(chanlist)[chanlist_mask]
    localist_first_filter = np.array(localist)[chanlist_mask]

        #### second filter selection after localist inspection
    df_loca = get_chanlocalist(sujet)
    locafile_mask_sig = df_loca['select'].values.astype('bool')
    locafile_mask_anat = (df_loca['Desikan_Killiany_jules'].values != 'Unknown')

    locafile_mask_final = locafile_mask_sig & locafile_mask_anat

    data_second_filter = [data_trial[locafile_mask_final] for _, data_trial in enumerate(data_first_filter)]
    chanlist_second_filter = chanlist_first_filter[locafile_mask_final]
    localist_second_filter = localist_first_filter[locafile_mask_final]

    path_export_chanlist = os.path.join(path_anatomy, f"{sujet}_chanlist_info", f"{sujet}_chanlist.txt")
    with open(path_export_chanlist, "w") as f:
        for _chan in chanlist_second_filter.tolist():
            f.write(str(_chan) + "\n")

    path_export_chanlist = os.path.join(path_anatomy, f"{sujet}_chanlist_info", f"{sujet}_localist.txt")
    with open(path_export_chanlist, "w") as f:
        for _chan in localist_second_filter.tolist():
            f.write(str(_chan) + "\n")

    #### modify trial names
    trial_list = trial_info['Type'].values
    trial_list_correspondance = {'baseline' : 'RB', 'HV' : 'HV', '(1)A0(1)R0C0' : 'AoRoCo', '(1/2)A0(6)R-C0' : 'AoR-Co', 
                                 '(2)A+(1)R0Cc': 'A+RoCc', '(1)A+(6)R-Cc': 'A+R-Cc', '(2)A+(1)R0C-': 'A+RoC-', '(1)A+(6)R-C-': 'A+R-C-'}
    trial_list_count = {'RB' : 0, 'HV' : 0, 'AoRoCo': 0, 'AoR-Co': 0, 'A+RoCc': 0, 'A+R-Cc': 0, 'A+RoC-': 0, 'A+R-C-': 0}
    trial_list_modified = []

    for trial_name in trial_list:

        cond_name = trial_list_correspondance[trial_name]
        trial_list_count[cond_name] += 1
        cond_count = trial_list_count[cond_name]
        trial_list_modified.append(f"{cond_name}_{cond_count}")


    #### export
    path_exportdata = os.path.join(path_prep, sujet)

    for cond_i, cond_name in enumerate(trial_list_modified):

        _data_export = data_second_filter[cond_i]
        path_namedataexport = os.path.join(path_exportdata, f"{sujet}_data_{cond_name}.npy")
        np.save(path_namedataexport, _data_export)

        _respi_export = resp_clean[cond_i]
        path_namedataexport = os.path.join(path_exportdata, f"{sujet}_respi_{cond_name}.npy")
        np.save(path_namedataexport, _respi_export)

        _co2_export = co2_corr[cond_i]
        path_namedataexport = os.path.join(path_exportdata, f"{sujet}_co2_{cond_name}.npy")
        np.save(path_namedataexport, _co2_export)

            

















################################
######## EXECUTE ########
################################

if __name__ == '__main__':

    sujet = 'NS217'
    sujet = 'LH018'

    generate_raw_loca_file(sujet)



