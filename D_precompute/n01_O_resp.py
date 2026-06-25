
from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *

import plotly.graph_objects as go
import seaborn as sns
import neo
import h5py

debug = False






################################
######## CO2 LAG ########
################################


#_respi, y = resp_clean_smooth, CO2
def find_best_lag(_respi, y, srate_sl):

    if debug:

        time_vec = np.arange(_respi.shape[0]) / srate_sl

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
    _peaks_respi, _ = scipy.signal.find_peaks(_respiz, height=_respiz.std(), distance=srate_sl*3)
    
    if debug:

        plt.plot(_respi)
        plt.scatter(_peaks_respi, _respi[_peaks_respi], color='r')
        plt.show()

    _peaks_y, _ = scipy.signal.find_peaks(y_det, height=y_det.std(), distance=srate_sl*3)

    if debug:

        time_vec = np.arange(_respi.shape[0]) / srate_sl

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

        time_vec = np.arange(_respi.shape[0]) / srate_sl

        _shifted_y = np.roll(y, _shift_val_y)
        plt.plot(time_vec, scipy.stats.zscore(_respi), label='respi')
        plt.plot(time_vec, scipy.stats.zscore(y), label=f'raw_y')
        plt.plot(time_vec, scipy.stats.zscore(_shifted_y), label=f'shift_y')
        plt.legend()
        plt.show()

    return _shift_val_y





################################
######## LOAD DATA ########
################################

#sujet = sujet_list[0]
def extract_respfeatures_and_mean_figures(sujet):
    
    #### load respi & CO2
    path_load_trial = os.path.join(path_prep, sujet, 'trial_exports')
    
    #cond = conditions[4]
    for cond in conditions:

        print(f"computing {cond}...")

        cond_trial_list = [f"{sujet}_{cond}_{trial_i+1}.nc" for trial_i in range(ntrail_dict_allpatient[sujet][cond])]
        alltrial_stretch_resp = []
        alltrial_stretch_CO2 = []

        #trial = cond_trial_list[0]
        for trial_i, trial in enumerate(cond_trial_list):

            trial_name = f"{trial.split('_')[1]}_{trial.split('_')[-1][0]}"

            xr_data = xr.load_dataarray(os.path.join(path_load_trial, trial))
            resp = xr_data.sel(chan='respi').data
            CO2 = xr_data.sel(chan='CO2').data

            if debug:

                plt.plot(scipy.stats.zscore(resp), label='resp')
                plt.plot(scipy.stats.zscore(CO2), label='CO2')
                plt.legend()
                plt.show()

            #### extract respfeatures
            resp_clean = physio.preprocess(resp, srate, band=25., btype='lowpass', ftype='bessel', order=5, normalize=False)
            resp_clean_smooth = physio.smooth_signal(resp_clean, srate, win_shape='gaussian', sigma_ms=40.0)

            if debug:
                plt.plot(resp, label='raw')
                plt.plot(resp_clean, label='resp_clean')
                plt.plot(resp_clean_smooth, label='resp_clean_smooth')
                plt.legend()
                plt.show()

            cycles = physio.detect_respiration_cycles(resp_clean_smooth, srate, method="crossing_baseline", baseline_mode='zero')

            resp_features = physio.compute_respiration_cycle_features(resp_clean_smooth, srate, cycles, baseline=None)

            resp_features.insert(0, 'sujet', [sujet]*resp_features.shape[0])
            resp_features.insert(1, 'cond', [cond]*resp_features.shape[0])
            resp_features.insert(2, 'trial', [trial_i+1]*resp_features.shape[0])
            resp_features.insert(3, 'cycle', np.arange(resp_features.shape[0]))

            time_vec = np.arange(resp.shape[0])/srate

            #### CO2 first correction
            CO2_shift_val = find_best_lag(resp_clean_smooth, CO2, srate)
            CO2_shifted = np.roll(CO2, CO2_shift_val)

            if debug:
                plt.plot(scipy.stats.zscore(resp_clean_smooth), label='resp')
                plt.plot(scipy.stats.zscore(CO2), label='CO2')
                plt.plot(scipy.stats.zscore(CO2_shifted), label='CO2_shifted')
                plt.legend()
                plt.show()
            
            #### stretch
            resp_stretch, _ = stretch_data(resp_features, stretch_point_TF, resp_clean_smooth, srate)
            CO2_stretch_before_2nd_correction, _ = stretch_data(resp_features, stretch_point_TF, CO2_shifted, srate)

            if debug:
                for cycle_i in range(resp_stretch.shape[0]):
                    plt.plot(resp_stretch[cycle_i])
                plt.show()

                plt.plot(resp_stretch.mean(axis=0))
                plt.show()

                for cycle_i in range(CO2_stretch_before_2nd_correction.shape[0]):
                    plt.plot(CO2_stretch_before_2nd_correction[cycle_i])
                plt.show()

                plt.plot(CO2_stretch_before_2nd_correction.mean(axis=0))
                plt.show()

            #### CO2 second correction
            CO2_2nd_shift = []

            for cycle_i in range(resp_stretch.shape[0]):

                _resp_cycle = resp_stretch[cycle_i]
                _CO2_cycle = CO2_stretch_before_2nd_correction[cycle_i] 

                if debug:
                    plt.plot(scipy.stats.zscore(_resp_cycle))
                    plt.plot(scipy.stats.zscore(_CO2_cycle))
                    plt.show()

                _shift_i = np.argmax(_resp_cycle) - np.argmax(_CO2_cycle)
                _CO2_cycle_shifted = np.roll(_CO2_cycle, _shift_i)

                if debug:
                    plt.plot(scipy.stats.zscore(_resp_cycle))
                    plt.plot(scipy.stats.zscore(_CO2_cycle_shifted))
                    plt.show()

                CO2_2nd_shift.append(_CO2_cycle_shifted)

            CO2_stretch = np.vstack(CO2_2nd_shift)

            if debug:

                for cycle_i in range(CO2_stretch.shape[0]):
                    plt.plot(CO2_stretch[cycle_i])
                plt.show()

                plt.plot(CO2_stretch.mean(axis=0))
                plt.show()

            #### append to respfeature CO2 amp
            CO2_etCO2 = [CO2_cycle[int(stretch_point_TF/2):].max() for CO2_cycle in CO2_stretch]
            resp_features['etCO2'] = CO2_etCO2

            #### add select col
            select_vec = np.ones((resp_features.index.shape[0]), dtype='int')
            resp_features.insert(resp_features.columns.shape[0], 'select', select_vec)

            #### fig final
            fig_final_resp, ax = plt.subplots(figsize=(18, 10))
            ax.plot(time_vec, resp)
            ax.scatter(cycles[:,0]/srate, resp[cycles[:,0]], color='g', label='inspi_selected')
            ax.scatter(cycles[:,1]/srate, resp[cycles[:,1]], color='c', label='expi_selected', marker='s')
            plt.legend()
            plt.title(f"{sujet} {trial_name} resp")
            plt.rcParams.update({'font.size': 12})
            # plt.show()
            plt.close()

            fig_final_CO2, ax = plt.subplots(figsize=(18, 10))
            ax.plot(time_vec, CO2_shifted)
            ax.scatter(cycles[:,0]/srate, CO2_shifted[cycles[:,0]], color='g', label='inspi_selected')
            ax.scatter(cycles[:,1]/srate, CO2_shifted[cycles[:,1]], color='c', label='expi_selected', marker='s')
            plt.legend()
            plt.title(f"{sujet} {trial_name} CO2")
            plt.rcParams.update({'font.size': 12})
            # plt.show()
            plt.close()

            #### export trial
            path_export_respfeatures = os.path.join(path_precompute, 'RESPI', 'respfeatures', sujet)
            resp_features.to_excel(os.path.join(path_export_respfeatures, f"{sujet}_{trial_name}_respfeatures_RAW.xlsx"))
            fig_final_resp.savefig(os.path.join(path_export_respfeatures, f"{sujet}_{trial_name}_resp_detect_RAW.png"))
            fig_final_CO2.savefig(os.path.join(path_export_respfeatures, f"{sujet}_{trial_name}_CO2_detect_RAW.png"))

            #### append all trial
            alltrial_stretch_resp.append(resp_stretch)
            alltrial_stretch_CO2.append(CO2_stretch)

        #### concat alltrials
        alltrial_stretch_resp = np.concat(alltrial_stretch_resp)
        alltrial_stretch_CO2 = np.concat(alltrial_stretch_CO2)
        
        #### export mean fig
        path_export_mean_fig = os.path.join(path_results, 'respi', 'mean')

        fig_mean, axs = plt.subplots(nrows=2, ncols=2, figsize=(10,8))
        for sig_type_i, sig_type in enumerate(['resp', 'CO2']):

            for plot_type_i, plot_type in enumerate(['mean', 'allcycle']):

                ax = axs[sig_type_i,plot_type_i]

                if sig_type == 'resp':
                    stretch_data_plot = alltrial_stretch_resp
                elif sig_type == 'CO2':
                    stretch_data_plot = alltrial_stretch_CO2

                if plot_type == 'mean':
                
                    _mean = stretch_data_plot.mean(axis=0)
                    _std = stretch_data_plot.std(axis=0)
                    (line,) = ax.plot(_mean)
                    ax.fill_between(range(_mean.size), _mean - _std, _mean + _std, color=line.get_color(), alpha=0.1)
                    ax.set_title(f"{sig_type} {plot_type}")
                    min, max = (_mean-_std).min(), (_mean+_std).max()
                
                elif plot_type == 'allcycle':

                    for cycle_i in range(stretch_data_plot.shape[0]):
                        ax.plot(stretch_data_plot[cycle_i])
                    ax.set_title(f"{sig_type} {plot_type}")
                    min, max = stretch_data_plot.min(), stretch_data_plot.max()
                
                ax.vlines(int(stretch_point_TF/2), ymin=min, ymax=max, color='r')         

        plt.suptitle(f'{sujet} {cond}')
        plt.tight_layout()
        plt.rcParams.update({'font.size': 12})
        # plt.show()

        fig_mean.savefig(os.path.join(path_export_mean_fig, f"{sujet}_{cond}_stretch_mean_RAW.png"))

    #### export respfeature allcond
    path_trigger = os.path.join(path_prep, sujet, f"{sujet}_trigger.xlsx")
    trial_info = pd.read_excel(path_trigger)
    trial_info['cond'] = [correspondance_cond[_cond] for _cond in trial_info['cond'].values]
    count_cond = {cond : 0 for cond in conditions}
    trial_list = []
    for trial in trial_info['cond']:

        count_cond[trial] += 1
        trial_list.append(f"{trial}_{count_cond[trial]}")

    load_list = [f"{sujet}_{file}_respfeatures_RAW.xlsx" for file in trial_list]
    respfeatures_allcond = [pd.read_excel(os.path.join(path_export_respfeatures, file)) for file in load_list]
    respfeatures_allcond = pd.concat(respfeatures_allcond).drop(columns=['Unnamed: 0'])
    respfeatures_allcond = respfeatures_allcond.sort_values(['cond', 'trial'])
    respfeatures_allcond = respfeatures_allcond.reset_index(drop=True)

    cycle_i_vec = []
    cycle_count = 0

    for row_i, row_val in respfeatures_allcond.iterrows():

        if row_i != 0:
            cond_row = row_val['cond']
            if cond_row == cond_prev:
                cycle_count += 1
            else:
                cycle_count = 0
        cycle_i_vec.append(cycle_count)
        cond_prev = row_val['cond']

    respfeatures_allcond['cycle'] = cycle_i_vec

    respfeatures_allcond.to_excel(os.path.join(path_export_respfeatures, f"respfeature_allcond_RAW.xlsx"))






#sujet = sujet_list[0]
def export_cond_relabeling_fig(sujet, calibrating_param='dynamic'):

    #### load respfeatures
    respfeatures_allcond = get_respfeatures_raw(sujet)

    #### params stats
    rf_metrics_fullcycle = ['cycle_duration', 'cycle_freq', 'total_amplitude', 'total_volume', 'etCO2']
    metrics_cycle_features = ['trial', 'cond', 'A_raw', 'R_raw', 'C_raw', 'rf_metric', 'rf_metric_val', 'A_relabel', 'R_relabel', 'C_relabel']
    col_3d_plot_raw = ['cycle_freq', 'total_amplitude', 'etCO2', 'trial', 'cond', 'A_raw', 'R_raw', 'C_raw']
    col_3d_plot_relabeled = ['cycle_freq', 'total_amplitude', 'etCO2', 'trial', 'cond_relabel', 'A_relabel', 'R_relabel', 'C_relabel']

    #### identify threshold
    # sd_factor = 0.5
    # Co_upT = np.median(endtidal_CO2.query(f"cond == 'AoRoCo'")['end_tidal'].values) + sd_factor * scipy.stats.median_abs_deviation(endtidal_CO2.query(f"cond == 'AoRoCo'")['end_tidal'].values)
    # Co_dwT = np.median(endtidal_CO2.query(f"cond == 'AoRoCo'")['end_tidal'].values) - sd_factor * scipy.stats.median_abs_deviation(endtidal_CO2.query(f"cond == 'AoRoCo'")['end_tidal'].values)
    # Ao_upT = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values) + sd_factor * scipy.stats.median_abs_deviation(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values)
    # Ao_dwT = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values) - sd_factor * scipy.stats.median_abs_deviation(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values)
    # Ro_upT = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values) + sd_factor * scipy.stats.median_abs_deviation(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values)
    # Ro_dwT = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values) - sd_factor * scipy.stats.median_abs_deviation(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values)

    # baseline_sel_percentile_range = 60
    # Ao_upT = np.percentile(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values, baseline_sel_percentile_range)
    # Ao_dwT = np.percentile(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values, 0)
    # Ro_upT = np.percentile(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values, 100)
    # Ro_dwT = np.percentile(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values, 100-baseline_sel_percentile_range)
    # Co_upT = np.percentile(respfeatures_allcond.query(f"cond == 'AoRoCo'")['etCO2'].values, 100) 
    # Co_dwT = np.percentile(respfeatures_allcond.query(f"cond == 'AoRoCo'")['etCO2'].values, 100-baseline_sel_percentile_range)

    A_sujet = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values)
    hA_sujet, dA_sujet = A_sujet/2, 2*A_sujet
    
    if calibrating_param == 'solid':
        A_upT = (dA_sujet - A_sujet)*0.50 + A_sujet
        A_dwT = (A_sujet - hA_sujet)*0.50 + hA_sujet
    elif calibrating_param == 'dynamic':
        A_dA_diff = (np.median(respfeatures_allcond.query(f"cond == 'A+RoCc'")['total_amplitude'].values) - np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values))
        A_upT = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values) + A_dA_diff / 2
        A_hA_diff = (np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values) - np.median(respfeatures_allcond.query(f"cond == 'AoR-Co'")['total_amplitude'].values))
        A_dwT = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values) - A_hA_diff / 2

    R_sujet = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values)
    if calibrating_param == 'solid':
        R_thresh = (R_sujet - 0.1)*0.50 + 0.1
    elif calibrating_param == 'dynamic':
        R_hR_diff = (np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values) - np.median(respfeatures_allcond.query(f"cond == 'AoR-Co'")['cycle_freq'].values))
        R_thresh = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values) - R_hR_diff / 2

    C_sujet = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['etCO2'].values)
    C_hypo = np.median(respfeatures_allcond.query(f"cond == 'A+RoC-'")['etCO2'].values)
    C_thresh = (C_sujet - C_hypo)*0.50 + C_hypo

    if sujet == 'NS217':
        A_sujet /= 2
        hA_sujet /= 2
        dA_sujet /= 2
        A_dwT /= 2
        A_upT /= 2

        C_thresh = np.median(respfeatures_allcond['etCO2'].values)

    if debug:

        med = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values)
        std = respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values.std()

        rf_metric = 'total_amplitude'
        data_plot = respfeatures_allcond.query(f"cond == 'AoRoCo'")[rf_metric].values
        count, _, _ = plt.hist(data_plot, bins=50, label='AoRoCo', alpha=0.5)
        data_plot = respfeatures_allcond.query(f"cond == 'A+RoCc'")[rf_metric].values
        count, _, _ = plt.hist(data_plot, bins=50, label='A+RoCc', alpha=0.5)
        data_plot = respfeatures_allcond.query(f"cond == 'AoR-Co'")[rf_metric].values
        count, _, _ = plt.hist(data_plot, bins=50, label='A+AoR-Cc', alpha=0.5)
        plt.vlines([hA_sujet, A_sujet, dA_sujet], ymin=0, ymax=count.max(), color='k', label='hA_sujet, A_sujet, dA_sujet')
        plt.vlines([A_dwT, A_upT], ymin=0, ymax=count.max(), color='k', linestyles='--', label='A_dwT, A_upT')
        plt.vlines([med-2*std, med+2*std], ymin=0, ymax=count.max(), color='r', linestyles='-', label='sd_dwT, sd_upT')
        plt.title(calibrating_param)
        plt.legend()
        plt.show()

        rf_metric = 'cycle_freq'
        data_plot = respfeatures_allcond.query(f"cond == 'AoRoCo'")[rf_metric].values
        count, _, _ = plt.hist(data_plot, bins=50, label='AoRoCo', alpha=0.5)
        data_plot = respfeatures_allcond.query(f"cond == 'AoR-Co'")[rf_metric].values
        count, _, _ = plt.hist(data_plot, bins=50, label='AoR-Co', alpha=0.5)
        plt.vlines([R_sujet], ymin=0, ymax=count.max(), color='g', label='R_sujet')
        plt.vlines([R_thresh], ymin=0, ymax=count.max(), color='g', linestyles='--', label='R_sujet')
        plt.legend()
        plt.show()

        rf_metric = 'etCO2'
        data_plot = respfeatures_allcond.query(f"cond == 'AoRoCo'")[rf_metric].values
        count, _, _ = plt.hist(data_plot, bins=50, label='AoRoCo', alpha=0.5)
        data_plot = respfeatures_allcond.query(f"cond == 'A+RoC-'")[rf_metric].values
        count, _, _ = plt.hist(data_plot, bins=50, label='A+RoC-', alpha=0.5)
        plt.vlines([C_sujet], ymin=0, ymax=count.max(), color='g', label='C_sujet')
        plt.vlines([C_thresh], ymin=0, ymax=count.max(), color='g', linestyles='--', label='C_thresh')
        plt.legend()
        plt.show()

    #### apply thresh
    relabel_R = (respfeatures_allcond['cycle_freq'] < R_thresh).replace({True : '-', False : 'o'}).values
    relabel_A = np.full((relabel_R.size), fill_value='', dtype='object')
    sel_05R = (respfeatures_allcond['total_amplitude'] < A_dwT).values
    relabel_A[sel_05R] = '05A'
    sel_2R = (respfeatures_allcond['total_amplitude'] > A_upT).values
    relabel_A[sel_2R] = '2A'
    sel_A = ~(sel_05R | sel_2R)
    relabel_A[sel_A] = 'A'
    relabel_C = (respfeatures_allcond['etCO2'] < C_thresh).replace({True : '-', False : 'o'}).values
    relabel_code = np.concat([relabel_A.reshape(-1,1), relabel_R.reshape(-1,1), relabel_C.reshape(-1,1)], axis=1)

    cond_code = {'AoRoCo' : ['A','o','o'], 'AoR-Co' : ['05A','-','o'], 'A+RoCc' : ['2A','o','o'], 
                 'A+R-Cc' : ['A','-','o'], 'A+RoC-' : ['2A','o','-'], 'A+R-C-' : ['A','-','-']}
    
    respfeatures_allcond_labeled = respfeatures_allcond.copy()
    respfeatures_allcond_labeled['A_relabel'] = relabel_A
    respfeatures_allcond_labeled['R_relabel'] = relabel_R
    respfeatures_allcond_labeled['C_relabel'] = relabel_C
    
    cond_relabeled = np.array([''] * respfeatures_allcond_labeled.shape[0], dtype='object')
    exclude_cycle_map = np.full((respfeatures_allcond_labeled.shape[0]), fill_value=False)

    for cond, label_code in cond_code.items():

        sel_vec = (relabel_code == label_code).all(axis=1)
        cond_relabeled[sel_vec] = cond
        exclude_cycle_map = exclude_cycle_map | sel_vec

    cond_relabeled[~exclude_cycle_map] = 'excluded'

    if debug:

        df_plot = pd.concat([
            pd.DataFrame({'cond': respfeatures_allcond_labeled['cond'], 'type': 'raw'}),
            pd.DataFrame({'cond': respfeatures_allcond_labeled['cond_relabel'], 'type': 'relabel'})
        ])

        plt.figure(figsize=(12, 6))

        sns.histplot(data=df_plot, x='cond', hue='type', multiple='dodge', shrink=.8)
        plt.axhline(30, color='red', linestyle='--')

        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.show()

    respfeatures_allcond_labeled['cond_relabel'] = cond_relabeled
    respfeatures_allcond_labeled['A_relabel'] = respfeatures_allcond_labeled['A_relabel'].replace({False : 'o', True : '+'})
    respfeatures_allcond_labeled['R_relabel'] = respfeatures_allcond_labeled['R_relabel'].replace({False : 'o', True : '-'})
    respfeatures_allcond_labeled['C_relabel'] = respfeatures_allcond_labeled['C_relabel'].replace({False : 'o', True : '-'})

    #### export respfeatures relabelled
    path_export_respfeatures = os.path.join(path_precompute, 'RESPI', 'respfeatures', sujet)
    respfeatures_allcond_labeled.to_excel(os.path.join(path_export_respfeatures, f"respfeature_allcond_relabel.xlsx"))

    ######## PLOT ########

    #### path
    path_export_label = os.path.join(path_results, 'respi', 'relabel')

    #### export relabel modification plot
    df_plot = pd.concat([
            pd.DataFrame({'cond': respfeatures_allcond_labeled['cond'], 'type': 'raw'}),
            pd.DataFrame({'cond': respfeatures_allcond_labeled['cond_relabel'], 'type': 'relabel'})
        ])

    fig_count = plt.figure(figsize=(10, 6))

    sns.histplot(data=df_plot, x='cond', hue='type', multiple='dodge', shrink=.8)
    plt.axhline(30, color='red', linestyle='--')

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.suptitle(calibrating_param)
    #plt.show()
    
    fig_count.savefig(os.path.join(path_export_label, f"relabeled_{sujet}.png"))

    #### export stats plot raw
    A_raw, R_raw, C_raw = [], [], []
    for _cond in respfeatures_allcond_labeled['cond']:
        if _cond in ['RB', 'HV']:            
            A_raw.append(_cond)
            R_raw.append(_cond)
            C_raw.append(_cond)
        else:
            A_raw.append(_cond[1])
            R_raw.append(_cond[3])
            C_raw.append(_cond[5])
        
    respfeatures_allcond_labeled['A_raw'], respfeatures_allcond_labeled['R_raw'], respfeatures_allcond_labeled['C_raw'] = A_raw, R_raw, C_raw

    df_plot = respfeatures_allcond_labeled.melt(id_vars=[c for c in respfeatures_allcond_labeled.columns if c not in rf_metrics_fullcycle],
                                    value_vars=rf_metrics_fullcycle, var_name="rf_metric", value_name="rf_metric_val")[metrics_cycle_features]
    
    A_raw_update = []
    for _cond in df_plot['cond']:
        if _cond in ['RB', 'HV']:
            A_raw_update.append(_cond)
        else:
            A_raw_update.append(cond_code[_cond][0])
    df_plot['A_raw'] = A_raw_update   
    
    df_plot_raw = df_plot.melt(id_vars=[c for c in df_plot.columns if c not in ['A_raw', 'R_raw', 'C_raw']],
                                    value_vars=['A_raw', 'R_raw', 'C_raw'], var_name="cond_label", value_name="cond_label_val")
    df_plot_raw['cond_label'] = df_plot_raw['cond_label'].replace({'A_raw' : 'A', 'R_raw' : 'R', 'C_raw' : 'C'})
    df_plot_raw['type'] = ['raw'] * df_plot_raw.shape[0]
    
    df_plot_rellabeled = df_plot.melt(id_vars=[c for c in df_plot.columns if c not in ['A_relabel', 'R_relabel', 'C_relabel']],
                                    value_vars=['A_relabel', 'R_relabel', 'C_relabel'], var_name="cond_label", value_name="cond_label_val")
    df_plot_rellabeled['cond_label'] = df_plot_rellabeled['cond_label'].replace({'A_relabel' : 'A', 'R_relabel' : 'R', 'C_relabel' : 'C'})
    df_plot_rellabeled['type'] = ['rellabeled'] * df_plot_rellabeled.shape[0]

    df_plot = pd.concat([df_plot_raw, df_plot_rellabeled]).reset_index(drop=True)
    df_plot = df_plot.query(f"cond_label_val not in ['RB', 'HV']")
    

    fig_raw_stats_fullcycle = sns.catplot(data=df_plot, kind='strip', x='cond_label_val', y='rf_metric_val', col='rf_metric', 
                                          row='cond_label', hue='type', sharey=False, dodge=True)
    fig_raw_stats_fullcycle.figure.suptitle(sujet)
    plt.tight_layout()
    # plt.show()

    fig_raw_stats_fullcycle.savefig(os.path.join(path_export_label, f"fullcycle_stats_{sujet}.png"))

    #### 3d params
    color_cond = {
    'AoRoCo': 'red',
    'AoR-Co': 'blue',
    'A+RoCc': 'black',
    'A+R-Cc': 'cyan',
    'A+RoC-': 'green',
    'A+R-C-': 'yellow'
    }

    #### 3d raw cond sel
    df_plot = respfeatures_allcond_labeled[col_3d_plot_raw]

    fig_raw = go.Figure()

    cond_sel = [cond for cond in conditions if cond not in ['RB', 'HV']]

    for cond in cond_sel:
        _sel_vec = df_plot['cond'] == cond
        _sel_CO2 = df_plot['etCO2'].values[_sel_vec]
        _sel_amp = df_plot['total_amplitude'].values[_sel_vec]
        _sel_freq = df_plot['cycle_freq'].values[_sel_vec]

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
    fig_raw.write_html(os.path.join(path_export_label, f"3d_{sujet}_allcond_raw.html"))

    #### 3d relabeled cond sel
    df_plot = respfeatures_allcond_labeled[col_3d_plot_relabeled]

    cond_sel = [cond for cond in conditions if cond not in ['RB', 'HV']]

    fig_relabeled = go.Figure()

    for cond in cond_sel:
        _sel_vec = df_plot['cond_relabel'] == cond
        _sel_CO2 = df_plot['etCO2'].values[_sel_vec]
        _sel_amp = df_plot['total_amplitude'].values[_sel_vec]
        _sel_freq = df_plot['cycle_freq'].values[_sel_vec]

        fig_relabeled.add_trace(go.Scatter3d(
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

    fig_relabeled.update_layout(
        scene=dict(
            xaxis_title="CO2",
            yaxis_title="AMP",
            zaxis_title="FREQ"
        ),
        legend=dict(title="Condition"),
        width=900,
        height=700,
        title={'text':"RELABELED", 'x':0.5, 'xanchor':"center", 'yanchor':"top", 'y':0.95}
    )

    # fig_relabeled.show()
    fig_relabeled.write_html(os.path.join(path_export_label, f"3d_{sujet}_allcond_relabeled.html"))







################################
######## SCALE ########
################################


def export_rf_cond_scaling(sujet):

    #### extract respfeatures
    respfeature = get_respfeatures_relabel_raw(sujet)

    #### extract ref
    _rf_baseline_VCCP = respfeature.query(f"cond_relabel == 'AoRoCo'")
    _rf_baseline_RB = respfeature.query(f"cond == 'RB'")
    A_VCCP_SD, R_VCCP_SD, C_VCCP_SD = _rf_baseline_VCCP['total_amplitude'].std(), _rf_baseline_VCCP['cycle_freq'].std(), _rf_baseline_VCCP['etCO2'].std() 
    A_RB_SD, R_RB_SD, C_RB_SD = _rf_baseline_RB['total_amplitude'].std(), _rf_baseline_RB['cycle_freq'].std(), _rf_baseline_RB['etCO2'].std() 
    
    #### get rf_scaled
    rf_SD_scaled_VCCP = respfeature.copy()
    rf_SD_scaled_RB = respfeature.copy()

    df_list = [rf_SD_scaled_VCCP, rf_SD_scaled_RB]
    df_type = ['VCCP', 'RB']
    df_scaled = []

    for _df_i, _df in enumerate(df_list):
        
        if df_type[_df_i] == 'VCCP':
            _A_baseline, _R_baseline, _C_baseline = A_VCCP_SD, R_VCCP_SD, C_VCCP_SD
        else:
            _A_baseline, _R_baseline, _C_baseline = A_RB_SD, R_RB_SD, C_RB_SD
    
        _df['A_scaled'] = _df['total_amplitude'] / _A_baseline
        _df['R_scaled'] = _df['cycle_freq'] / _R_baseline
        _df['C_scaled'] = _df['etCO2'] / _C_baseline

        df_scaled.append(_df)

    rf_SD_scaled_VCCP_export, rf_SD_scaled_RB_export = df_scaled[0], df_scaled[1]

    #### export
    path_export = os.path.join(path_precompute, 'RESPI', 'respfeatures', sujet)
    rf_SD_scaled_VCCP_export.to_excel(os.path.join(path_export, f"{sujet}_df_scaled_VCCP.xlsx"))
    rf_SD_scaled_RB_export.to_excel(os.path.join(path_export, f"{sujet}_df_scaled_RB.xlsx"))





################################
######## EXECUTE ########
################################

if __name__ == '__main__':

    #sujet = 'LH018'
    for sujet in sujet_list:

        extract_respfeatures_and_mean_figures(sujet)
        export_cond_relabeling_fig(sujet)
        export_rf_cond_scaling(sujet)















