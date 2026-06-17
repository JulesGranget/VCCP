
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
    trial_list = os.listdir(path_load_trial)

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

            resp_features.insert(0, 'cond', [cond]*resp_features.shape[0])
            resp_features.insert(1, 'trial', [trial_i+1]*resp_features.shape[0])
            select_vec = np.ones((resp_features.index.shape[0]), dtype='int')
            resp_features.insert(resp_features.columns.shape[0], 'select', select_vec)

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
            resp_features.to_excel(os.path.join(path_export_respfeatures, f"{sujet}_{trial_name}_respfeatures.xlsx"))
            fig_final_resp.savefig(os.path.join(path_export_respfeatures, f"{sujet}_{trial_name}_resp_detect.png"))
            fig_final_CO2.savefig(os.path.join(path_export_respfeatures, f"{sujet}_{trial_name}_CO2_detect.png"))

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

        plt.suptitle(f'{sujet}')
        plt.tight_layout()
        plt.rcParams.update({'font.size': 12})
        # plt.show()

        fig_mean.savefig(os.path.join(path_export_mean_fig, f"{sujet}_stretch_mean.png"))

    #### export respfeature allcond
    listdir_respfeatures = os.listdir(path_export_respfeatures)
    load_list = [file for file in listdir_respfeatures if file.find('.xlsx') != -1]
    respfeatures_allcond = [pd.read_excel(os.path.join(path_export_respfeatures, file)) for file in load_list]
    respfeatures_allcond = pd.concat(respfeatures_allcond).drop(columns=['Unnamed: 0'])
    respfeatures_allcond = respfeatures_allcond.reset_index(drop=True)
    respfeatures_allcond.to_excel(os.path.join(path_export_respfeatures, f"respfeature_allcond.xlsx"))






#sujet = sujet_list[0]
def export_cond_relabeling_fig(sujet):

    #### load respfeatures
    respfeatures_allcond = get_respfeatures(sujet)

    #### get cond labels
    A_label = []
    R_label = []
    C_label = []

    for cond in respfeatures_allcond_raw_import['cond']:

        A_label.append(cond[1])
        R_label.append(cond[3])
        C_label.append(cond[5])
        
    cond_label = pd.DataFrame({'A_raw' : A_label, 'R_raw' : R_label, 'C_raw' : C_label})
    respfeatures_allcond_labeled_raw = pd.concat([respfeatures_allcond_raw_import.reset_index(), cond_label], axis=1).drop(columns=['index', 'Unnamed: 0'])

    #### params stats
    rf_metrics_bothphase = ['cycle_duration', 'inspi_duration', 'expi_duration', 'cycle_freq', 'cycle_ratio', 'inspi_volume',
       'expi_volume', 'total_amplitude', 'inspi_amplitude', 'expi_amplitude','total_volume', 'CO2_amp']
    rf_metrics_fullcycle = ['cycle_duration', 'cycle_freq', 'total_amplitude', 'total_volume', 'CO2_amp']
    info_columns_raw = ['trial', 'cond', 'A_raw', 'R_raw', 'C_raw', 'rf_metric', 'rf_metric_val']
    info_columns_relabel = ['trial', 'rf_metric', 'rf_metric_val', 'cond_relabel', 'A_relabel', 'R_relabel', 'C_relabel']
    col_3d_plot_raw = ['cycle_freq', 'total_amplitude', 'CO2_amp', 'trial', 'cond', 'A_raw', 'R_raw', 'C_raw']
    col_3d_plot_relabeled = ['cycle_freq', 'total_amplitude', 'CO2_amp', 'trial', 'cond_relabel', 'A_relabel', 'R_relabel', 'C_relabel']

    #### relabelling
    # sd_factor = 0.5
    # Co_upT = np.median(endtidal_CO2.query(f"cond == 'AoRoCo'")['end_tidal'].values) + sd_factor * scipy.stats.median_abs_deviation(endtidal_CO2.query(f"cond == 'AoRoCo'")['end_tidal'].values)
    # Co_dwT = np.median(endtidal_CO2.query(f"cond == 'AoRoCo'")['end_tidal'].values) - sd_factor * scipy.stats.median_abs_deviation(endtidal_CO2.query(f"cond == 'AoRoCo'")['end_tidal'].values)
    # Ao_upT = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values) + sd_factor * scipy.stats.median_abs_deviation(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values)
    # Ao_dwT = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values) - sd_factor * scipy.stats.median_abs_deviation(respfeatures_allcond.query(f"cond == 'AoRoCo'")['total_amplitude'].values)
    # Ro_upT = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values) + sd_factor * scipy.stats.median_abs_deviation(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values)
    # Ro_dwT = np.median(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values) - sd_factor * scipy.stats.median_abs_deviation(respfeatures_allcond.query(f"cond == 'AoRoCo'")['cycle_freq'].values)

    # baseline_sel_percentile_range = 60
    # Ao_upT = np.percentile(respfeatures_allcond_labeled_raw.query(f"cond == 'AoRoCo'")['total_amplitude'].values, baseline_sel_percentile_range)
    # Ao_dwT = np.percentile(respfeatures_allcond_labeled_raw.query(f"cond == 'AoRoCo'")['total_amplitude'].values, 0)
    # Ro_upT = np.percentile(respfeatures_allcond_labeled_raw.query(f"cond == 'AoRoCo'")['cycle_freq'].values, 100)
    # Ro_dwT = np.percentile(respfeatures_allcond_labeled_raw.query(f"cond == 'AoRoCo'")['cycle_freq'].values, 100-baseline_sel_percentile_range)
    # Co_upT = np.percentile(respfeatures_allcond_labeled_raw.query(f"cond == 'AoRoCo'")['CO2_amp'].values, 100) 
    # Co_dwT = np.percentile(respfeatures_allcond_labeled_raw.query(f"cond == 'AoRoCo'")['CO2_amp'].values, 100-baseline_sel_percentile_range)

    A_sujet = np.median(respfeatures_allcond_labeled_raw.query(f"cond == 'AoRoCo'")['total_amplitude'].values)
    hA_sujet, dA_sujet = A_sujet/2, 2*A_sujet
    A_upT = (dA_sujet - A_sujet)*0.50 + A_sujet
    A_dwT = (A_sujet - hA_sujet)*0.50 + hA_sujet

    R_sujet = np.median(respfeatures_allcond_labeled_raw.query(f"cond == 'AoRoCo'")['cycle_freq'].values)
    R_thresh = (R_sujet - 0.1)*0.50 + 0.1

    C_sujet = np.median(respfeatures_allcond_labeled_raw.query(f"cond == 'AoRoCo'")['CO2_amp'].values)
    C_hypo = np.median(respfeatures_allcond_labeled_raw.query(f"cond == 'A+RoC-'")['CO2_amp'].values)
    C_thresh = (C_sujet - C_hypo)*0.50 + C_hypo

    if sujet == 'NS217':
        A_sujet /= 2
        hA_sujet /= 2
        dA_sujet /= 2
        A_dwT /= 2
        A_upT /= 2

        C_thresh = np.median(respfeatures_allcond_labeled_raw['CO2_amp'].values)

    if debug:

        rf_metric = 'total_amplitude'
        data_plot = respfeatures_allcond_labeled_raw.query(f"cond == 'AoRoCo'")[rf_metric].values
        count, _, _ = plt.hist(data_plot, bins=50, label='AoRoCo', alpha=0.5)
        data_plot = respfeatures_allcond_labeled_raw.query(f"cond == 'A+RoCc'")[rf_metric].values
        count, _, _ = plt.hist(data_plot, bins=50, label='A+RoCc', alpha=0.5)
        plt.vlines([hA_sujet, A_sujet, dA_sujet], ymin=0, ymax=count.max(), color='g', label='hA_sujet, A_sujet, dA_sujet')
        plt.vlines([A_dwT, A_upT], ymin=0, ymax=count.max(), color='g', linestyles='--', label='A_dwT, A_upT')
        plt.legend()
        plt.show()

        rf_metric = 'cycle_freq'
        data_plot = respfeatures_allcond_labeled_raw.query(f"cond == 'AoRoCo'")[rf_metric].values
        count, _, _ = plt.hist(data_plot, bins=50, label='AoRoCo', alpha=0.5)
        data_plot = respfeatures_allcond_labeled_raw.query(f"cond == 'AoR-Co'")[rf_metric].values
        count, _, _ = plt.hist(data_plot, bins=50, label='AoR-Co', alpha=0.5)
        plt.vlines([R_sujet], ymin=0, ymax=count.max(), color='g', label='R_sujet')
        plt.vlines([R_thresh], ymin=0, ymax=count.max(), color='g', linestyles='--', label='R_sujet')
        plt.legend()
        plt.show()

        rf_metric = 'CO2_amp'
        data_plot = respfeatures_allcond_labeled_raw.query(f"cond == 'AoRoCo'")[rf_metric].values
        count, _, _ = plt.hist(data_plot, bins=50, label='AoRoCo', alpha=0.5)
        data_plot = respfeatures_allcond_labeled_raw.query(f"cond == 'A+RoC-'")[rf_metric].values
        count, _, _ = plt.hist(data_plot, bins=50, label='A+RoC-', alpha=0.5)
        plt.vlines([C_sujet], ymin=0, ymax=count.max(), color='g', label='C_sujet')
        plt.vlines([C_thresh], ymin=0, ymax=count.max(), color='g', linestyles='--', label='C_thresh')
        plt.legend()
        plt.show()

    
    relabel_R = (respfeatures_allcond_labeled_raw['cycle_freq'] < R_thresh).replace({True : '-', False : 'o'}).values
    relabel_A = np.full((relabel_R.size), fill_value='', dtype='object')
    sel_05R = (respfeatures_allcond_labeled_raw['total_amplitude'] < A_dwT).values
    relabel_A[sel_05R] = '05A'
    sel_2R = (respfeatures_allcond_labeled_raw['total_amplitude'] > A_upT).values
    relabel_A[sel_2R] = '2A'
    sel_A = ~(sel_05R | sel_2R)
    relabel_A[sel_A] = 'A'
    relabel_C = (respfeatures_allcond_labeled_raw['CO2_amp'] < C_thresh).replace({True : '-', False : 'o'}).values
    relabel_code = np.concat([relabel_A.reshape(-1,1), relabel_R.reshape(-1,1), relabel_C.reshape(-1,1)], axis=1)

    cond_code = {'AoRoCo' : ['A','o','o'], 'AoR-Co' : ['05A','-','o'], 'A+RoCc' : ['2A','o','o'], 
                 'A+R-Cc' : ['A','-','o'], 'A+RoC-' : ['2A','o','-'], 'A+R-C-' : ['A','-','-']}
    
    respfeatures_allcond_labeled = respfeatures_allcond_labeled_raw.copy()
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

    respfeatures_allcond_labeled['cond_relabel'] = cond_relabeled
    respfeatures_allcond_labeled['A_relabel'] = respfeatures_allcond_labeled['A_relabel'].replace({False : 'o', True : '+'})
    respfeatures_allcond_labeled['R_relabel'] = respfeatures_allcond_labeled['R_relabel'].replace({False : 'o', True : '-'})
    respfeatures_allcond_labeled['C_relabel'] = respfeatures_allcond_labeled['C_relabel'].replace({False : 'o', True : '-'})

    os.chdir(os.path.join(path_results, 'respi', 'respfeatures', sujet))
    respfeatures_allcond_labeled.to_excel(f"{sujet}_allcond_respfeatures_relabel.xlsx")

    ######## PLOT ########
    #### export relabel modification plot
    counts_cond = respfeatures_allcond_labeled["cond"].value_counts().rename("count").reset_index()
    counts_cond_relabel = respfeatures_allcond_labeled["cond_relabel"].value_counts().rename("count").reset_index()
    counts_cond_relabel = counts_cond_relabel.rename(columns={'cond_relabel' : 'cond'})

    counts_cond["source"] = "cond"
    counts_cond_relabel["source"] = "cond_relabel"

    df_counts = pd.concat([counts_cond, counts_cond_relabel], ignore_index=True)

    g = sns.catplot(df_counts, kind='bar', x='cond', y='count', hue='source')
    plt.title(sujet)
    plt.tight_layout()
    # plt.show()
    
    os.chdir(os.path.join(path_results, 'respi', 'stats_cycle_select'))
    g.savefig(f"relabeled_{sujet}.png")

    #### export stats plot raw
    df_plot = respfeatures_allcond_labeled.melt(id_vars=[c for c in respfeatures_allcond_labeled.columns if c not in rf_metrics_fullcycle],
                                    value_vars=rf_metrics_fullcycle, var_name="rf_metric", value_name="rf_metric_val")[info_columns_raw]
    df_plot = df_plot.melt(id_vars=[c for c in df_plot.columns if c not in ['A_raw', 'R_raw', 'C_raw']],
                                    value_vars=['A_raw', 'R_raw', 'C_raw'], var_name="cond_label", value_name="cond_label_val")

    fig_raw_stats_fullcycle = sns.catplot(data=df_plot, kind='swarm', x='cond_label_val', y='rf_metric_val', col='rf_metric', row='cond_label', sharey=False)
    fig_raw_stats_fullcycle.figure.suptitle(sujet)
    plt.tight_layout()
    # plt.show()

    os.chdir(os.path.join(path_results, 'respi', 'stats_cycle_select'))
    fig_raw_stats_fullcycle.savefig(f"fullcycle_stats_{sujet}_raw.png")

    #### export stats plot relabel
    df_plot = respfeatures_allcond_labeled.melt(id_vars=[c for c in respfeatures_allcond_labeled.columns if c not in rf_metrics_fullcycle],
                                    value_vars=rf_metrics_fullcycle, var_name="rf_metric", value_name="rf_metric_val")[info_columns_relabel]
    df_plot = df_plot.melt(id_vars=[c for c in df_plot.columns if c not in ['A_relabel', 'R_relabel', 'C_relabel']],
                                    value_vars=['A_relabel', 'R_relabel', 'C_relabel'], var_name="cond_label", value_name="cond_label_val")

    fig_relabel_stats_fullcycle = sns.catplot(data=df_plot, kind='swarm', x='cond_label_val', y='rf_metric_val', col='rf_metric', row='cond_label', sharey=False)
    fig_relabel_stats_fullcycle.figure.suptitle(sujet)
    plt.tight_layout()
    # plt.show()

    os.chdir(os.path.join(path_results, 'respi', 'stats_cycle_select'))
    fig_relabel_stats_fullcycle.savefig(f"fullcycle_stats_{sujet}_relabel.png")

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

    for cond in conditions:
        _sel_vec = df_plot['cond'] == cond
        _sel_CO2 = df_plot['CO2_amp'].values[_sel_vec]
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
    os.chdir(os.path.join(path_results, 'respi', 'stats_cycle_select'))
    fig_raw.write_html(f"3d_{sujet}_allcond_raw.html")

    #### 3d relabeled cond sel
    df_plot = respfeatures_allcond_labeled[col_3d_plot_relabeled]

    fig_relabeled = go.Figure()

    for cond in conditions:
        _sel_vec = df_plot['cond_relabel'] == cond
        _sel_CO2 = df_plot['CO2_amp'].values[_sel_vec]
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
    os.chdir(os.path.join(path_results, 'respi', 'stats_cycle_select'))
    fig_relabeled.write_html(f"3d_{sujet}_allcond_relabeled.html")










################################
######## EXECUTE ########
################################

if __name__ == '__main__':

    sujet = 'NS217'
    sujet = 'LH018'

    extract_respfeatures_and_mean_figures(sujet)
    export_cond_relabeling_fig(sujet)















