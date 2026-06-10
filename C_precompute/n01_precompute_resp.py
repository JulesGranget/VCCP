
from config.n00_O_config_params import *
from config.n00bis_O_config_analysis_functions import *
from config.n01bis_patient_info import *

import plotly.graph_objects as go
import seaborn as sns
import neo
import h5py

debug = False



################################
######## LOAD DATA ########
################################


def extract_respfeatures_and_mean_figures(sujet):
    
    #### load data
    os.chdir(os.path.join(path_prep, f"{sujet}_exports"))

    cond_list = []
    resp = []
    CO2 = []
    for cond in conditions:

        for trial_i in range(ntrail_dict_allpatient[sujet][cond]):

            _resp = np.load(f"{sujet}_{cond}0{trial_i+1}_aux.npy")[np.where(np.array(aux_chanlist) == 'resp')[0][0]]
            _CO2 = np.load(f"{sujet}_{cond}0{trial_i+1}_aux.npy")[np.where(np.array(aux_chanlist) == 'co2')[0][0]]
            cond_list.append(f"{cond}0{trial_i+1}")
            resp.append(_resp)
            CO2.append(_CO2)

    #### extract respfeatures
    respfeatures = []
    fig_resp_detect_resp = []
    fig_resp_detect_co2 = []
    
    for resp_i, resp_sig in enumerate(resp):
        
        cycles = physio.detect_respiration_cycles(resp_sig, srate, method="crossing_baseline", baseline_mode='zero')

        _resp_features = physio.compute_respiration_cycle_features(resp_sig, srate, cycles, baseline=None)

        select_vec = np.ones((_resp_features.index.shape[0]), dtype='int')
        _resp_features.insert(_resp_features.columns.shape[0], 'select', select_vec)

        time_vec = np.arange(resp_sig.shape[0])/srate

        #### fig final
        fig_final_resp, ax = plt.subplots(figsize=(18, 10))
        ax.plot(time_vec, resp_sig)
        ax.scatter(cycles[:,0]/srate, resp_sig[cycles[:,0]], color='g', label='inspi_selected')
        ax.scatter(cycles[:,1]/srate, resp_sig[cycles[:,1]], color='c', label='expi_selected', marker='s')
        plt.legend()
        plt.title(f"{sujet} {cond_list[resp_i]} resp")
        plt.rcParams.update({'font.size': 12})
        # plt.show()
        plt.close()

        _CO2 = CO2[resp_i]
        fig_final_CO2, ax = plt.subplots(figsize=(18, 10))
        ax.plot(time_vec, _CO2)
        ax.scatter(cycles[:,0]/srate, _CO2[cycles[:,0]], color='g', label='inspi_selected')
        ax.scatter(cycles[:,1]/srate, _CO2[cycles[:,1]], color='c', label='expi_selected', marker='s')
        plt.legend()
        plt.title(f"{sujet} {cond_list[resp_i]} CO2")
        plt.rcParams.update({'font.size': 12})
        # plt.show()
        plt.close()
        
        respfeatures.append(_resp_features)
        fig_resp_detect_resp.append(fig_final_resp)
        fig_resp_detect_co2.append(fig_final_CO2)

    #### stretch
    resp_stretch_allcond = {}
    co2_stretch_allcond = {}

    co2_amp = []

    for cond in conditions:

        stretch_respi = []
        stretch_co2 = []

        for trial_i in range(ntrail_dict_allpatient[sujet][cond]):

            _trial_name = f"{cond}0{trial_i+1}"
            _trial_i = np.where(np.array(cond_list) == _trial_name)[0][0]
            _resp, _respfeature = resp[_trial_i], respfeatures[_trial_i]
            _co2 = CO2[_trial_i]

            _resp_stretch, _ = stretch_data(_respfeature, nb_point_by_cycle, _resp, srate)
            _co2_stretch, _ = stretch_data(_respfeature, nb_point_by_cycle, _co2, srate)

            if debug:
                plt.plot(_resp_stretch.mean(axis=0))
                plt.show()

                plt.plot(_co2_stretch.mean(axis=0))
                plt.show()

            stretch_respi.append(_resp_stretch)
            stretch_co2.append(_co2_stretch)

            co2_amp.append(np.abs(_co2_stretch.min(axis=1)) + np.abs(_co2_stretch.max(axis=1)))

        resp_stretch_allcond[cond] = np.vstack(stretch_respi)
        co2_stretch_allcond[cond] = np.vstack(stretch_co2)

    #### append to respfeature co2 amp
    respfeatures_with_co2 = []

    for respfeature_trial_i, respfeatures_trial in enumerate(respfeatures):

        _co2_amp = co2_amp[respfeature_trial_i]
        respfeatures_trial['CO2_amp'] = _co2_amp
        respfeatures_with_co2.append(respfeatures_trial)
        
    #### export
    for trial_i, cond in enumerate(cond_list):

        os.chdir(os.path.join(path_results, 'respi', 'respfeatures', sujet))
        respfeatures[trial_i].to_excel(f"{sujet}_{cond}_respfeatures.xlsx")

        os.chdir(os.path.join(path_results, 'respi', 'fig_detect', sujet))
        fig_resp_detect_resp[trial_i].savefig(f"{sujet}_{cond}_resp.png")
        fig_resp_detect_co2[trial_i].savefig(f"{sujet}_{cond}_co2.png")

    #### export mean fig
    os.chdir(os.path.join(path_results, 'respi', 'summary', 'patient_wise'))

    fig_mean_resp, ax = plt.subplots()
    for cond in conditions:

        _data_stretch = resp_stretch_allcond[cond]
        _mean = _data_stretch.mean(axis=0)
        _std = _data_stretch.std(axis=0)
        (line,) = ax.plot(_mean, label=f"{cond} n:{_data_stretch.shape[0]}")
        ax.fill_between(range(_mean.size), _mean - _std, _mean + _std, color=line.get_color(), alpha=0.1)
    plt.legend()
    plt.title(f'{sujet} RESP')
    plt.rcParams.update({'font.size': 12})
    # plt.show()

    fig_mean_resp.savefig(f"stretch_mean_resp_{sujet}_STRICT.jpeg")

    fig_mean_CO2, ax = plt.subplots()
    for cond in conditions:

        _data_stretch = co2_stretch_allcond[cond]
        _mean = _data_stretch.mean(axis=0)
        _std = _data_stretch.std(axis=0)
        (line,) = ax.plot(_mean, label=f"{cond} n:{_data_stretch.shape[0]}")
        ax.fill_between(range(_mean.size), _mean - _std, _mean + _std, color=line.get_color(), alpha=0.1)
    plt.legend()
    plt.title(f'{sujet} CO2')
    plt.rcParams.update({'font.size': 12})
    # plt.show()

    fig_mean_CO2.savefig(f"stretch_mean_CO2_{sujet}_STRICT.jpeg")









#sujet = 'NS217'
def export_cond_relabeling_fig(sujet):

    ######## LOAD ########
    #### load respfeatures
    os.chdir(os.path.join(path_results, 'respi', 'respfeatures', sujet))

    cond_list = []
    respfeatures = []
    
    for cond in conditions:

        for trial_i in range(ntrail_dict_allpatient[sujet][cond]):

            _respfeatures = pd.read_excel(f"{sujet}_{cond}0{trial_i+1}_respfeatures.xlsx")
            _respfeatures['trial'] = [f"{cond}0{trial_i+1}"] * _respfeatures.shape[0]
            _respfeatures['cond'] = [cond] * _respfeatures.shape[0]
            
            cond_list.append(f"{cond}0{trial_i+1}")
            respfeatures.append(_respfeatures)

    respfeatures_allcond_raw_import = pd.concat(respfeatures)

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

    extract_respfeatures_and_mean_figures(sujet)
    export_cond_relabeling_fig(sujet)















