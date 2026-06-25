

import joblib


from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *

import joblib
import statsmodels.api as sm
from io import StringIO
from statsmodels.tools.sm_exceptions import ValueWarning
import warnings





################################
######## EXPORT DF ########
################################



#sujet = sujet_list[0]
def export_Pxx_cond_df_sujet(sujet):

    for band in freq_band_dict:

        print(band)
        
        #### load
        path_load_Pxx = os.path.join(path_precompute, 'Pxx')

        _xr_raw = xr.open_dataarray(os.path.join(path_load_Pxx, f"{sujet}_xr_Pxx_RAW.nc")).sel(band=band).rename({'phase' : 'phase_cycle'})
        _xr_relabel = xr.open_dataarray(os.path.join(path_load_Pxx, f"{sujet}_xr_Pxx_RELABEL.nc")).sel(band=band).rename({'phase' : 'phase_cycle'})

        #### generate df
        df_band_Pxx_sujet_raw = _xr_raw.to_dataframe().reset_index()
        df_band_Pxx_sujet_relabel = _xr_relabel.to_dataframe().reset_index()
        
        #### clean
        nan_exclude_vec = ~df_band_Pxx_sujet_raw['Pxx'].isna()
        df_band_Pxx_sujet_raw = df_band_Pxx_sujet_raw[nan_exclude_vec]
        df_band_Pxx_sujet_raw = df_band_Pxx_sujet_raw.query(f"ROI in {ROI_short_list}").reset_index(drop=True)

        nan_exclude_vec = ~df_band_Pxx_sujet_relabel['Pxx'].isna()
        df_band_Pxx_sujet_relabel = df_band_Pxx_sujet_relabel[nan_exclude_vec]
        df_band_Pxx_sujet_relabel = df_band_Pxx_sujet_relabel.query(f"ROI in {ROI_short_list}").reset_index(drop=True)

        #### prepare cond
        df_band_Pxx_sujet_attention = df_band_Pxx_sujet_raw.query(f"cond in ['RB', 'AoRoCo']")
        df_band_Pxx_sujet_HV = df_band_Pxx_sujet_raw.query(f"cond in ['RB', 'HV']")
        df_band_Pxx_sujet_VCCP = df_band_Pxx_sujet_relabel

        A_cond, R_cond, C_cond = [], [], []
        for row_i, row_val in df_band_Pxx_sujet_VCCP.iterrows():

            A_cond.append(row_val['cond'][1])
            R_cond.append(row_val['cond'][3])
            C_cond.append(row_val['cond'][5])

        df_band_Pxx_sujet_VCCP['A'], df_band_Pxx_sujet_VCCP['R'], df_band_Pxx_sujet_VCCP['C'] = A_cond, R_cond, C_cond
        df_band_Pxx_sujet_VCCP = df_band_Pxx_sujet_VCCP.drop(columns=['cond'])

        df_band_Pxx_sujet_attention['sujet'] = [sujet] * df_band_Pxx_sujet_attention.shape[0]
        df_band_Pxx_sujet_HV['sujet'] = [sujet] * df_band_Pxx_sujet_HV.shape[0]
        df_band_Pxx_sujet_VCCP['sujet'] = [sujet] * df_band_Pxx_sujet_VCCP.shape[0]

        #### calibrate for different protocols
        df_band_Pxx_sujet_VCCP['C'] = df_band_Pxx_sujet_VCCP['C'].replace({'c' : 'o'})
        df_band_Pxx_sujet_attention['cond'] = df_band_Pxx_sujet_attention['cond'].replace({'AoRoCo' : 'CB'})

        #### export
        path_export_df = os.path.join(path_precompute, 'EXPORT_DF', 'patient_wise', 'Pxx')

        df_band_Pxx_sujet_attention.to_excel(os.path.join(path_export_df, f"{sujet}_df_Pxx_{band}_attention_R.xlsx"))
        df_band_Pxx_sujet_HV.to_excel(os.path.join(path_export_df, f"{sujet}_df_Pxx_{band}_HV_R.xlsx"))
        df_band_Pxx_sujet_VCCP.to_excel(os.path.join(path_export_df, f"{sujet}_df_Pxx_{band}_VCCP_R.xlsx"))




#sujet = sujet_list[0]
def export_Pxx_reg_df_sujet(sujet):

    #### config
    sel_vec_rf_metric_scaled = ['sujet', 'cycle', 'A_relabel', 'R_relabel', 'C_relabel', 'cond',
       'A_scaled', 'R_scaled', 'C_scaled']

    #band = 'theta'
    for band in freq_band_dict:

        print(band)
        
        #### load
        chanlist, _, _, localist = get_chanlist_localist(sujet)
        path_load_Pxx = os.path.join(path_precompute, 'Pxx')
        path_load_rf_scaled = os.path.join(path_precompute, 'RESPI', 'respfeatures', sujet)
        
        _xr_raw = xr.open_dataarray(os.path.join(path_load_Pxx, f"{sujet}_xr_Pxx_RAW.nc")).sel(band=band).rename({'phase' : 'phase_cycle'})
        _xr_relabel = xr.open_dataarray(os.path.join(path_load_Pxx, f"{sujet}_xr_Pxx_RELABEL.nc")).sel(band=band).rename({'phase' : 'phase_cycle'})

        #### generate df
        df_band_Pxx_sujet_raw = _xr_raw.to_dataframe().reset_index()
        df_band_Pxx_sujet_relabel = _xr_relabel.to_dataframe().reset_index()
        
        #### clean
        nan_exclude_vec = ~df_band_Pxx_sujet_raw['Pxx'].isna()
        df_band_Pxx_sujet_raw = df_band_Pxx_sujet_raw[nan_exclude_vec]
        df_band_Pxx_sujet_raw = df_band_Pxx_sujet_raw.query(f"ROI in {ROI_short_list}").reset_index(drop=True)

        nan_exclude_vec = ~df_band_Pxx_sujet_relabel['Pxx'].isna()
        df_band_Pxx_sujet_relabel = df_band_Pxx_sujet_relabel[nan_exclude_vec]
        df_band_Pxx_sujet_relabel = df_band_Pxx_sujet_relabel.query(f"ROI in {ROI_short_list}").reset_index(drop=True)

        #### prepare cond
        df_band_Pxx_sujet_attention = df_band_Pxx_sujet_raw.query(f"cond in ['RB', 'AoRoCo']")
        df_band_Pxx_sujet_attention['cond'] = df_band_Pxx_sujet_attention['cond'].replace({'AoRoCo' : 'CB'})
        df_band_Pxx_sujet_VCCP = df_band_Pxx_sujet_relabel

        #### load Pxx and rf scaled
        rf_SD_scaled_VCCP = pd.read_excel(os.path.join(path_load_rf_scaled, f"{sujet}_df_scaled_VCCP.xlsx")).query(f"cond_relabel != 'excluded'")
        rf_SD_scaled_RB = pd.read_excel(os.path.join(path_load_rf_scaled, f"{sujet}_df_scaled_RB.xlsx"))

        #### index identification VCCP
        rf_SD_scaled_VCCP = rf_SD_scaled_VCCP.sort_values(['cond_relabel', 'trial']).reset_index(drop=True)

        cycle_i_vec = []
        cycle_count = 0

        for row_i, row_val in rf_SD_scaled_VCCP.iterrows():

            if row_i != 0:
                cond_row = row_val['cond_relabel']
                if cond_row == cond_prev:
                    cycle_count += 1
                else:
                    cycle_count = 0
            cycle_i_vec.append(cycle_count)
            cond_prev = row_val['cond_relabel']

        rf_SD_scaled_VCCP['cycle'] = cycle_i_vec

        verif_count_scaled = rf_SD_scaled_VCCP.groupby('cond_relabel').count().reset_index(names='cond_relabel')[['cond_relabel', 'cycle']]
        verif_count_df = df_band_Pxx_sujet_VCCP.query(f"phase_cycle == 'inspi' and ROI == '{ROI_short_list[0]}'").groupby('cond').count().reset_index(names='cond_relabel')[['cond_relabel', 'cycle']]

        if ~(verif_count_scaled == verif_count_df).values.all():
            raise ValueError(f"NOT THE SAME NUMBER OF CYCLES")

        rf_SD_scaled_VCCP = rf_SD_scaled_VCCP.drop(columns=['cond']).rename(columns={'cond_relabel' : 'cond'})
        rf_SD_scaled_VCCP = rf_SD_scaled_VCCP[sel_vec_rf_metric_scaled]

        #### prepare attention
        rf_SD_scaled_attention = rf_SD_scaled_RB.query(f"cond in ['RB', 'AoRoCo']")
        rf_SD_scaled_attention['cond'] = rf_SD_scaled_attention['cond'].replace({'AoRoCo' : 'CB'})
        rf_SD_scaled_attention = rf_SD_scaled_attention[sel_vec_rf_metric_scaled]

        #### merge
        df_reg_VCCP = pd.merge(df_band_Pxx_sujet_VCCP, rf_SD_scaled_VCCP, on=['sujet', 'cond', 'cycle'], how='left')
        df_reg_attention = pd.merge(df_band_Pxx_sujet_attention, rf_SD_scaled_attention, on=['sujet', 'cond', 'cycle'], how='left')

        #### calibrate for different protocols
        df_reg_VCCP['C_relabel'] = df_reg_VCCP['C_relabel'].replace({'c' : 'o'})

        #### melt
        col_melt = ['A_scaled', 'R_scaled', 'C_scaled']
        df_reg_VCCP = pd.melt(df_reg_VCCP, id_vars=[col for col in df_reg_VCCP.columns if col not in col_melt], value_vars=col_melt, var_name='factor', value_name='val_scaled')
        df_reg_VCCP['factor'] = df_reg_VCCP['factor'].replace({'A_scaled' : 'A', 'R_scaled' : 'R', 'C_scaled' : 'C'})

        df_reg_attention = pd.melt(df_reg_attention, id_vars=[col for col in df_reg_attention.columns if col not in col_melt], value_vars=col_melt, var_name='factor', value_name='val_scaled')
        df_reg_attention['factor'] = df_reg_attention['factor'].replace({'A_scaled' : 'A', 'R_scaled' : 'R', 'C_scaled' : 'C'})

        #### export
        path_export_df = os.path.join(path_precompute, 'EXPORT_DF', 'patient_wise', 'Pxx')

        df_reg_attention.to_excel(os.path.join(path_export_df, f"{sujet}_df_reg_Pxx_{band}_attention_R.xlsx"))
        df_reg_VCCP.to_excel(os.path.join(path_export_df, f"{sujet}_df_reg_Pxx_{band}_VCCP_R.xlsx"))




#sujet = sujet_list[0]
def export_Cxy_cond_df_sujet(sujet):
        
    #### load
    path_load_Cxy = os.path.join(path_precompute, 'Cxy')

    chanlist, _, _, localist = get_chanlist_localist(sujet)
    _xr_Cxy = xr.open_dataarray(os.path.join(path_load_Cxy, f"{sujet}_Cxy.nc"))

    ### extract Cxy around respi
    respfeature = get_respfeatures(sujet)

    trial_list = get_alltrials_sujet(sujet)

    resp_med_allcond = {}

    #cond = conditions[-1]
    for cond in conditions:

        trial_list_cond = [trial for trial in trial_list if trial.find(cond) != -1]

        _resp_med_cond = []

        for trial_i, trial in enumerate(trial_list_cond):

            _respfeature_trial = respfeature.query(f"cond == '{cond}' and trial == {trial_i+1}")
            _resp_med_cond.append(_respfeature_trial['cycle_freq'].median())

        resp_med_allcond[cond] = np.median(_resp_med_cond)

    data_Cxy_med = []
    hzCxy = _xr_Cxy['freq'].values

    for cond in conditions:

        sel_vec_hzCxy = (hzCxy < (resp_med_allcond[cond] + Cxy_extraction_range)) & ((resp_med_allcond[cond] - Cxy_extraction_range) < hzCxy) 
        data_Cxy_med.append(_xr_Cxy.sel(cond=cond, freq=sel_vec_hzCxy).median('freq').values)

    data_Cxy_med = np.stack(data_Cxy_med)

    xr_coords = {'cond' : conditions, 'chan' : chanlist, 'ROI' : ('chan', localist['loca'].values)}
    _xr_Cxy_med = xr.DataArray(data_Cxy_med, dims=['cond', 'chan'], coords=xr_coords) 

    #### generate df
    df_Cxy_sujet_raw = _xr_Cxy_med.to_dataframe(name='Cxy').reset_index()
    
    #### clean
    df_Cxy_sujet = df_Cxy_sujet_raw.query(f"ROI in {ROI_short_list}").reset_index(drop=True)

    #### prepare cond
    df_Cxy_sujet_attention = df_Cxy_sujet.query(f"cond in ['RB', 'AoRoCo']")
    df_Cxy_sujet_HV = df_Cxy_sujet.query(f"cond in ['RB', 'HV']")
    df_Cxy_sujet_protocol = df_Cxy_sujet.query(f"cond not in ['RB', 'HV']")

    A_cond, R_cond, C_cond = [], [], []
    for row_i, row_val in df_Cxy_sujet_protocol.iterrows():

        A_cond.append(row_val['cond'][1])
        R_cond.append(row_val['cond'][3])
        C_cond.append(row_val['cond'][5])

    df_Cxy_sujet_protocol['A'], df_Cxy_sujet_protocol['R'], df_Cxy_sujet_protocol['C'] = A_cond, R_cond, C_cond
    df_Cxy_sujet_protocol = df_Cxy_sujet_protocol.drop(columns=['cond'])

    df_Cxy_sujet_attention['sujet'] = [sujet] * df_Cxy_sujet_attention.shape[0]
    df_Cxy_sujet_HV['sujet'] = [sujet] * df_Cxy_sujet_HV.shape[0]
    df_Cxy_sujet_protocol['sujet'] = [sujet] * df_Cxy_sujet_protocol.shape[0]

    #### calibrate for different protocols
    df_Cxy_sujet_protocol['C'] = df_Cxy_sujet_protocol['C'].replace({'c' : 'o'})
    df_Cxy_sujet_attention['cond'] = df_Cxy_sujet_attention['cond'].replace({'AoRoCo' : 'CB'})

    #### export
    path_export_df = os.path.join(path_precompute, 'EXPORT_DF', 'patient_wise', 'Cxy')

    df_Cxy_sujet_attention.to_excel(os.path.join(path_export_df, f"{sujet}_df_Cxy_attention_R.xlsx"))
    df_Cxy_sujet_HV.to_excel(os.path.join(path_export_df, f"{sujet}_df_Cxy_HV_R.xlsx"))
    df_Cxy_sujet_protocol.to_excel(os.path.join(path_export_df, f"{sujet}_df_Cxy_protocol_R.xlsx"))





def export_Pxx_df():
        
    for band in freq_band_dict:

        print(band)
        
        #### load
        os.chdir(os.path.join(path_precompute, 'TF', 'session'))

        df_band_Pxx_allsujet = []

        for sujet in sujet_list:

            _xr = xr.open_dataarray(f"{sujet}_xr_Pxx.nc").sel(pre_post='post', band=band).rename({'phase' : 'phase_cycle'})
            _xr_whole = _xr.median('phase_cycle').expand_dims(phase_cycle=['whole'])
            _xr_allphase = xr.concat([_xr, _xr_whole], dim='phase_cycle')

            df_band_Pxx_allsujet.append(_xr_allphase.to_dataframe().reset_index().dropna())

        df_band_Pxx_allsujet = pd.concat(df_band_Pxx_allsujet)

        df_band_Pxx_allsujet[['resp', 'state']] = df_band_Pxx_allsujet['cond'].str.split('_', expand=True).rename(columns={0: 'resp', 1: 'state'})
        df_band_Pxx_allsujet = df_band_Pxx_allsujet.drop(columns=['cond'])

        df_band_Pxx_allsujet = df_band_Pxx_allsujet.query(f"ROI in {ROI_short_list}")

        #### export
        os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_R'))
        df_band_Pxx_allsujet.to_excel(f"df_R_{band}.xlsx")

    








################################
######## PERM ########
################################



def spearman_permutation_test(X, Y, n_perm_spearman, percentile_thresh_list):

    rho_obs, _ = scipy.stats.spearmanr(X, Y)

    #### ensure that there is no double
    Y_raw_i = np.arange(Y.size)
    seen_perm = set()     
    rand_perm_i = []       

    for i in range(n_perm_spearman):

        _Y_perm_i = np.random.permutation(Y_raw_i)
        _perm_key = tuple(_Y_perm_i)

        if _perm_key not in seen_perm:
            seen_perm.add(_perm_key)
            rand_perm_i.append(_Y_perm_i)

    #### compute perm
    null_dist = []
    for _perm_i in rand_perm_i:

        _rho, _ = scipy.stats.spearmanr(X, Y[_perm_i])
        null_dist.append(_rho)

    null_dist = np.array(null_dist)

    lower, upper = np.percentile(null_dist, percentile_thresh_list)
    signi = rho_obs < lower or rho_obs > upper

    return rho_obs, signi













################################
######## REG ########
################################





def get_reg_allsujet():

    if os.path.exists(os.path.join(path_precompute, 'TF', 'session', 'df_reg', f"df_reg.xlsx")):
        
        print('ALREADY COMPUTED')
    
    #### load Pxx
    das = []
    chan_list_shape = []

    for sujet in sujet_list:
        
        fp = os.path.join(path_precompute, "TF", "session", f"{sujet}_xr_Pxx.nc")
        da = xr.load_dataarray(fp)

        da = da.expand_dims(sujet=[sujet])

        chan_list_shape.append(da.shape[2])
        das.append(da)

    max_nchan = np.max(np.array(chan_list_shape))

    das_padded = []
    for da in das:

        nchan = da.sizes["chan"]
        
        _da_padded = da.assign_coords(chan=np.arange(nchan)).reindex(chan=np.arange(max_nchan)).assign_coords(
            chan_label=("chan", list(da["chan"].values) + [np.nan] * (max_nchan - nchan)),
            ROI=("chan", list(da["ROI"].values) + [np.nan] * (max_nchan - nchan)))

        das_padded.append(_da_padded)

    da_Pxx_allsujet_phase_cycle = xr.concat(das_padded, dim="sujet")
    da_Pxx_allsujet_phase_cycle = da_Pxx_allsujet_phase_cycle.rename({'phase' : 'phase_cycle'})
    da_whole = da_Pxx_allsujet_phase_cycle.median('phase_cycle').expand_dims(phase_cycle=["whole"])
    da_Pxx_allphase_cycle = xr.concat([da_Pxx_allsujet_phase_cycle, da_whole],dim="phase_cycle")
    
    #### load resp features
    rf_metrics_raw = ['cycle_duration', 'inspi_duration', 'expi_duration', 'cycle_freq', 'inspi_volume',
       'expi_volume', 'total_amplitude', 'inspi_amplitude', 'expi_amplitude',
       'total_volume']
    phase_cycle_list = ['whole', 'inspi', 'expi']
    mapping_pase_cycle_list = {'whole' : ['cycle_duration', 'cycle_freq', 'total_amplitude', 'total_volume', 'oc_ratio', 'oc_val'], 
                               'inspi' : ['inspi_duration', 'inspi_cycle_freq', 'inspi_amplitude', 'inspi_volume', 'oc_ratio', 'oc_val'], 
                               'expi' : ['expi_duration', 'expi_cycle_freq', 'expi_amplitude', 'expi_volume', 'oc_ratio', 'oc_val']}
    rf_metric_allphase_cycle = ['duration', 'cycle_freq', 'amplitude', 'volume', 'oc_ratio', 'oc_val']
    
    resp_vars = {}

    for sujet in sujet_list:

        # load respiration features
        os.chdir(os.path.join(path_precompute, 'RESP', 'respfeatures')) 
        _oc_respfeatures = pd.read_excel(f'{sujet}_df_oc.xlsx')

        _respfeatures = pd.read_excel(f'{sujet}_respfeatures_cleaned_label.xlsx')
        _respfeatures_pre = pd.read_excel(f'{sujet}_respfeatures_cleaned_label_pre.xlsx')

        _cond_da = {}
        
        for cond_i, cond in enumerate(conditions):

            resp_pre  = _respfeatures_pre.query(f"cond == '{cond}'").drop(columns=['Unnamed: 0'])[rf_metrics_raw]
            resp_post = _respfeatures.query(f"cond == '{cond}'").drop(columns=['Unnamed: 0'])[rf_metrics_raw]

            resp_pre[['oc_ratio', 'oc_val']] = _oc_respfeatures.query(f"cond == '{cond}'")[['oc_ratio', 'oc_val']]
            resp_post[['oc_ratio', 'oc_val']] = _oc_respfeatures.query(f"cond == '{cond}'")[['oc_ratio', 'oc_val']]

            resp_pre['inspi_cycle_freq'], resp_pre['expi_cycle_freq'] = 1/resp_pre['inspi_duration'], 1/resp_pre['expi_duration']
            resp_post['inspi_cycle_freq'], resp_post['expi_cycle_freq'] = 1/resp_post['inspi_duration'], 1/resp_post['expi_duration']
            
            allmetric_resp_pre = []
            allmetric_resp_post = []

            for _phase_cycle in phase_cycle_list:

                allmetric_resp_pre.append(resp_pre[mapping_pase_cycle_list[_phase_cycle]].to_numpy())
                allmetric_resp_post.append(resp_post[mapping_pase_cycle_list[_phase_cycle]].to_numpy())

            # convert to numpy
            pre_vals  = np.stack(allmetric_resp_pre)
            post_vals = np.stack(allmetric_resp_post)

            ds_resp = xr.Dataset(
                data_vars={
                    "respfeatures": (("phase_group", "phase_cycle", "cycle", "feature"),
                                    np.stack([pre_vals, post_vals], axis=0))
                },
                coords={
                    "phase_group": ["pre", "post"],
                    "phase_cycle": phase_cycle_list,
                    "cycle": np.arange(pre_vals.shape[1]),
                    "feature": rf_metric_allphase_cycle,
                    "sujet": sujet,
                    "cond": cond,
                }
            ).expand_dims(("sujet", "cond"))

            # collect
            _cond_da[cond] = ds_resp

        resp_vars[sujet] = xr.concat(list(_cond_da.values()), dim="cond").sortby("cond")

    ds_respfeatures_all = xr.concat(list(resp_vars.values()), dim="sujet").sortby("sujet")
    ds_respfeatures_all = ds_respfeatures_all.rename({"phase_group": "pre_post"})

    #### merge
    ds_whole = xr.merge([da_Pxx_allphase_cycle, ds_respfeatures_all])

    #### inspect
    if debug:

        time_phase_list = ['pre', 'post']

        for sujet in sujet_list:

            chanlist_sujet = get_chanlist(sujet)[0]

            for chan_i, chan in enumerate(chanlist_sujet):

                for cond in conditions:

                    for band in freq_band_dict:

                        for phase_protocol in time_phase_list:

                            for cycle_phase in phase_cycle_list:

                                for feature in rf_metric_allphase_cycle:

                                    _ds_Pxx = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='pre', phase_cycle=cycle_phase, feature=feature)['Pxx']
                                    X_pre = _ds_Pxx.values[~np.isnan(_ds_Pxx.values)]

                                    _ds_Pxx = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature=feature)['Pxx']
                                    X_post = _ds_Pxx.values[~np.isnan(_ds_Pxx.values)]

                                    _ds_Pxx = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='pre', phase_cycle=cycle_phase, feature='oc_ratio')['Pxx']
                                    X_pre = _ds_Pxx.values[~np.isnan(_ds_Pxx.values)]

                                    _ds_Pxx = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature='oc_ratio')['Pxx']
                                    X_post = _ds_Pxx.values[~np.isnan(_ds_Pxx.values)]

                                    _df_rf = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='pre', phase_cycle=cycle_phase, feature='oc_ratio')['respfeatures']
                                    Y_pre = _df_rf.values[~np.isnan(_df_rf.values)]

                                    _df_rf = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature='oc_ratio')['respfeatures']
                                    Y_post = _df_rf.values[~np.isnan(_df_rf.values)]

                                    _df_rf = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='pre', phase_cycle=cycle_phase, feature=feature)['respfeatures']
                                    Y_pre = _df_rf.values[~np.isnan(_df_rf.values)]

                                    _df_rf = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature=feature)['respfeatures']
                                    Y_post = _df_rf.values[~np.isnan(_df_rf.values)]

    #### compute OLS
    percentile_thresh_list = [2.5, 97.5]
    n_perm_spearman = 500
    time_phase_list = ['pre', 'post']

    df_loca_allsujet = get_df_loca_allsujet()

    reg_allsujet_data = []

    for sujet in sujet_list:

        print(sujet)

        # chan_i, chan = 0, chanlist_sujet[0]
        def generate_reg_patient(sujet, chan_i, chan):
        # for chan_i, chan in enumerate(chanlist_sujet):

            chanlist_sujet = get_chanlist(sujet)[0]
            print_advancement(chan_i, len(chanlist_sujet), [25,50,75])

            reg_chan = []

            for cond in conditions:

                for band in freq_band_dict:

                    for phase_protocol in time_phase_list:

                        for cycle_phase in phase_cycle_list:

                            for feature in rf_metric_allphase_cycle:

                                _ds_Pxx = ds_whole.sel(sujet=sujet, chan=chan_i, cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature=feature)['Pxx']
                                X = _ds_Pxx.values[~np.isnan(_ds_Pxx.values)]

                                _df_rf = ds_whole.sel(sujet=sujet, chan=chan_i, cond=cond, band=band, pre_post=phase_protocol, phase_cycle=cycle_phase, feature=feature)['respfeatures']
                                Y = _df_rf.values[~np.isnan(_df_rf.values)]

                                _rho, _signi = spearman_permutation_test(X, Y, n_perm_spearman, percentile_thresh_list)

                                X = sm.add_constant(X)

                                _mdl = sm.OLS(Y, X).fit()

                                with warnings.catch_warnings():
                                    warnings.simplefilter("ignore", category=UserWarning)
                                    warnings.simplefilter("ignore", category=ValueWarning)

                                    table = pd.read_html(
                                        StringIO(_mdl.summary().tables[1].as_html()),
                                        header=0,
                                        index_col=0
                                    )[0]

                                coeff     = table["coef"].values[-1]
                                pval      = table["P>|t|"].values[-1]

                                _ROI = _ds_Pxx['ROI'].data.tolist()

                                _df = pd.DataFrame({'sujet' : [sujet], 'chan' : [chan], 'ROI' : [_ROI], 'cond' : [cond], 'band' : [band], 'phase_protocol' : [phase_protocol], 'rf_metric' : [feature],
                                                    'phase_cycle' : [cycle_phase], 'rho' : [_rho], 'rho_signi' : [_signi*1], 'OLS_a' : coeff, 'OLS_p' : pval})

                                reg_chan.append(_df)

            return pd.concat(reg_chan)

        chanlist_sujet = get_chanlist(sujet)[0]
        reg_patient_chan = joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(generate_reg_patient)(sujet, chan_i, chan) for chan_i, chan in enumerate(chanlist_sujet))

        df_reg_patient = pd.concat(reg_patient_chan)

        df_reg_patient = df_reg_patient.query(f"ROI in {ROI_short_list}")

        os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_reg'))
        df_reg_patient.to_excel(f"df_reg_{sujet}.xlsx")

        #### inspect
        if debug:

            _df = generate_reg_patient(sujet, chan_i, chan)
            _df.query(f"cond=='{cond}' and band=='{band}' and phase_protocol=='pre' and phase_cycle=='{cycle_phase}' and rf_metric=='oc_ratio'")
            _df.query(f"cond=='{cond}' and band=='{band}' and phase_protocol=='post' and phase_cycle=='{cycle_phase}' and rf_metric=='oc_ratio'")

            _df.query(f"cond=='{cond}' and band=='{band}' and phase_protocol=='pre' and phase_cycle=='{cycle_phase}' and rf_metric=='duration'")
            _df.query(f"cond=='{cond}' and band=='{band}' and phase_protocol=='post' and phase_cycle=='{cycle_phase}' and rf_metric=='duration'")


        #### add df
        reg_allsujet_data.append(df_reg_patient)

    reg_allsujet_data = pd.concat(reg_allsujet_data)

    os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_reg'))
    reg_allsujet_data.to_excel(f"df_reg.xlsx")

    reg_allsujet_data_R = reg_allsujet_data.copy()
    reg_allsujet_data_R[['resp', 'state']] = reg_allsujet_data['cond'].str.split('_', expand=True).rename(columns={0: 'resp', 1: 'state'})
    reg_allsujet_data_R = reg_allsujet_data_R.drop(columns=['cond'])

    os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_R'))
    reg_allsujet_data_R.to_excel(f"df_reg_R.xlsx")








def get_reg_allROI():

    if os.path.exists(os.path.join(path_precompute, 'TF', 'session', 'df_reg', f"df_reg_allROI.xlsx")):
        
        print('ALREADY COMPUTED')
    
    #### load Pxx
    das = []
    chan_list_shape = []

    for sujet in sujet_list:
        
        fp = os.path.join(path_precompute, "TF", "session", f"{sujet}_xr_Pxx.nc")
        da = xr.load_dataarray(fp)

        da = da.expand_dims(sujet=[sujet])

        chan_list_shape.append(da.shape[2])
        das.append(da)

    max_nchan = np.max(np.array(chan_list_shape))

    das_padded = []
    for da in das:

        nchan = da.sizes["chan"]
        
        _da_padded = da.assign_coords(chan=np.arange(nchan)).reindex(chan=np.arange(max_nchan)).assign_coords(
            chan_label=("chan", list(da["chan"].values) + [np.nan] * (max_nchan - nchan)),
            ROI=("chan", list(da["ROI"].values) + [np.nan] * (max_nchan - nchan)))

        das_padded.append(_da_padded)

    da_Pxx_allsujet_phase_cycle = xr.concat(das_padded, dim="sujet")
    da_Pxx_allsujet_phase_cycle = da_Pxx_allsujet_phase_cycle.rename({'phase' : 'phase_cycle'})
    da_whole = da_Pxx_allsujet_phase_cycle.median('phase_cycle').expand_dims(phase_cycle=["whole"])
    da_Pxx_allphase_cycle = xr.concat([da_Pxx_allsujet_phase_cycle, da_whole],dim="phase_cycle")
    
    #### load resp features
    rf_metrics_raw = ['cycle_duration', 'inspi_duration', 'expi_duration', 'cycle_freq', 'inspi_volume',
       'expi_volume', 'total_amplitude', 'inspi_amplitude', 'expi_amplitude',
       'total_volume']
    phase_cycle_list = ['whole', 'inspi', 'expi']
    mapping_pase_cycle_list = {'whole' : ['cycle_duration', 'cycle_freq', 'total_amplitude', 'total_volume', 'oc_ratio', 'oc_val'], 
                               'inspi' : ['inspi_duration', 'inspi_cycle_freq', 'inspi_amplitude', 'inspi_volume', 'oc_ratio', 'oc_val'], 
                               'expi' : ['expi_duration', 'expi_cycle_freq', 'expi_amplitude', 'expi_volume', 'oc_ratio', 'oc_val']}
    rf_metric_allphase_cycle = ['duration', 'cycle_freq', 'amplitude', 'volume', 'oc_ratio', 'oc_val']
    
    resp_vars = {}

    for sujet in sujet_list:

        # load respiration features
        os.chdir(os.path.join(path_precompute, 'RESP', 'respfeatures')) 
        _oc_respfeatures = pd.read_excel(f'{sujet}_df_oc.xlsx')

        _respfeatures = pd.read_excel(f'{sujet}_respfeatures_cleaned_label.xlsx')
        _respfeatures_pre = pd.read_excel(f'{sujet}_respfeatures_cleaned_label_pre.xlsx')

        _cond_da = {}
        
        for cond_i, cond in enumerate(conditions):

            resp_pre  = _respfeatures_pre.query(f"cond == '{cond}'").drop(columns=['Unnamed: 0'])[rf_metrics_raw]
            resp_post = _respfeatures.query(f"cond == '{cond}'").drop(columns=['Unnamed: 0'])[rf_metrics_raw]

            resp_pre[['oc_ratio', 'oc_val']] = _oc_respfeatures.query(f"cond == '{cond}'")[['oc_ratio', 'oc_val']]
            resp_post[['oc_ratio', 'oc_val']] = _oc_respfeatures.query(f"cond == '{cond}'")[['oc_ratio', 'oc_val']]

            resp_pre['inspi_cycle_freq'], resp_pre['expi_cycle_freq'] = 1/resp_pre['inspi_duration'], 1/resp_pre['expi_duration']
            resp_post['inspi_cycle_freq'], resp_post['expi_cycle_freq'] = 1/resp_post['inspi_duration'], 1/resp_post['expi_duration']
            
            allmetric_resp_pre = []
            allmetric_resp_post = []

            for _phase_cycle in phase_cycle_list:

                allmetric_resp_pre.append(resp_pre[mapping_pase_cycle_list[_phase_cycle]].to_numpy())
                allmetric_resp_post.append(resp_post[mapping_pase_cycle_list[_phase_cycle]].to_numpy())

            # convert to numpy
            pre_vals  = np.stack(allmetric_resp_pre)
            post_vals = np.stack(allmetric_resp_post)

            ds_resp = xr.Dataset(
                data_vars={
                    "respfeatures": (("phase_group", "phase_cycle", "cycle", "feature"),
                                    np.stack([pre_vals, post_vals], axis=0))
                },
                coords={
                    "phase_group": ["pre", "post"],
                    "phase_cycle": phase_cycle_list,
                    "cycle": np.arange(pre_vals.shape[1]),
                    "feature": rf_metric_allphase_cycle,
                    "sujet": sujet,
                    "cond": cond,
                }
            ).expand_dims(("sujet", "cond"))

            # collect
            _cond_da[cond] = ds_resp

        resp_vars[sujet] = xr.concat(list(_cond_da.values()), dim="cond").sortby("cond")

    ds_respfeatures_all = xr.concat(list(resp_vars.values()), dim="sujet").sortby("sujet")
    ds_respfeatures_all = ds_respfeatures_all.rename({"phase_group": "pre_post"})

    #### merge
    ds_whole = xr.merge([da_Pxx_allphase_cycle, ds_respfeatures_all])

    #### inspect
    if debug:

        time_phase_list = ['pre', 'post']

        for sujet in sujet_list:

            chanlist_sujet = get_chanlist(sujet)[0]

            for chan_i, chan in enumerate(chanlist_sujet):

                for cond in conditions:

                    for band in freq_band_dict:

                        for phase_protocol in time_phase_list:

                            for cycle_phase in phase_cycle_list:

                                for feature in rf_metric_allphase_cycle:

                                    _ds_Pxx = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='pre', phase_cycle=cycle_phase, feature=feature)['Pxx']
                                    X_pre = _ds_Pxx.values[~np.isnan(_ds_Pxx.values)]

                                    _ds_Pxx = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature=feature)['Pxx']
                                    X_post = _ds_Pxx.values[~np.isnan(_ds_Pxx.values)]

                                    _ds_Pxx = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='pre', phase_cycle=cycle_phase, feature='oc_ratio')['Pxx']
                                    X_pre = _ds_Pxx.values[~np.isnan(_ds_Pxx.values)]

                                    _ds_Pxx = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature='oc_ratio')['Pxx']
                                    X_post = _ds_Pxx.values[~np.isnan(_ds_Pxx.values)]

                                    _df_rf = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='pre', phase_cycle=cycle_phase, feature='oc_ratio')['respfeatures']
                                    Y_pre = _df_rf.values[~np.isnan(_df_rf.values)]

                                    _df_rf = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature='oc_ratio')['respfeatures']
                                    Y_post = _df_rf.values[~np.isnan(_df_rf.values)]

                                    _df_rf = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='pre', phase_cycle=cycle_phase, feature=feature)['respfeatures']
                                    Y_pre = _df_rf.values[~np.isnan(_df_rf.values)]

                                    _df_rf = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature=feature)['respfeatures']
                                    Y_post = _df_rf.values[~np.isnan(_df_rf.values)]

    #### identify ROI extarction patient
    df_loca_allsujet = get_df_loca_allsujet()
    ROI_list = []
    for _ROI in df_loca_allsujet['loca'].unique():
        
        if isinstance(_ROI, str):
        
            if _ROI.find('UNSORTED') == -1:

                ROI_list.append(_ROI)

    ROI_id_extraction = {}

    for ROI in ROI_short_list:

        ROI_id_extraction[ROI] = {'sujet_list' : [], 'chan_i_list' : {}}

        _df_ROI = df_loca_allsujet.query(f"loca == '{ROI}'")

        for sujet in _df_ROI['sujet'].unique():

            ROI_id_extraction[ROI]['sujet_list'].append(sujet)
            _chan_list_sel = _df_ROI.query(f"sujet == '{sujet}'")['chan'].values

            _chanlist_sujet_all, _ = get_chanlist(sujet)
            _chan_i_list = [np.where(_chanlist_sujet_all == _chan)[0][0] for _chan in _chan_list_sel]
            ROI_id_extraction[ROI]['chan_i_list'][sujet] = _chan_i_list


    #### example OLS
    percentile_thresh_list = [2.5, 97.5]
    n_perm_spearman = 500
    time_phase_list = ['pre', 'post']

    #ROI_i, ROI = 0, ROI_list[0]
    for ROI_i, ROI in enumerate(ROI_short_list):

        if ROI not in ROI_short_list:
            continue

        cond = 'rsp_ctrl'
        feature = 'oc_ratio'
        phase_protocol = 'post'

        for band in ['theta', 'gamma']:

            for cycle_phase in ['inspi', 'expi']:

                X = []
                Y = []
                _ds_Pxx = ds_whole.sel(cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature=feature)
                _df_rf = ds_whole.sel(cond=cond, band=band, pre_post=phase_protocol, phase_cycle=cycle_phase, feature=feature)

                for sujet in ROI_id_extraction[ROI]['sujet_list']:

                    _ds_ROI_extraction_Pxx = _ds_Pxx.sel(sujet=sujet, chan=ROI_id_extraction[ROI]['chan_i_list'][sujet])['Pxx'].values
                    valid_cols = ~np.all(np.isnan(_ds_ROI_extraction_Pxx), axis=0)
                    _data_X = np.median(_ds_ROI_extraction_Pxx[:,valid_cols], axis=0)
                    X.append(_data_X.reshape(-1))

                    _ds_ROI_extraction_rf = _df_rf.sel(sujet=sujet, chan=ROI_id_extraction[ROI]['chan_i_list'][sujet])['respfeatures'].values
                    valid_cols = ~np.isnan(_ds_ROI_extraction_rf)
                    _data_Y = _ds_ROI_extraction_rf[valid_cols]
                    Y.append(_data_Y.reshape(-1))

                for sujet_i, sujet in enumerate(ROI_id_extraction[ROI]['sujet_list']):

                    plt.scatter(Y[sujet_i],X[sujet_i], label=sujet)

                plt.xlabel(feature)
                plt.ylabel('Pxx')
                plt.title(f"{ROI} {band} {cycle_phase}")
                plt.legend()
                # plt.show()

                filename = os.path.join(path_results, 'Pxx', 'reg_with_RF', 'ctrl', 'examples', f"{ROI}_{band}_{cycle_phase}.png")
                plt.savefig(filename)

                plt.close('all')

                
    
    #### compute OLS
    # ROI = ROI_list[0]
    def generate_reg_ROI(ROI_i, ROI):
    # for chan_i, chan in enumerate(chanlist_sujet):

        print_advancement(ROI_i, len(ROI_list), [25,50,75])

        reg_ROI = []

        for cond in conditions:

            for band in freq_band_dict:

                for phase_protocol in time_phase_list:

                    for cycle_phase in phase_cycle_list:

                        for feature in rf_metric_allphase_cycle:

                            X = []
                            Y = []
                            _ds_Pxx = ds_whole.sel(cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature=feature)
                            _df_rf = ds_whole.sel(cond=cond, band=band, pre_post=phase_protocol, phase_cycle=cycle_phase, feature=feature)

                            for sujet in ROI_id_extraction[ROI]['sujet_list']:

                                _ds_ROI_extraction_Pxx = _ds_Pxx.sel(sujet=sujet, chan=ROI_id_extraction[ROI]['chan_i_list'][sujet])['Pxx'].values
                                valid_cols = ~np.all(np.isnan(_ds_ROI_extraction_Pxx), axis=0)
                                _data_X = np.median(_ds_ROI_extraction_Pxx[:,valid_cols], axis=0)
                                X.append(_data_X.reshape(-1))

                                _ds_ROI_extraction_rf = _df_rf.sel(sujet=sujet, chan=ROI_id_extraction[ROI]['chan_i_list'][sujet])['respfeatures'].values
                                valid_cols = ~np.isnan(_ds_ROI_extraction_rf)
                                _data_Y = _ds_ROI_extraction_rf[valid_cols]
                                Y.append(_data_Y.reshape(-1))

                            if debug:

                                for sujet_i, sujet in enumerate(ROI_id_extraction[ROI]['sujet_list']):

                                    plt.scatter(Y[sujet_i],X[sujet_i], label=sujet)

                                plt.xlabel(feature)
                                plt.ylabel('Pxx')
                                plt.title(ROI)
                                plt.legend()
                                plt.show()

                            #### make data linear
                            X_linear = np.concatenate(X)
                            Y_linear = np.concatenate(Y)

                            if debug:

                                plt.scatter(Y_linear,X_linear)
                                plt.xlabel(feature)
                                plt.ylabel('Pxx')
                                plt.title(ROI)
                                plt.legend()
                                plt.show()


                            _rho, _signi = spearman_permutation_test(X_linear, Y_linear, n_perm_spearman, percentile_thresh_list)

                            X_linear = sm.add_constant(X_linear)

                            _mdl = sm.OLS(Y_linear, X_linear).fit()

                            with warnings.catch_warnings():
                                warnings.simplefilter("ignore", category=UserWarning)
                                warnings.simplefilter("ignore", category=ValueWarning)

                                table = pd.read_html(
                                    StringIO(_mdl.summary().tables[1].as_html()),
                                    header=0,
                                    index_col=0
                                )[0]

                            coeff     = table["coef"].values[-1]
                            pval      = table["P>|t|"].values[-1]

                            _df = pd.DataFrame({'ROI' : [ROI], 'cond' : [cond], 'band' : [band], 'phase_protocol' : [phase_protocol], 'rf_metric' : [feature],
                                                'phase_cycle' : [cycle_phase], 'rho' : [_rho], 'rho_signi' : [_signi*1], 'OLS_a' : coeff, 'OLS_p' : pval})

                            reg_ROI.append(_df)

        return pd.concat(reg_ROI)

    reg_allROI_data = joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(generate_reg_ROI)(ROI_i, ROI) for ROI_i, ROI in enumerate(ROI_short_list))
    df_reg_allROI = pd.concat(reg_allROI_data)

    os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_reg'))
    df_reg_allROI.to_excel(f"df_reg_ALLROI.xlsx")

    #### generate df_R_allROI
    rf_metric_allphase_cycle_short = ['amplitude', 'oc_ratio', 'oc_val']

    def generate_reg_ROI_data_R(ROI_i, ROI):
    # for chan_i, chan in enumerate(chanlist_sujet):

        print_advancement(ROI_i, len(ROI_list), [25,50,75])

        reg_ROI = []

        for cond in conditions:

            for band in ['theta', 'beta', 'gamma']:

                for cycle_phase in ['inspi', 'expi']:

                    for feature in rf_metric_allphase_cycle_short:

                        _ds_Pxx = ds_whole.sel(cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature=feature)
                        _df_rf = ds_whole.sel(cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature=feature)

                        for sujet in ROI_id_extraction[ROI]['sujet_list']:

                            _ds_ROI_extraction_Pxx = _ds_Pxx.sel(sujet=sujet, chan=ROI_id_extraction[ROI]['chan_i_list'][sujet])['Pxx'].values
                            valid_cols = ~np.all(np.isnan(_ds_ROI_extraction_Pxx), axis=0)
                            _data_X = np.median(_ds_ROI_extraction_Pxx[:,valid_cols], axis=0)

                            _ds_ROI_extraction_rf = _df_rf.sel(sujet=sujet, chan=ROI_id_extraction[ROI]['chan_i_list'][sujet])['respfeatures'].values
                            valid_cols = ~np.isnan(_ds_ROI_extraction_rf)
                            _data_Y = _ds_ROI_extraction_rf[valid_cols]

                            _nrows = _data_X.shape[0]

                            _df = pd.DataFrame({'sujet' : [sujet] * _nrows, 'ROI' : [ROI] * _nrows, 'cond' : [cond] * _nrows, 'band' : [band] * _nrows, 
                                                'rf_metric' : [feature] * _nrows, 'phase_cycle' : [cycle_phase] * _nrows, 'cycle_i' : range(_nrows), 'rf_metric_val' : _data_Y,
                                                'Pxx' : _data_X})

                            reg_ROI.append(_df)

        return pd.concat(reg_ROI)

    reg_allROI_alldata = joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(generate_reg_ROI_data_R)(ROI_i, ROI) for ROI_i, ROI in enumerate(ROI_short_list))
    df_reg_allROI_alldata_R = pd.concat(reg_allROI_alldata)

    filename = os.path.join(path_precompute, 'TF', 'session', 'df_R', f"df_reg_ALLROI_ALLDATA_R.xlsx")
    df_reg_allROI_alldata_R.to_excel(filename)
    





################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    for sujet in sujet_list:

        export_Pxx_cond_df_sujet(sujet)
        export_Cxy_cond_df_sujet(sujet)





    export_Pxx_df()
    get_reg_allsujet()

    get_reg_allROI()



                        