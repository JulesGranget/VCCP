

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import mne
import pandas as pd

from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *

debug = False



########################################
######## TO DELETE LATER ########
########################################

#resp = respi
def debugged_detect_respiration_cycles(resp, srate, baseline_mode='manual', baseline=None, 
                              epsilon_factor1=10, epsilon_factor2=5, inspiration_adjust_on_derivative=False):
    """
    Detect respiration cycles based on:
      * crossing zeros (or crossing baseline)
      * some cleanning with euristicts

    Parameters
    ----------
    resp: np.array
        Preprocess traces of respiratory signal.
    srate: float
        Sampling rate
    baseline_mode: 'manual' / 'zero' / 'median' / 'mode'
        How to compute the baseline for zero crossings.
    baseline: float or None
        External baseline when baseline_mode='manual'
    inspration_ajust_on_derivative: bool (default False)
        For the inspiration detection, the zero crossing can be refined to auto detect the inflection point.
        This can be usefull when expiration ends with a long plateau.
    Returns
    -------
    cycles: np.array
        Indices of inspiration and expiration. shape=(num_cycle, 3)
        with [index_inspi, index_expi, index_next_inspi]
    """

    baseline = physio.get_respiration_baseline(resp, srate, baseline_mode=baseline_mode, baseline=baseline)

    #~ q90 = np.quantile(resp, 0.90)
    q10 = np.quantile(resp, 0.10)
    epsilon = (baseline - q10) / 100.

    baseline_dw = baseline - epsilon * epsilon_factor1
    baseline_insp = baseline - epsilon * epsilon_factor2

    resp0 = resp[:-1]
    resp1 = resp[1:]

    ind_dw, = np.nonzero((resp0 >= baseline_dw) & (resp1 < baseline_dw))
    
    ind_insp, = np.nonzero((resp0 >= baseline_insp) & (resp1 < baseline_insp))
    ind_insp_no_clean = ind_insp.copy()
    keep_inds = np.searchsorted(ind_insp, ind_dw, side='left')
    keep_inds = keep_inds[keep_inds > 0]
    ind_insp = ind_insp[keep_inds - 1]
    ind_insp = np.unique(ind_insp)

    ind_exp, = np.nonzero((resp0 < baseline) & (resp1 >= baseline))
    keep_inds = np.searchsorted(ind_exp, ind_insp, side='right')
    keep_inds = keep_inds[keep_inds<ind_exp.size]
    ind_exp = ind_exp[keep_inds]
    
    # this is tricky to read but quite simple in concept
    # this remove ind_exp assigned to the same ind_insp
    bad, = np.nonzero(np.diff(ind_exp) == 0)
    keep = np.ones(ind_insp.size, dtype='bool')
    keep[bad + 1] = False
    ind_insp = ind_insp[keep]
    keep = np.ones(ind_exp.size, dtype='bool')
    keep[bad + 1] = False
    ind_exp = ind_exp[keep]

    #~ import matplotlib.pyplot as plt
    # fig, ax = plt.subplots()
    # ax.plot(resp)
    # ax.scatter(ind_insp_no_clean, resp[ind_insp_no_clean], color='m', marker='*', s=100)
    # ax.scatter(ind_dw, resp[ind_dw], color='orange', marker='o', s=30)
    # ax.scatter(ind_insp, resp[ind_insp], color='g', marker='o')
    # ax.scatter(ind_exp, resp[ind_exp], color='r', marker='o')
    # ax.axhline(baseline, color='r')
    # ax.axhline(baseline_insp, color='g')
    # ax.axhline(baseline_dw, color='orange')
    # ax.axhline(q10, color='k')
    # plt.show()


    if ind_insp.size == 0:
        print('no cycle dettected')
        return


    mask = (ind_exp > ind_insp[0]) & (ind_exp < ind_insp[-1])
    ind_exp = ind_exp[mask]

    if inspiration_adjust_on_derivative:
        # lets find local minima on second derivative
        # this can be slow
        delta_ms = 10.
        delta = int(delta_ms * srate / 1000.)
        derivate1 = np.gradient(resp)
        derivate2 = np.gradient(derivate1)
        for i in range(ind_exp.size):
            i0, i1 = ind_insp[i], ind_exp[i]
            i0 = max(0, i0 - delta)
            i1 = i0 + np.argmin(resp[i0:i1])
            d1 = derivate1[i0:i1]
            i1 = i0 + np.argmin(d1)
            if (i1 - i0) >2:
                # find the last crossing zeros in this this short segment
                d2 = derivate2[i0:i1]
                i1 = i0 + np.argmin(d2)
                if (i1 - i0) >2:
                    d2 = derivate2[i0:i1]
                    mask = (d2[:-1] >=0) & (d2[1:] < 0)
                    if np.any(mask):
                        ind_insp[i] = i0 + np.nonzero(mask)[0][-1]

    if ind_insp[-1] >= ind_exp[-1] and ind_insp[-2] >= ind_exp[-1]:
        ind_insp = ind_insp[:-1]
    
    cycles = np.zeros((ind_insp.size - 1, 3), dtype='int64')
    cycles[:, 0] = ind_insp[:-1]
    cycles[:, 1] = ind_exp
    cycles[:, 2] = ind_insp[1:]

    return cycles




########################################
######## COMPUTE RESPI FEATURES ########
########################################



#cycles_init = cycles
def exclude_bad_cycles(respi, cycles_init, srate, exclusion_coeff=1):

    if debug:
        plt.plot(respi, label='respi')
        plt.scatter(cycles_init[:,0], respi[cycles_init[:,0]], color='g', label='inspi')
        plt.scatter(cycles_init[:,1], respi[cycles_init[:,1]], color='r', label='expi')
        plt.scatter(cycles_init[:,2], respi[cycles_init[:,2]], color='b', label='next_inspi')
        plt.legend()
        plt.show()

    #### compute average durations and dispertion for inspi and expi together
    all_inspi_expi = np.diff(cycles_init, axis=1).reshape(-1)
    duration_med = int(np.median(all_inspi_expi))
    duration_mad = int(np.median(np.abs(all_inspi_expi - duration_med)) / 0.6744897501960817)

    #### compute average inspi/expi cycle for cross correlation
    time_vec_mean_cycle = np.arange(-duration_med, duration_med)
    cycle_average_inspi = np.zeros((cycles_init.shape[0], time_vec_mean_cycle.size))
    cycle_average_expi = np.zeros((cycles_init.shape[0], time_vec_mean_cycle.size)) 

    for _inspi_start_i, _inspi_start in enumerate(cycles_init[:,0]):
        _start, _stop = _inspi_start - duration_med, _inspi_start + duration_med  
        if _start < 0 or  _stop < 0 or _start > respi.size or _stop > respi.size:
            continue
        cycle_average_inspi[_inspi_start_i, :] = respi[_start:_stop]

    for _expi_start_i, _expi_start in enumerate(cycles_init[:,1]):
        _start, _stop = _expi_start - duration_med, _expi_start + duration_med  
        if _start < 0 or _stop < 0 or _start > respi.size or _stop > respi.size:
            continue
        cycle_average_expi[_expi_start_i, :] = respi[_start:_stop]

    sel_i_inspi = np.where(np.sum(cycle_average_inspi, axis=1) != 0)[0]
    sel_i_expi = np.where(np.sum(cycle_average_expi, axis=1) != 0)[0]

    cycle_average_inspi = np.median(cycle_average_inspi[sel_i_inspi, :], axis=0)
    cycle_average_expi = np.median(cycle_average_expi[sel_i_expi, :], axis=0)

    if debug:
        plt.plot(time_vec_mean_cycle/srate, cycle_average_inspi)
        plt.plot(time_vec_mean_cycle/srate, cycle_average_expi)
        plt.show()

    #### cross correlate respi with average inspi/expi cycle and identify max
    correlate_inspi = scipy.signal.correlate(respi, cycle_average_inspi, mode='same')
    correlate_expi = scipy.signal.correlate(respi, cycle_average_expi, mode='same')
    
    peak_thresh_inspi = np.std(correlate_inspi)*0.5
    cross_corr_inspi = scipy.signal.find_peaks(correlate_inspi, distance=duration_med, height=peak_thresh_inspi)[0]

    peak_thresh_expi = np.std(correlate_expi)*0.5
    cross_corr_expi = scipy.signal.find_peaks(correlate_expi, distance=duration_med, height=peak_thresh_expi)[0]

    if debug:
        plt.plot(correlate_inspi)
        plt.scatter(cross_corr_inspi, correlate_inspi[cross_corr_inspi], color='r')
        plt.hlines(peak_thresh_inspi, xmin=0, xmax=correlate_inspi.size, color='r')
        plt.show()

        respi_zscore = zscore(respi)
        correlate_inspi_zscore = zscore(correlate_inspi)
        correlate_expi_zscore = zscore(correlate_expi)

        plt.plot(respi_zscore, label='respi')
        plt.plot(correlate_inspi_zscore, label='cross_corr_inspi')
        plt.plot(correlate_expi_zscore, label='cross_corr_expi')
        plt.scatter(cross_corr_inspi, respi_zscore[cross_corr_inspi], color='r', label='cross_corr')
        plt.scatter(cross_corr_expi, respi_zscore[cross_corr_expi], color='r')
        plt.scatter(cycles_init[:,1], respi_zscore[cycles_init[:,1]], color='b', label='cycles_init')
        plt.scatter(cycles_init[:,0], respi_zscore[cycles_init[:,0]], color='b')
        plt.legend()
        plt.show()

    #### correct inspi/expi position based on cross correlation with dispertion coeff
    exclusion_thresh = int(exclusion_coeff * duration_mad)

    inspi_corrected = []

    for inspi_cross_corr_i, inspi_cross_corr in enumerate(cross_corr_inspi):

        _up, _down = inspi_cross_corr + exclusion_thresh, inspi_cross_corr - exclusion_thresh
        mask_thresh = (cycles_init[:,0] <= _up) & (cycles_init[:,0] >= _down)
        if mask_thresh.sum() != 0:
            _inspi_init = int(np.median(cycles_init[:,0][mask_thresh]))
            inspi_corrected.append(_inspi_init)
        else:
            inspi_corrected.append(inspi_cross_corr)

    inspi_corrected = np.array(inspi_corrected)

    if debug:
        plt.plot(respi)
        for inspi_cross_corr in cross_corr_inspi:
            _up, _down = inspi_cross_corr + exclusion_thresh, inspi_cross_corr - exclusion_thresh
            plt.vlines([_up, _down], ymin=respi.min(), ymax=respi.max(), color='r')
        plt.scatter(cycles_init[:,0], respi[cycles_init[:,0]], color='b', label='cycles_init_inspi')
        plt.scatter(cross_corr_inspi, respi[cross_corr_inspi], color='r', label='cross_corr_inspi')
        plt.scatter(inspi_corrected, respi[inspi_corrected], color='g', label='corrected_inspi', marker='x', s=100)
        plt.legend()
        plt.show()

    expi_corrected = []

    for expi_cross_corr_i, expi_cross_corr in enumerate(cross_corr_expi):

        _up, _down = expi_cross_corr + exclusion_thresh, expi_cross_corr - exclusion_thresh
        mask_thresh = (cycles_init[:,1] <= _up) & (cycles_init[:,1] >= _down)
        if mask_thresh.sum() != 0:
            _expi_init = int(np.median(cycles_init[:,1][mask_thresh]))
            expi_corrected.append(_expi_init)
        else:
            expi_corrected.append(expi_cross_corr)

    expi_corrected = np.array(expi_corrected)

    if debug:
        plt.plot(respi)
        for expi_cross_corr in cross_corr_expi:
            _up, _down = expi_cross_corr + exclusion_thresh, expi_cross_corr - exclusion_thresh
            plt.vlines([_up, _down], ymin=respi.min(), ymax=respi.max(), color='r')
        plt.scatter(cycles_init[:,1], respi[cycles_init[:,1]], color='b', label='cycles_init_expi')
        plt.scatter(cross_corr_expi, respi[cross_corr_expi], color='r', label='cross_corr_expi')
        plt.scatter(expi_corrected, respi[expi_corrected], color='g', label='corrected_expi', marker='x', s=100)
        plt.legend()
        plt.show()

    if debug:
        plt.plot(respi)
        plt.scatter(inspi_corrected, respi[inspi_corrected], color='g', label='inspi')
        plt.scatter(expi_corrected, respi[expi_corrected], color='r', label='expi')
        plt.legend()
        plt.show()

    #### verify that detection start and stop by inspi 
    if expi_corrected[0] < inspi_corrected[0]:
        inspi_corrected = np.insert(inspi_corrected, 0, cycles_init[0,0])

    if expi_corrected[-1] > inspi_corrected[-1]:
        expi_corrected = expi_corrected[:-1]

    #### clean by altercating inspi and expi
    inspi_cleaned = np.array([], dtype='int')
    
    for inspi_i, inspi in enumerate(inspi_corrected):

        if inspi_i == 0 or inspi == inspi_corrected[-1]:
            inspi_cleaned = np.append(inspi_cleaned, inspi)
            continue

        mask_expi = expi_corrected > inspi
        next_expi = expi_corrected[mask_expi][0]

        mask_inspi = inspi_corrected < next_expi
        correct_inspi = inspi_corrected[mask_inspi][-1]

        if inspi != correct_inspi:
            continue

        else:
            inspi_cleaned = np.append(inspi_cleaned, inspi)

    expi_cleaned = np.array([], dtype='int')

    for expi_i, expi in enumerate(expi_corrected):

        mask_inspi = expi < inspi_corrected
        next_inspi = inspi_corrected[mask_inspi][0]

        mask_expi = expi_corrected < next_inspi
        correct_expi = expi_corrected[mask_expi][-1]

        if expi != correct_expi:
            continue

        else:
            expi_cleaned = np.append(expi_cleaned, expi)

    if inspi_cleaned.size != expi_cleaned.size+1:
        raise ValueError('!!! bad detection !!!')
    
    if debug:

        plt.plot(respi)
        plt.scatter(inspi_cleaned, respi[inspi_cleaned], color='g', label='inspi')
        plt.scatter(expi_cleaned, respi[expi_cleaned], color='r', label='expi')
        plt.legend()
        plt.show()
    
    cycles_cleaned = np.concatenate((inspi_cleaned[:-1].reshape(-1,1), expi_cleaned.reshape(-1,1), inspi_cleaned[1:].reshape(-1,1)), axis=1)

    if debug:

        plt.plot(respi)
        plt.scatter(cycles_cleaned[:,0], respi[cycles_cleaned[:,0]], color='g', label='inspi')
        plt.scatter(cycles_cleaned[:,1], respi[cycles_cleaned[:,1]], color='r', label='expi')
        plt.scatter(cycles_cleaned[:,2], respi[cycles_cleaned[:,2]], color='b', label='next_inspi')
        plt.legend()
        plt.show()

    #### export fig
    time_vec = np.arange(respi.shape[0])/srate

    inspi_removed = np.array([inspi for inspi in cycles_init[:,0] if inspi not in cycles_cleaned[:,0]])
    expi_removed = np.array([expi for expi in cycles_init[:,1] if expi not in cycles_cleaned[:,1]])
    
    fig_respi_exclusion, ax = plt.subplots(figsize=(18, 10))
    ax.plot(time_vec, respi)

    ax.scatter(cycles_init[:,0]/srate, respi[cycles_init[:,0]], color='g', label='inspi')
    if inspi_removed.size != 0:
        ax.scatter(inspi_removed/srate, respi[inspi_removed], color='k', marker='x', s=100)
    ax.scatter(cycles_cleaned[:,0]/srate, respi[cycles_cleaned[:,0]], color='g')
    ax.scatter(cycles_cleaned[:,-1]/srate, respi[cycles_cleaned[:,-1]], color='g')

    ax.scatter(cycles_init[:,1]/srate, respi[cycles_init[:,1]], color='r', label='expi')
    if expi_removed.size != 0:
        ax.scatter(expi_removed/srate, respi[expi_removed], color='k', marker='x', s=100)
    ax.scatter(cycles_cleaned[:,1]/srate, respi[cycles_cleaned[:,1]], color='r')

    plt.legend()
    # plt.show()

    fig_final, ax = plt.subplots(figsize=(18, 10))
    ax.plot(time_vec, respi)

    ax.scatter(cycles_cleaned[:,0]/srate, respi[cycles_cleaned[:,0]], color='g', label='inspi')
    ax.scatter(cycles_cleaned[:,-1]/srate, respi[cycles_cleaned[:,-1]], color='g')
    ax.scatter(cycles_cleaned[:,1]/srate, respi[cycles_cleaned[:,1]], color='r', label='expi')

    plt.legend()
    # plt.show()

    return cycles_cleaned, fig_respi_exclusion, fig_final


    




############################
######## LOAD DATA ########
############################


def load_respfeatures(sujet):

    #### load data
    os.chdir(path_prep)

    respfeatures_allcond = {}
    respi_allcond = {}

    #cond = 'CHARGE'
    for cond in cond_list:

        respi = load_data_sujet(sujet,cond)[np.where(chan_list == 'pression')[0][0],:]
        respi = scipy.signal.detrend(respi, type='linear')
        respi_allcond[cond] = respi

        params = physio.get_respiration_parameters('human_airflow')
        params['cycle_clean']['low_limit_log_ratio'] = 6
        # params['cycle_detection']['inspiration_adjust_on_derivative'] = True

        diff1 = np.diff(respi)*-1
        peaks, _ = scipy.signal.find_peaks(diff1, prominence=(diff1).std()*3)

        baseline_val = []

        for peak_i, peak in enumerate(peaks):

            if peak_i == peaks.size-1:
                continue

            backward_signal = diff1[peaks[peak_i]:peaks[peak_i+1]][::-1]
            backward_val = np.where(np.diff(backward_signal) > 0)[0][0]
            baseline_val.append(respi[peak-backward_val])

        baseline = np.median(baseline_val)

        if debug:

            plt.plot(respi)
            plt.hlines(np.median(baseline_val), xmin=0, xmax=respi.size, color='r')
            plt.show()

        params['baseline']['baseline_mode'] = 'manual'
        params['baseline']['baseline'] = baseline

        respi_clean, resp_features_i = physio.compute_respiration(raw_resp=respi, srate=srate, parameters=params)

        respi -= baseline

        time_vec = np.arange(respi.size)/srate

        fig_final, ax = plt.subplots(figsize=(18, 10))
        ax.plot(time_vec, respi)

        ax.scatter(resp_features_i['inspi_index'].values/srate, respi[resp_features_i['inspi_index'].values], color='g', label='inspi')
        ax.scatter(resp_features_i['expi_index'].values/srate, respi[resp_features_i['expi_index'].values], color='r', label='expi')
        ax.scatter(resp_features_i['next_inspi_index'].values/srate, respi[resp_features_i['next_inspi_index'].values], color='g', label='inspi')

        plt.legend()
        # plt.show()

        # cycles = physio.detect_respiration_cycles(respi, srate, baseline_mode='median',
        #                                           baseline=None, epsilon_factor1=10, epsilon_factor2=5, inspiration_adjust_on_derivative=False)
        
        # cycles = debugged_detect_respiration_cycles(respi, srate, baseline_mode='median',
        #                                             baseline=None, epsilon_factor1=10, epsilon_factor2=5, inspiration_adjust_on_derivative=False)
        
        # if debug:

        #     fig, ax = plt.subplots()
        #     ax.plot(respi)
        #     ax.scatter(cycles[:,0], respi[cycles[:,0]], color='g')
        #     plt.show()

        # cycles, fig_respi_exclusion, fig_final = exclude_bad_cycles(respi, cycles, srate, exclusion_coeff=1)
            
        # if debug:

        #     fig, ax = plt.subplots()
        #     ax.plot(respi)
        #     ax.scatter(cycles[:,0], respi[cycles[:,0]], color='r')
        #     plt.show()

        # #### get resp_features
        # resp_features_i = physio.compute_respiration_cycle_features(respi, srate, cycles, baseline=None)

        # select_vec = np.ones((resp_features_i.index.shape[0]), dtype='int')
        # resp_features_i.insert(resp_features_i.columns.shape[0], 'select', select_vec)
        
        respfeatures_allcond[cond] = [resp_features_i, fig_final]

    return respi_allcond, respfeatures_allcond




####################################
######## PLOT MEAN RESPI ########
####################################

def plot_mean_respi(sujet):

    time_vec = np.arange(stretch_point_ERP)
    colors_respi = {'VS' : 'b', 'CHARGE' : 'r'}
    colors_respi_sem = {'VS' : 'c', 'CHARGE' : 'm'}

    respi_allcond, respfeatures = load_respfeatures(sujet)
    sem_allcond = {}
    lim = {'min' : np.array([]), 'max' : np.array([])} 

    #### load
    #cond = 'VS'
    for cond in cond_list:

        respi_stretch = stretch_data(respfeatures[cond][0], stretch_point_ERP, respi_allcond[cond], srate)[0]
        respi_allcond[cond] = respi_stretch.mean(axis=0)
        sem_allcond[cond] = respi_stretch.std(axis=0)/np.sqrt(respi_stretch.shape[0])
        lim['min'], lim['max'] = np.append(lim['min'], respi_allcond[cond].min()-sem_allcond[cond]), np.append(lim['max'], respi_allcond[cond].max()+sem_allcond[cond])

    #### plot
    fig, ax = plt.subplots()

    #cond = 'VS'
    for cond in cond_list:

        ax.plot(time_vec, respi_allcond[cond], color=colors_respi[cond], label=cond)
        ax.fill_between(time_vec, respi_allcond[cond]+sem_allcond[cond], respi_allcond[cond]-sem_allcond[cond], alpha=0.25, color=colors_respi_sem[cond])

    ax.vlines(stretch_point_ERP/2, ymin=lim['min'].min(), ymax=lim['max'].max(), color='r')
    plt.ylim(lim['min'].min(), lim['max'].max())
    plt.title(sujet)
    plt.legend()
    # plt.show()

    os.chdir(os.path.join(path_results, 'RESPI', 'plot'))
    plt.savefig(f"{sujet}_respi_mean.png")

    plt.close('all')















############################
######## EXECUTE ########
############################



if __name__ == '__main__':

    
    # sujet_list = ['01NM_MW', '02NM_OL', '03NM_MC', '04NM_LS', '05NM_JS', '06NM_HC', '07NM_YB', '08NM_CM', '09NM_CV', '10NM_VA', '11NM_LC', '12NM_PS', '13NM_JP', '14NM_LD',
    #           '15PH_JS',  '16PH_LP',  '17PH_SB',  '18PH_TH',  '19PH_VA',  '20PH_VS',
    #           '21IL_NM', '22IL_DG', '23IL_DM', '24IL_DJ', '25IL_DC', '26IL_AP', '27IL_SL', '28IL_LL', '29IL_VR', '30IL_LC', '31IL_MA', '32IL_LY', '33IL_BA', '34IL_CM', '35IL_EA', '36IL_LT',
    #           '37DL_05', '38DL_06', '39DL_07', '40DL_08', '41DL_11', '42DL_12', '43DL_13', '44DL_14', '45DL_15', '46DL_16', '47DL_17', '48DL_18', '49DL_19', '50DL_20', '51DL_21', '52DL_22',
    #           '53DL_23', '54DL_24', '55DL_25', '56DL_26', '57DL_27', '58DL_28', '59DL_29', '60DL_30', '61DL_31', '62DL_32', '63DL_34',
    #           ]

    sujet = '49DL_19'

    for sujet in sujet_list:

        if os.path.exists(os.path.join(path_results, 'RESPI', 'count', f'{sujet}_count_cycles.xlsx')):
            print(f"{sujet} ALREADY COMPUTED")
            continue
        else:
            print(sujet)
        
        respi_allcond, respfeatures_allcond = load_respfeatures(sujet)

        ########################################
        ######## VERIF RESPIFEATURES ########
        ########################################
        
        if debug == True :

            cond = 'VS'
            cond = 'CHARGE' 

            respfeatures_allcond[cond][1].show()
            # respfeatures_allcond[cond][2].show()

        ########################################
        ######## EDIT CYCLES SELECTED ########
        ########################################

        # respfeatures_allcond = edit_df_for_sretch_cycles_deleted(respi_allcond, respfeatures_allcond)

        #### generate df
        df_count_cycle = pd.DataFrame(columns={'sujet' : [], 'cond' : [], 'count' : []})

        for cond in cond_list:
            
            data_i = {'sujet' : [sujet], 'cond' : [cond], 'count' : [respfeatures_allcond[cond][0].shape[0]]}
            df_i = pd.DataFrame(data_i, columns=data_i.keys())
            df_count_cycle = pd.concat([df_count_cycle, df_i])

        #### export
        os.chdir(os.path.join(path_results, 'RESPI', 'count'))
        df_count_cycle.to_excel(f'{sujet}_count_cycles.xlsx')

        if debug :

            for cond in cond_list:

                _respi = respi_allcond[cond]
                plt.plot(_respi)
                plt.vlines(respfeatures_allcond[cond][0]['inspi_index'].values, ymin=_respi.min(), ymax=_respi.max(), color='r')
                plt.title(f"{cond}")
                plt.show()

        ################################
        ######## SAVE FIG ########
        ################################

        for cond in cond_list:

            os.chdir(os.path.join(path_results, 'RESPI', 'respfeatures'))
            respfeatures_allcond[cond][0].to_excel(f"{sujet}_{cond}_respfeatures.xlsx")
            
            os.chdir(os.path.join(path_results, 'RESPI', 'detection'))
            respfeatures_allcond[cond][1].savefig(f"{sujet}_{cond}_fig0.jpeg")
            # respfeatures_allcond[cond][2].savefig(f"{sujet}_{cond}_fig1.jpeg")

        plt.close('all')

        ################################
        ######## PLOT MEAN RESPI ########
        ################################

        plot_mean_respi(sujet)





    ####################################
    ######## AGGREGATES COUNT ########
    ####################################
        
    df_count_cycle = pd.DataFrame(columns={'sujet' : [], 'cond' : [], 'count' : []})

    os.chdir(os.path.join(path_results, 'RESPI', 'count'))

    for sujet in sujet_list:

        _df_count_cycle = pd.read_excel(f'{sujet}_count_cycles.xlsx')
        df_count_cycle = pd.concat([df_count_cycle, _df_count_cycle], axis=0)

    df_count_cycle = df_count_cycle.drop(columns='Unnamed: 0')

    ncount_selsujet = []

    for ncount in np.arange(40, 300, 1):
        df_sel_ncycle = df_count_cycle.query(f"count >= {ncount}")
        nsujet = [sujet for sujet in sujet_list if (df_sel_ncycle['sujet'].values == sujet).sum() == 2] 
        ncount_selsujet.append(len(nsujet))

    plt.plot(np.arange(40, 300, 1), ncount_selsujet)
    plt.show()

    df_sel_ncycle = df_count_cycle.query(f"count >= 100")
    sujet_list_fc = [sujet for sujet in sujet_list if (df_sel_ncycle['sujet'].values == sujet).sum() == 2]
    nsujet = len(sujet_list_fc)

    plt.hist(df_count_cycle.query(f"cond == 'VS'")['count'].values, bins=20, label='VS', alpha=0.7)
    plt.hist(df_count_cycle.query(f"cond == 'CHARGE'")['count'].values, bins=20, label='CHARGE', alpha=0.7)
    plt.xlabel('count')
    plt.ylabel('n_sujet')
    plt.title('cycle_count')
    plt.legend()
    plt.show()

    plt.savefig('allsujet_cycle_count.png')