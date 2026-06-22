



from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *


import plotly.express as px
from plotly.subplots import make_subplots









################################
######## GENERATE ALLDF ########
################################

def generate_all_df_Pxx_wholecycle_patient_wise():

    phase_cycle_list = ['inspi', 'expi']
    factor_select = ['A+', 'R-', 'C-']

    path_import_res_Lmm = os.path.join(path_precompute, 'EXPORT_DF', 'patient_wise', 'LMM', 'df_res')

    df_R_Pxx = []

    for phase_cycle in phase_cycle_list:

        for band in freq_band_dict:

            for ROI in ROI_short_list:

                _df = pd.read_excel(os.path.join(path_import_res_Lmm, f"RES_VCCP_{band}_{phase_cycle}_{ROI}_Pxx.xlsx"))
                _df = _df.query(f"term in {factor_select}").reset_index(drop=True)
                _df_add = pd.concat([pd.DataFrame({'phase_cycle' : [phase_cycle]*_df.shape[0], 'band' : [band]*_df.shape[0], 'ROI' : [ROI]*_df.shape[0]}), _df[['term', 'estimate', 'p.value']]], axis=1)
                
                df_R_Pxx.append(_df_add)

    df_R_Pxx = pd.concat(df_R_Pxx)
    df_R_Pxx = df_R_Pxx.rename(columns={'p.value' : 'pvalue'}).reset_index(drop=True)

    return df_R_Pxx









def p_to_stars(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return ''












########################
######## PXX ########
########################




def export_res_Pxx_patient_wise(sujet):

    #### config
    phase_cycle_list = ['inspi', 'expi']
    path_export_res = os.path.join(path_results, 'Pxx', 'LMM', 'VCCP', 'patientwise')

    #### load df
    df_R_Pxx = generate_all_df_Pxx_wholecycle_patient_wise()

    chanlist, chanlist_aux, chanlist_all, localist = get_chanlist_localist(sujet)

    label_name_nplot = {}
    
    for ROI in ROI_short_list:

        _nplot = localist.query(f"loca == '{ROI}'").size
        label_name_nplot[ROI] = f"{ROI} c({_nplot})"

    for ROI in ROI_short_list:
        df_R_Pxx['ROI'] = df_R_Pxx['ROI'].replace({ROI: label_name_nplot[ROI]})

    ROI_plot_short_list = df_R_Pxx['ROI'].unique().tolist()

    #### inspi/expi
    color_map = {
        "inspi": "#1f77b4",
        "expi":  "#d62728" 
    }

    LMM_param_list = ['A+', 'R-', 'C-']
    
    #band = 'theta'
    for band in freq_band_dict:


        fig = make_subplots(
            rows=len(LMM_param_list),
            cols=1,
            shared_xaxes=True,
            subplot_titles=LMM_param_list
        )

        y_max = df_R_Pxx.query(f"band == '{band}' and ROI in {ROI_plot_short_list}").copy()["estimate"].abs().max()

        for r, LMM_param in enumerate(LMM_param_list, start=1):

            _df_plot = df_R_Pxx.query(f"band == '{band}' and term == '{LMM_param}' and ROI in {ROI_plot_short_list}").copy()

            _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

            for phase_cycle in _df_plot["phase_cycle"].unique():

                df_term = _df_plot[_df_plot["phase_cycle"] == phase_cycle]

                fig.add_bar(
                    x=df_term["ROI"],
                    y=df_term["estimate"],
                    name=phase_cycle,
                    text=df_term["sig"],
                    textposition="outside",
                    marker_color=color_map[phase_cycle],
                    row=r,
                    col=1
                )

                margin = y_max * 0.50   # 15% extra space

                fig.update_yaxes(
                    range=[-y_max - margin, y_max + margin],
                    row=r,
                    col=1
                )

                # fig.update_xaxes(
                #     categoryorder="array",
                #     categoryarray=ROI_order,
                #     row=r,
                #     col=1
                # )

        fig.update_layout(
            title=f"{band} all ROI",
            template="simple_white",
            barmode="group",
            xaxis_tickangle=-45,
            yaxis_title="estimate",
            legend_title="phase_cycle",
            height=320 * len(phase_cycle_list),  # scale height
            width=500,
        )

        # fig.show()

        fig.write_html(
            os.path.join(path_export_res, f"VCCP_{band}_threshROI_LMM.html"),
            include_plotlyjs="cdn"
        )








def export_res_Pxx_allpatient():
 
    phase_cycle_list = ['inspi', 'expi']

    #### load df
    df_R_Pxx = generate_all_df_Pxx_wholecycle()

    df_loca_allsujet = get_df_loca_allsujet()

    label_name_nsujet = {}
    
    for ROI in ROI_short_list:

        _nsujet = df_loca_allsujet.query(f"loca == '{ROI}'")['sujet'].unique().size
        label_name_nsujet[ROI] = f"{ROI} s({_nsujet})"

    for ROI in ROI_short_list:
        df_R_Pxx['ROI'] = df_R_Pxx['ROI'].replace({ROI: label_name_nsujet[ROI]})

    df_ROI_scount = df_loca_allsujet.groupby('loca').nunique('sujet')['sujet'].reset_index(name='count')

    for ROI in ROI_short_list:
        df_ROI_scount['loca'] = df_ROI_scount['loca'].replace({ROI: label_name_nsujet[ROI]})

    localist_thresh = df_ROI_scount.query(f"count >= {sujet_thresh}")['loca'].values.tolist()

    ROI_plot_short_list = []

    for _ROI_thresh in localist_thresh:
        for _ROI_short in ROI_short_list:
            if _ROI_thresh.find(_ROI_short) != -1:
                ROI_plot_short_list.append(_ROI_thresh)

    #### plot
        #### allROI
    # phase_cycle = 'whole'

    # for band in freq_band_dict:
        
    #     _df_plot = df_R_Pxx.query(f"band == '{band}' and term != '(Intercept)' and phase_cycle == '{phase_cycle}'").copy()
    #     _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

    #     fig = px.bar(_df_plot, x="ROI", y="estimate", color="term", barmode="group", text="sig", title=f"{band} {phase_cycle} allROI")

    #     fig.update_layout(xaxis_tickangle=-45, yaxis_title="estimate", xaxis_title="ROI", legend_title="term")
    #     fig.update_traces(textposition="outside", cliponaxis=False)

    #     # fig.show()

    #     outdir = os.path.join(path_results, "LMM", "fig", 'Pxx')
    #     fig.write_html(os.path.join(outdir, f"{phase_cycle}_{band}_allROI_barplot_LMM.html"),
    #                 include_plotlyjs="cdn")

    #     #### signiROI 
    # phase_cycle = 'whole'
    
    # for band in freq_band_dict:
        
    #     ROI_signi_sel = df_R_Pxx.query(f"band == '{band}' and term == 'respoc:statechl' and pvalue < 0.05 and phase_cycle == '{phase_cycle}'")["ROI"].values
    #     _df_plot = df_R_Pxx.query(f"band == '{band}' and term != '(Intercept)' and ROI in {ROI_signi_sel.tolist()}").copy()

    #     _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

    #     fig = px.bar(_df_plot, x="ROI", y="estimate", color="term", barmode="group", text="sig", title=f"{band} {phase_cycle} signiROI")

    #     fig.update_layout(xaxis_tickangle=-45, yaxis_title="estimate", xaxis_title="ROI", legend_title="term")
    #     fig.update_traces(textposition="outside", cliponaxis=False)

    #     # fig.show()

    #     outdir = os.path.join(path_results, "LMM", "fig", 'Pxx')
    #     fig.write_html(os.path.join(outdir, f"{phase_cycle}_{band}_signiROI_barplot_LMM.html"),
    #             include_plotlyjs="cdn")

    #     #### threshROI
    # phase_cycle = 'whole'

    # for band in freq_band_dict:
        
    #     _df_plot = df_R_Pxx.query(f"ROI in {localist_thresh} and band == '{band}' and term != '(Intercept)' and phase_cycle == '{phase_cycle}'").copy()
    #     _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

    #     fig = px.bar(_df_plot, x="ROI", y="estimate", color="term", barmode="group", text="sig", title=f"{band} {phase_cycle} threshROI")

    #     fig.update_layout(xaxis_tickangle=-45, yaxis_title="estimate", xaxis_title="ROI", legend_title="term")
    #     fig.update_traces(textposition="outside", cliponaxis=False)

    #     # fig.show()

    #     outdir = os.path.join(path_results, "LMM", "fig", 'Pxx')
    #     fig.write_html(os.path.join(outdir, f"{phase_cycle}_{band}_threshROI_barplot_LMM.html"),
    #                 include_plotlyjs="cdn")


        #### inspi/expi
    color_map = {
        "inspi": "#1f77b4",
        "expi":  "#d62728" 
    }

    LMM_param_list = ['respoc', 'statechl', 'respoc:statechl']
    
    ROI_order = ['Amygdala s(8)', 'Hippocampus s(7)', 'insula-ant s(5)', 'insula-pos s(4)', 'lateralorbitofrontal s(5)', 'medialorbitofrontal s(5)', 'postcentral s(3)', 'precentral s(4)']

    #band = 'theta'
    for band in freq_band_dict:


        fig = make_subplots(
            rows=len(LMM_param_list),
            cols=1,
            shared_xaxes=True,
            subplot_titles=LMM_param_list
        )

        y_max = df_R_Pxx.query(f"band == '{band}' and term != '(Intercept)' and phase_cycle != 'whole' and ROI in {ROI_plot_short_list}").copy()["estimate"].abs().max()

        for r, LMM_param in enumerate(LMM_param_list, start=1):

            _df_plot = df_R_Pxx.query(f"band == '{band}' and term == '{LMM_param}' and phase_cycle != 'whole' and ROI in {ROI_plot_short_list}").copy()

            _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

            for phase_cycle in _df_plot["phase_cycle"].unique():

                df_term = _df_plot[_df_plot["phase_cycle"] == phase_cycle]

                fig.add_bar(
                    x=df_term["ROI"],
                    y=df_term["estimate"],
                    name=phase_cycle,
                    text=df_term["sig"],
                    textposition="outside",
                    marker_color=color_map[phase_cycle],
                    row=r,
                    col=1
                )

                margin = y_max * 0.50   # 15% extra space

                fig.update_yaxes(
                    range=[-y_max - margin, y_max + margin],
                    row=r,
                    col=1
                )

                fig.update_xaxes(
                    categoryorder="array",
                    categoryarray=ROI_order,
                    row=r,
                    col=1
                )

        fig.update_layout(
            title=f"{band} all ROI",
            template="simple_white",
            barmode="group",
            xaxis_tickangle=-45,
            yaxis_title="estimate",
            legend_title="phase_cycle",
            height=320 * len(phase_cycle_list),  # scale height
            width=500,
        )

        # fig.show()

        outdir = os.path.join(path_results, "LMM", "fig", "Pxx", "summary")
        fig.write_html(
            os.path.join(outdir, f"{band}_threshROI_LMM.html"),
            include_plotlyjs="cdn"
        )

    #### export df data
    df_export = df_R_Pxx.query(f"phase_cycle != 'whole' and ROI in {ROI_plot_short_list} and term != '(Intercept)'")
    
    output_df = os.path.join(path_results, "LMM", "fig", "Pxx", "summary")
    filename = os.path.join(output_df, "df_Pxx_LMM.xlsx")
    df_export.to_excel(filename)
















################################
######## EXECUTE ########
################################


if __name__ == '__main__':


    for sujet in sujet_list:
        export_res_Pxx_patient_wise(sujet)

    export_res_Pxx_patient_wise()
    
    export_res_Pxx()

        


    























