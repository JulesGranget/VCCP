
clear all; close all

%restoredefaultpath
addpath("/home/david/MATLAB Add-Ons/Collections/FieldTrip")
ft_defaults

addpath('/home/david/Documents/VCCP/HBML_JH_iEEG_preproc')

addpath /home/david/Downloads/scripts
addpath /home/david/Downloads/scripts/fieldtrip-20230118/external/dmlt/external/murphy/KPMtools
addpath /home/david/Documents/RRET/BU/fct_jose/
addpath /home/david/Documents/EPIPE/Matlab/EXTERNAL/TDT/MatlabSDK/TDTSDK/TDTbin2mat
addpath(genpath('/home/david/Documents/RRET/BU/EPIPROC-master'))
addpath('/home/david/Documents/RRET/BU/FT_private/')
addpath('/home/david/Documents/RRET/BU')
addpath '/home/david/Documents/RRET/latestscripts_JA'
addpath('/home/david/Documents/RRET/samruddhi_scripts')
addpath('/home/david/Documents/Projects/VCCP/HBML_JH_iEEG_preproc')
addpath('/usr/local/MATLAB/R2025a/toolbox')



%% -----------------  CONFIG ---------------
pat = 'LH018'; %LH018; %NS215; %NS213; %NS131_02; %'NS211';%'LH013';%'NS200';%'NS057_3';%'LH013';%'NS189';%'LH014_3';%'NS191';%'LH014_3'; 'NS057_3'

plot_display = false;

path_base = '/home/david/Documents/Projects/VCCP'; 
path_raw_data = fullfile(path_base, "data_raw");
path_out_extraction = fullfile(path_base, "data_extraction");

fs = 512;
fs_resampled = 500;
fs_sl = 100;

params_preproc = [];
params_preproc.filter_notch    = [60 120 180];
params_preproc.filter_bandpass = [0.01 200];

% patient information
allpatient_et_data_presence = containers.Map( ...
    {'LH018', 'test'}, ...
    [1,       0] );

allpatient_rec_system = containers.Map( ...
    {'LH018', 'test'}, ...
    {'NATUS',       'TDT'} );


allpatient_SLpath = containers.Map( ...
    {'LH018', 'test'}, ...
    {'SmartLab_2-12-2026_LH018_VCCL.csv',       'TDT'} );




if strcmp(pat,'LH018')
    
    filename_data_raw = fullfile(path_raw_data, pat, "tdt", "LH018_VCCL.edf");

    ttl_chan = 259;

    pre_chunk_sync_SL = round(1.012 * 1e4);
    post_chunk_sync_SL = round(3.9968 * 1e5);

    pre_chunk_sync_TDT = 127500;
    post_chunk_sync_TDT = round(2.074 * 1e6);

    pat_correction_SL = [218850 219050;
                  ];

    pat_correction_TDT = [218850 219050;
                  ];

    chanlist_aux = {'respi' 'O2' 'CO2'};
    
    if allpatient_et_data_presence(pat)
    
        currentBlockRecs = 1:16; % rec# (examples in /home/david/Documents/RRET/patients/NS206/logs/et/preproc)
        currentBlockTrs = 1:16; % tr# (examples in /home/david/Documents/RRET/patients/NS206/logs/et/preproc)
        if length(currentBlockRecs)~=length(currentBlockTrs)
            error('ET recs and trs input incorrectly')
        end
        
        params.et = [];
        params.et.device = 'neon';
        params.et.dir = '/home/david/Documents/RRET/patients/LH018/vccp/logs/et/preproc/';
        params.et.fs = 200;
        params.et.currentBlockRecs = currentBlockRecs;
        params.et.currentBlockTrs = currentBlockTrs;
        params.et.sinktimes = [[fullfile(path_raw_data, pat, "tdt") '/LH018_Experiment_Log.xlsx','/LH018_VCCP_sink'] '_et'];
    
    else
        params.et = 'none';
    end
end


%% ----------------- EXTRACT ---------------------------------------------------------------------


if strcmp(allpatient_rec_system(pat),'TDT')
    [ecog] = TDT2ecog(params.paramsfile);
    chanlist = extract_data.label;
elseif strcmp(allpatient_rec_system(pat),'NATUS')
    [extract_data, eeg_raw] = LIJ_edfread(filename_data_raw);
    extract_data.eeg = eeg_raw;
    save(fullfile(path_out_extraction, pat, sprintf('%s_data_raw.mat', pat)), 'extract_data', '-v7.3')  
    chanlist = extract_data.label;
end


%% --------------------VISUALIZATION ------------------------------------------------------------------------
Visualization = 0;
if Visualization

    cfg =[];
    cfg.channel = 'all';%'LIa*';%{'LFo*'};%'all';
    eeg = eeg_raw;

    % notch
    selchan = 20;

    eeg = fillmissing(eeg, 'linear', 2); %added to handle NaNs in eeg
    eeg = double(eeg); %added due to some precision errors during filtering (inf/nans)

    [eeg1] = BU_filternoise_harmonics(eeg,fs,[],selchan);
    [eeg1] = BU_filternoise_harmonics(eeg,fs,params_preproc.filter_notch,selchan);

    % filter eeg (bandpas)
    n = 2; Wn = [params_preproc.filter_bandpass]/(fs/2);[b,a] = butter(n,Wn);
    eeg2 = filtfilt(b,a,double(eeg1'))';

    figure
    plot(demean(eeg(5,fs:fs*10)),'b');hold on
    plot(demean(eeg1(5,fs:fs*10)),'m');hold on
    plot(demean(eeg2(5,fs:fs*10)),'r');hold on

    % remove mean of each trial
    n = size(eeg, 1);
    for a=1:n
       mn=mean(eeg(a,:));
       eeg(a,:)=eeg(a,:)-repmat(mn,size(eeg(a,:),1),1);
    end
    
    cfg.continuous='yes';
    cfg.viewmode = 'vertical';
    cfg.position = [0.5 0.5 900 600];
    cfg.blocksize = 30;
    cfg.preproc.demean = 'yes';

    data_visu = [];
    data_visu.label = extract_data.label;          
    data_visu.fsample = fs;                
    data_visu.trial = {eeg};                
    data_visu.time = {(0:size(eeg,2)-1)/fs};
    
    artf=ft_databrowser(cfg,data_visu);
end


%% ------------------------ PRE-PROCESSING ------------------------------------------------
% filter (band-pass,notch) and downsamples
% make all signals equal length (iEEG, resp, ekg,...)
% adds smartlab data

% band-pass filtering
figure;plot(eeg_raw(1,fs:fs*10));hold on;
n = 2; Wn = [params_preproc.filter_bandpass]/(fs/2);[b,a] = butter(n,Wn);%jos
eeg_prep_BP = filtfilt(b,a,double(eeg_raw)')';
plot(eeg_prep_BP(1,fs:fs*10),'r');hold on;
filter_band = params_preproc.filter_bandpass;
if filter_band(1) == 0;n = 2; wn = filter_band(2)/(fs/2);[b,a] = butter(n,wn,'low');end%%%low-pass filter <200Hz
if filter_band(1) ~= 0;n = 2; wn = [filter_band]/(fs/2);[b,a] = butter(n,wn);end%%%band-pass

% notch
selchan_plot = 1;
plotting=1;
filternoise_harmonics = [params_preproc.filter_notch];
eeg_prep_BP_Notch = BU_filternoise_harmonics(eeg_prep_BP,fs,filternoise_harmonics,selchan_plot,plotting);

% downsampling eeg
[p, q] = rat(fs_resampled/fs);
eeg_prep_BP_Notch_DS = resample(eeg_prep_BP_Notch.', p, q); % transpose matrix, resample use 1st dimension
eeg_prep_BP_Notch_DS = eeg_prep_BP_Notch_DS.'; % transpose back 
clear newfs oldfs

eeg_prep = eeg_prep_BP_Notch_DS;


%% ---------------- SMARTLAB EXTRACT  ------------------------


SL_raw_sig = readtable(fullfile(path_raw_data, pat, "logs", "sl", allpatient_SLpath(pat)));
respi_SL_raw = table2array(SL_raw_sig(:, {'Module2_RawPressure'})).';
O2 = table2array(SL_raw_sig(:, {'InternalSensor_O2Concentration'})).';
CO2 = table2array(SL_raw_sig(:, {'SerialPort1_CO2Concentration'})).';
SL_TTL = table2array(SL_raw_sig(:, {'DigitalI_O_Input1'})).';

respi_SL_raw = respi_SL_raw(~isnan(respi_SL_raw));
O2 = O2(~isnan(O2));
CO2 = CO2(~isnan(CO2));
SL_TTL = SL_TTL(~isnan(SL_TTL));

if debug:
    figure
    plot(pressure, 'b')
end

%% ---------------- TDT TTL EXTRACT  ------------------------


TDT_TTL = eeg_prep(ttl_chan,:);% 'demean' brings signal to zero


%% ---------------- SYNCHRONIZE SL  ------------------------


debug = false;

mod1 = table2array(SL_raw_sig(:, {'Module1_RawFlow'})).';
mod4 = table2array(SL_raw_sig(:, {'Module4_RawFlow'})).';

respi_SL = -mod4(2:end) + mod1(2:end);

if debug
    figure;
    plot(mod1(2:end)); hold on;
    plot(mod4(2:end));
    legend({'mod1', 'mod4'});

    figure;
    plot(respi_SL);
end

if debug
    figure;
    plot(zscore(TDT_TTL));
    title('TTL');

    figure;
    plot(zscore(SL_TTL));
    title('SL');
end


%% config

sr_TDT = fs_resampled;
sr_SL = fs_sl;


%% chunk and clean SL_TTL

if debug
    figure;
    plot(zscore(respi_SL)); hold on;
    plot(zscore(SL_TTL));
    legend({'respi SL', 'SL TTL'});
end

if debug
    pre_chunk_sync_SL = round(1.012 * 1e4);
    post_chunk_sync_SL = round(3.9968 * 1e5);
end

respi_chunk = respi_SL(pre_chunk_sync_SL:post_chunk_sync_SL);
SL_TTL_chunk = SL_TTL(pre_chunk_sync_SL:post_chunk_sync_SL);

if debug
    figure;
    plot(zscore(respi_chunk)); hold on;
    plot(zscore(SL_TTL_chunk));
    legend({'respi chunk', 'SL TTL chunk'});
end

if debug
    pat_correction_SL = [218850 219050;];
end

for k = 1:size(pat_correction_SL,1)
    SL_TTL_chunk(pat_correction_SL(k,1):pat_correction_SL(k,2)) = 0;
end

if debug:
    figure;
    plot(zscore(respi_chunk)); hold on;
    plot(zscore(SL_TTL_chunk));
    legend({'respi chunk', 'SL TTL corrected'});
end

SL_TTL_chunk_diff = diff(SL_TTL_chunk) * -1;
SL_TTL_chunk_diff(SL_TTL_chunk_diff < 0) = 0;
SL_TTL_clean = [SL_TTL_chunk_diff 0];

[~, peaks_SL_raw] = findpeaks(SL_TTL_clean, 'MinPeakHeight', 0.25 * max(SL_TTL_clean), ...
    'MinPeakDistance', sr_SL * 20);

peaks_SL = peaks_SL_raw - peaks_SL_raw(1);
SL_TTL_clean_corrected = SL_TTL_clean(peaks_SL_raw(1):peaks_SL_raw(end));

SL_correction_factor = [pre_chunk_sync_SL + peaks_SL_raw(1), pre_chunk_sync_SL + peaks_SL_raw(end)];

peaks_SL_ordered = reshape(peaks_SL, 2, []).';

if debug
    figure;
    plot(SL_TTL_clean_corrected); hold on;
    xline(peaks_SL_ordered(:,1), 'g');
    xline(peaks_SL_ordered(:,2), 'r');
end

%% chunk and clean TDT_TTL

if debug
    figure;
    plot(zscore(TDT_TTL));
    title('TDT TTL');
    
    figure;
    plot(zscore(SL_TTL_chunk));
    title('SL TTL chunk');
end

if debug
    pre_chunk_sync_TDT = 127500;
    post_chunk_sync_TDT = round(2.074 * 1e6);
end

TDT_TTL_chunk = TDT_TTL(pre_chunk_sync_TDT:post_chunk_sync_TDT);

if debug
    figure;
    plot(TDT_TTL_chunk);
end

if debug
    pat_correction_TDT = [round(1.091 * 1e6) round(1.093 * 1e6);];
end

for k = 1:size(pat_correction_TDT,1)
    TDT_TTL_chunk(pat_correction_TDT(k,1):pat_correction_TDT(k,2)) = 0;
end

[~, peaks_TDT_raw] = findpeaks(TDT_TTL_chunk, 'MinPeakHeight', 0.25 * max(TDT_TTL), ...
    'MinPeakDistance', sr_TDT * 20);

time_vec_TDT = linspace(0, 1, numel(TDT_TTL_chunk));

if debug
    figure;
    plot(time_vec_TDT, zscore(TDT_TTL_chunk)); hold on;
    scatter(time_vec_TDT(peaks_TDT_raw), zscore(TDT_TTL_chunk(peaks_TDT_raw)), 'r');
    legend({'TDT', 'peaks'});
end

TDT_TTL_clean = zeros(size(TDT_TTL_chunk));
TDT_TTL_clean(peaks_TDT_raw) = 1;
TDT_TTL_clean_corrected = TDT_TTL_clean(peaks_TDT_raw(1):peaks_TDT_raw(end));
peaks_TDT = peaks_TDT_raw - peaks_TDT_raw(1);
peaks_TDT_ordered = reshape(peaks_TDT, 2, []).';

TDT_correction_factor = [pre_chunk_sync_TDT + peaks_TDT_raw(1), pre_chunk_sync_TDT + peaks_TDT_raw(end)];

if debug
    figure;
    plot(TDT_TTL_clean_corrected); hold on;
    xline(peaks_TDT_ordered(:,1), 'g');
    xline(peaks_TDT_ordered(:,2), 'r');
end


%% prepare data to sync

idx_start = SL_correction_factor(1);
idx_stop = SL_correction_factor(end);

data_tosync = [
    respi_SL(idx_start:idx_stop);
    O2(idx_start:idx_stop);
    CO2(idx_start:idx_stop);
];

debug = false;
data_synced = sync_SL_on_TDT(data_tosync, peaks_TDT, peaks_SL, debug);

%% ---------------- EXPORT  ------------------------

eeg_export = [eeg_prep(:, TDT_correction_factor(1):TDT_correction_factor(end)); data_synced];
save(fullfile(path_out_extraction, pat, sprintf('%s_data_prep.mat', pat)), 'eeg_export', '-v7.3')

condlist = readtable(fullfile(path_base, "data_raw", pat, sprintf('%s_condlist.xlsx', pat)));
condlist = array2table(condlist.Type);
condlist.Properties.VariableNames = {'cond'};

trigger_export = array2table(reshape(peaks_TDT,2,[]).');
trigger_export.Properties.VariableNames = {'start','stop'};
trigger_export = [condlist trigger_export];
writetable(trigger_export,fullfile(path_out_extraction, pat, sprintf('%s_trigger.xlsx', pat)));

chanlist_export = cell2table([chanlist chanlist_aux].');
chanlist_export.Properties.VariableNames = {'chan'};
writetable(chanlist_export,fullfile(path_out_extraction, pat, sprintf('%s_chanlist.xlsx', pat)));
