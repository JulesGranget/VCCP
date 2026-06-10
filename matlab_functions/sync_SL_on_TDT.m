function data_synced = sync_SL_on_TDT(data_tosync, peaks_TDT, peaks_SL, debug)

    if debug
        figure;
        plot(data_tosync(1,:));
    end

    time_vec_SL_cell = {};
    time_vec_TDT_cell = {};
    
    for seg_i = 1:(numel(peaks_SL)-1)
    
        n_SL  = peaks_SL(seg_i+1) - peaks_SL(seg_i) + 1;
        n_TDT = peaks_TDT(seg_i+1) - peaks_TDT(seg_i) + 1;
    
        time_vec_SL_cell{end+1} = linspace(seg_i-1, seg_i, n_SL);
        time_vec_TDT_cell{end+1} = linspace(seg_i-1, seg_i, n_TDT);
    
        % remove last sample to avoid overlap between segments
        time_vec_SL_cell{end}(end) = [];
        time_vec_TDT_cell{end}(end) = [];
    
    end
    
    time_vec_SL  = [time_vec_SL_cell{:}];
    time_vec_TDT = [time_vec_TDT_cell{:}];    
    time_vec_SL  = [time_vec_SL numel(peaks_SL)-1];
    time_vec_TDT = [time_vec_TDT numel(peaks_SL)-1];    

    if debug
        figure;
        plot(time_vec_SL); hold on;
        plot(time_vec_TDT);
        legend({'SL', 'TDT'});
    end

    data_synced = [];

    for data_i = 1:size(data_tosync, 1)

        data = data_tosync(data_i, :);

        synced = interp1( ...
            time_vec_SL, ...
            data, ...
            time_vec_TDT, ...
            'linear', ...
            'extrap' ...
        );

        data_synced(data_i, :) = synced;
    end

    if debug
        sig_plot = data_synced(1,:);

        figure;
        plot(1:numel(sig_plot), sig_plot); hold on;
        xline(peaks_TDT, 'r');

        figure;
        plot(1:numel(sig_plot), sig_plot); hold on;
        xline(peaks_SL, 'r');
    end
end

