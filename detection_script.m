% Prepare workspace
clear
clc
close all

% REALTIME_DATASET_MODE=1;

% BOARD_CONNECT=1;
board_ip = "192.168.137.223";

% Initial settings
bad_dataname = 'exp_16_output';
% good_dataname = 'pig41exp2';
channel = 2;

ar_threshold = 500;
buffer_size = 1000;
buffer_offset = buffer_size/2;
% buffer_offset = 100;
close_prox_act_remove_threshold = 500;

if exist("bad_dataname", "var")
    close_prox_act_remove_threshold = 500;
elseif exist("good_dataname", "var")
    close_prox_act_remove_threshold = 400;
end

% close_prox_act_remove_threshold = 500; % override threshold
ed_scalar_val = 59;
actd_scalar_val = 5.9;

% MODIFICATION = 0;

% Export data settings
% data_save = 1;
% ch_data_save = 1;
% ds_data_save = 1;
% lpf_data_save = 1;
% hpf_data_save = 1;
% ad_data_save = 1;
% neo_data_save = 1;
% maf_data_save = 1;
% ed_data_save = 1;
% actdpre_data_save = 1;
% actd_data_save = 1;

% Plot
whole_signal_plot = 1;
buffer_plot = 1;
buffer_num = [1,3];

lowpass_plot = 1; 
highpass_plot = 1;
ar_plot = 1;
neo_transform_plot = 1;
neo_maf_plot = 1;
ed_plot = 1;

plot_count = 1;
if exist("lowpass_plot", "var")
    plot_count = plot_count + 1;
end
if exist("highpass_plot", "var")
    plot_count = plot_count + 1;
end
if exist("ar_plot", "var")
    plot_count = plot_count + 1;
end
if exist("neo_maf_plot", "var")
    plot_count = plot_count + 1;
end
if exist("ed_plot", "var")
    plot_count = plot_count + 1;
end

% ----------------------------------------------------------------------- %
%                                CONFIG                                   %
% ----------------------------------------------------------------------- %

% --------------------------- Import data ------------------------------- %

time = 0;
lap = 0;

% if use_bad_signal_data > 0
if exist('bad_dataname', 'var')
    % -- Bad signal -- %
    dataname = bad_dataname;
    importfile([dataname, '.mat']);
    initial_sampling_f = 512;
    channel_select = channel; 
    signal = data_struct(1:end,channel_select);
    
    if exist('data_save', 'var') && data_save > 0
        % Saving verify file for C program 
        verify_filename = strjoin({dataname, '_', num2str(initial_sampling_f), '.csv'}, '');
        writematrix(data_struct, verify_filename);
        disp(['Saved verify file as: ', verify_filename]);
    end

    if exist('ch_data_save', 'var') && ch_data_save > 0
        % Save channel data to CSV
        channel_data = signal;
        ver_ch_filename = strjoin({'ver_chdata_', dataname, '_', num2str(initial_sampling_f), '_ch', num2str(channel_select), '.csv'}, '');
        writematrix(channel_data, ver_ch_filename);
        disp(['Downsampled signal saved as: ', ver_ch_filename]);
    end

     % visualise_signal(signal,sampling_f,"Original Data, Sampled At 512");
    signal = downsample(signal,16); % Downsample to 32Hz
    signal = signal.'; % transpose
    sampling_f = 32;
    % visualise_signal(signal,sampling_f,"Downsampled Data At 32Hz");
    
    if exist('ds_data_save', 'var') && ds_data_save > 0
        % Save downsampled signal to CSV
        ds_signal = signal.';
        signal = signal.'; % transpose
        ver_ds_filename = strjoin({'ver_ds_', dataname, '_', num2str(initial_sampling_f), '_ch', num2str(channel_select), '.csv'}, '');
        writematrix(ds_signal, ver_ds_filename);
        disp(['Downsampled signal saved as: ', ver_ds_filename]);
    end
    
elseif exist('good_dataname', 'var')
    % -- Good signal -- %
    dataname = good_dataname;
    importfile([dataname, '.mat']);
    initial_sampling_f = 30;
    sampling_f = 30;
    channel_select = channel;
    signal = sigRawSelect(channel_select,1:end);
    data_struct = sigRawSelect.';
    % % visualise_signal(signal,sampling_f,"Original Data At 512 Hz");
    % % Saving verify file for C program 
    % verify_filename = strjoin({dataname, '_', num2str(initial_sampling_f), '.csv'}, '');
    % writematrix(data_struct, verify_filename);
    % disp(['Saved verify file as: ', verify_filename]);
    
    if exist('data_save', 'var') && data_save > 0
        % Saving verify file for C program 
        verify_filename = strjoin({dataname, '_', num2str(initial_sampling_f), '.csv'}, '');
        writematrix(data_struct, verify_filename);
        disp(['Saved verify file as: ', verify_filename]);
    end
    if exist('ch_data_save', 'var') && ch_data_save > 0  
        % Save channel data to CSV

        channel_data = signal.';
        ver_ch_filename = strjoin({'ver_chdata_', dataname, '_', num2str(initial_sampling_f), '_ch', num2str(channel_select), '.csv'}, '');
        writematrix(channel_data, ver_ch_filename);
        disp(['Channel signal saved as: ', ver_ch_filename]);
    end
    if exist('ds_data_save', 'var') && ds_data_save > 0
        % Save downsampled signal to CSV
        ds_signal = signal.'; % transpose
        ver_ds_filename = strjoin({'ver_ds_', dataname, '_', num2str(initial_sampling_f), '_ch', num2str(channel_select), '.csv'}, '');
        writematrix(ds_signal, ver_ds_filename);
        disp(['Downsampled signal saved as: ', ver_ds_filename]);
    end
    signal = signal.';
else
    error("Error: Not specified good data or bad data.")
end

if exist('REALTIME_DATASET_MODE', 'var')
    Ts = 1/sampling_f
    filename = strjoin({dataname, '_', num2str(initial_sampling_f)}, '');
    mode = REALTIME_DATASET_MODE;

    % Connect to TCP server
    % if exist('BOARD_CONNECT', 'var')
    %     client = tcpclient(board_ip, 8081);
    % else
    %     client = tcpclient("127.0.0.1", 8081);
    % end
    % client = tcpclient("127.0.0.1", 8081);
    % 
    %  %Send mode identifier (1 byte)
    % write(client, uint8(REALTIME_DATASET_MODE), "uint8");
    % pause(0.1);  % Allow server time to process mode byte
    % 
    % % Send filename as uint8 character array
    % write(client, uint8(filename), "uint8");
    % pause(0.05);  % Allow server time to process filename
    % 
    % % Send sampling frequency (4 bytes)
    % write(client, int32(sampling_f), "int32");
    % pause(0.05);  % Allow server time to process sampling frequency
    % 
    % % Send channel number (4 bytes)
    % write(client, int32(channel), "int32");
    % pause(0.05);  % Allow server time to process channel number
    % 
    % fprintf("Waiting for go signal from server...\n");
    %   % Wait until relay server sends a start signal (1 byte, uint8)
    % while client.NumBytesAvailable < 1
    %     pause(0.01);  % prevent busy-waiting
    % end
    % 
    % start_signal = read(client, 1, "uint8");
    % 
    % if start_signal == 1
    %     disp('Start signal received. Beginning data transmission.');
    % else
    %     warning('Unexpected start signal received: %d', start_signal);
    % end

    relay_server(mode, filename, sampling_f, channel, signal);

    % Clean up
    fprintf("Finished sending signal from Realtime Dataset Mode. Exiting program.");
    % clear client;
    return;
end


locs = [];
i = 1;
j = i + buffer_size;
shift = 0;

% verification file variables
lpf_result = [];
hpf_result = [];
% ar_result = [];
% neo_result = [];

% plot vars
lpf_signal = [];
hpf_signal = [];
ar_signal = [];
neo_signal = [];
maf_signal = [];
ed_signal = [];

col = 1; % Column index for storing results


while j < length(signal)
    % read in signal input then transfer to buffer 
    buffer = signal(i:j);
    tic;

    % % ------------------ Visualize the Original Buffer ------------------ %
    if exist('buffer_plot', 'var') && buffer_plot > 0 && (shift+1 >= buffer_num(1) && shift+1 <= buffer_num(2))
        figure;
        set(gcf, 'Position', [100 100 400 700]);
        subplot(7,1,1); % First subplot for original signal
        plot(buffer);
        title(['Original Signal Before Low-Pass Filtering (Iteration ', num2str(i), ')']);
        xlabel('Samples');
        ylabel('Amplitude');
        grid on;
    end
    % figure;
    % subplot(2,1,1); % First subplot for original signal
    % plot(buffer);
    % title(['Original Signal Before Low-Pass Filtering (Iteration ', num2str(i), ')']);
    % xlabel('Samples');
    % ylabel('Amplitude');
    % grid on;
    
    % ------------------------------------------------------------------- %
    %                             PROCESSING                              %
    % ------------------------------------------------------------------- %
    
    % This stage involves preparing the data and filtering the signal. 
    % Firstly a high pass filter is used to remove baseline drift, then a 
    % lowpass filter to remove high frequency noise. A NEO transform is 
    % then used to accentuate the steep gradients which is then filtered 
    % with a moving average filter.
    
    % ------------------------- Low pass filter ------------------------- %
    padded_buffer = ones(1, length(buffer) + 120) * buffer(1);
    padded_buffer(61:end-60) = buffer; % padded signal
    padded_buffer(end-60:end) = buffer(end);

    lowpass_signal = lowpass(padded_buffer,3,sampling_f); % filtering
    lowpass_signal = lowpass_signal(61:end-60); % removing padding

   

    % % ------------------------- Low pass filter ------------------------- %
    % padded_buffer = ones(1, length(buffer) + 120) * buffer(1);
    % padded_buffer(61:end-60) = buffer; % padded signal
    % padded_buffer(end-60:end) = buffer(end);
    % 
    % Fs = 32;  % Sampling frequency
    % Wpass = 3; 
    % Wstop = 5.0410;
    % Ap = 0.1; % Passband ripple in dB
    % Ast = 60; % Stopband attenuation in dB
    % 
    % % ✅ Generate FIR lowpass filter manually (to get coefficients)
    % d = designfilt('lowpassfir', 'PassbandFrequency', Wpass, ...
    %            'StopbandFrequency', Wstop, 'SampleRate', Fs, ...
    %            'PassbandRipple', Ap, 'StopbandAttenuation', Ast);
    % 
    % b = d.Coefficients; % Get FIR filter coefficients
    % disp(b);
    % 
    % % ✅ Check if `D` is a valid digital filter object
    % if isa(d, 'digitalFilter')
    %     disp('Extracted FIR Coefficients:');
    %     disp(d.Coefficients);  % Only access if D is a valid filter
    %     writematrix(d.Coefficients, 'fir_coeffs_new.txt'); % Save for C implementation
    % else
    %     error('Filter design failed! `D` is not a valid digital filter object.');
    % end
    % 
    % % ✅ Apply actual lowpass filtering
    % lowpass_signal = lowpass(padded_buffer, 3, sampling_f);
    % lowpass_signal = lowpass_signal(61:end-60); % removing padding


    % % ------------------ Visualize the Low-Passed Signal ------------------ %
    if exist('buffer_plot', 'var') && buffer_plot > 0 && exist('lowpass_plot', 'var') && lowpass_plot > 0 && (shift+1 >= buffer_num(1) && shift+1 <= buffer_num(2))
        subplot(7,1,2); % Second subplot for low-pass filtered signal
        plot(lowpass_signal, 'r'); % Red color to differentiate
        title(['Low-Pass Filtered Signal (Iteration ', num2str(i), ')']);
        xlabel('Samples');
        ylabel('Amplitude');
        grid on;
    end
   
        
    % Store filtered result in a new column
    lpf_result(:, col) = lowpass_signal(:); 

    if(shift == 0)
        segment = lowpass_signal(1:buffer_size);
        lpf_signal = segment(:);  % Ensure column vector
    else
        segment = lowpass_signal((buffer_size - buffer_offset) + 1 : buffer_size);
        lpf_signal = [lpf_signal; segment(:)];
    end

    % % Save low pass filtered results to CSV file
    % ver_lpf_filename = strjoin({'ver_lpf_1stbuff', dataname, '_', num2str(initial_sampling_f), '_ch', num2str(channel_select), '.csv'}, '');
    % writematrix(lpf_result, ver_lpf_filename);
    % disp(['Downsampled signal saved as: ', ver_lpf_filename]);
    

    % ------------------------- High pass filter ------------------------ %
    % Filter info: w = 0.32Hz, windowed sinc FIR
    padded_lowpass = ones(1, length(lowpass_signal) + 100) * lowpass_signal(1);
    padded_lowpass(51:end - 50) = lowpass_signal;
    padded_lowpass(end - 50:end) = lowpass_signal(end); % padded signal

    importfile("fir_51.mat"); % Highpass kernel

    highpass_signal = conv(padded_lowpass, fir_high_pass); % filtering
    highpass_signal = highpass_signal(51:end - 50); % removing padding

    % % ------------------ Visualize the High-Passed Signal ------------------ %
    if exist('buffer_plot', 'var') && buffer_plot > 0 && exist('highpass_plot', 'var') && highpass_plot > 0 && (shift+1 >= buffer_num(1) && shift+1 <= buffer_num(2))
        % visualise_signal(highpass_signal, sampling_f, "High pass filtered signal");
        subplot(7,1,3); % Second subplot for low-pass filtered signal
        plot(highpass_signal, 'r'); % Red color to differentiate
        title(['High-Pass Filtered Signal (Iteration ', num2str(i), ')']);
        xlabel('Samples');
        ylabel('Amplitude');
        grid on;
    end

    % Store filtered result in a new column
    hpf_result(:, col) = highpass_signal(:); 

    hpf_seg = highpass_signal(26:end-25);
    if(shift == 0)
        segment = hpf_seg(1:buffer_size);
        hpf_signal = segment(:);  % Ensure column vector
    else
        segment = hpf_seg((buffer_size - buffer_offset) + 1 : buffer_size);
        hpf_signal = [hpf_signal; segment(:)];
    end

    % ------------------------- Artifact removal ------------------------ %

    artifacts_removed = highpass_signal;
    window_width = 100;
    ii = 1;
    jj = ii + window_width;

    art_size = 60;

    while jj < length(artifacts_removed)
        window = artifacts_removed(ii:jj);
        loc = artifact_detect(window, ar_threshold);

        if (isnan(loc) ~= 1)
        % artifact detected
        % disp(['Artifact detected at ', num2str((ii-1) + loc)]);
            if (loc + (art_size/2) <= length(window)) && (loc - (art_size/2) > 1)
            % window contains the whole artifact
                removed = artifact_remove(window, loc - art_size/2, loc + art_size/2);
                artifacts_removed(ii:jj) = removed;
            end
        else
        % artifact not detected or not fully in window
        % do nothing
        end

        ii = ii + window_width / 2;
        jj = ii + window_width; 
    end

    % visualise_signal(artifacts_removed,sampling_f,"art rem");
    if exist('buffer_plot', 'var') && buffer_plot > 0 && exist('ar_plot', 'var') && ar_plot > 0 && (shift+1 >= buffer_num(1) && shift+1 <= buffer_num(2))
        % visualise_signal(artifacts_removed,sampling_f,"art rem");

        subplot(7,1,4); % 4th subplot for Artifact Removed Signal
        plot(artifacts_removed, 'r'); % Red color to differentiate
        title(['Artifact Removed Signal (Iteration ', num2str(i), ')']);
        xlabel('Samples');
        ylabel('Amplitude');
        grid on;
    end


    % % Store artifact removed result in a new column
    ar_result(:, col) = artifacts_removed(:); 

    ar_seg = artifacts_removed(26:end-25);
    if(shift == 0)
        segment = ar_seg(1:buffer_size);
        ar_signal = segment(:);  % Ensure column vector
    else
        segment = ar_seg((buffer_size - buffer_offset) + 1 : buffer_size);
        ar_signal = [ar_signal; segment(:)];
    end

    % ----------------- Non Linear Energy (NEO) Transform --------------- %

    test_neo = NEO_transform(artifacts_removed);
    % visualise_signal(test_neo, sampling_f, "NEO Signal");
    neo_result(:, col) = test_neo(:);

    if exist('buffer_plot', 'var') && buffer_plot > 0 && exist('neo_transform_plot', 'var') && neo_maf_plot > 0 && (shift+1 >= buffer_num(1) && shift+1 <= buffer_num(2))
        % visualise_signal(neo_filtered, sampling_f, "Filtered NEO Signal");

        subplot(7,1,5); % 4th subplot for Filtered NEO Signal
        plot(test_neo, 'r'); % Red color to differentiate
        title(['NEO Transform Signal (Iteration ', num2str(i), ')']);
        xlabel('Samples');
        ylabel('Amplitude');
        grid on;
    end
     

    % ---------------------- Moving average filter ---------------------- %

    neo_filtered = moving_average_1s_window(test_neo, sampling_f);

    if exist('buffer_plot', 'var') && buffer_plot > 0 && exist('neo_maf_plot', 'var') && neo_maf_plot > 0 && (shift+1 >= buffer_num(1) && shift+1 <= buffer_num(2))
        % visualise_signal(neo_filtered, sampling_f, "Filtered NEO Signal");

        subplot(7,1,6); % 4th subplot for Filtered NEO Signal
        plot(neo_filtered, 'r'); % Red color to differentiate
        title(['Filtered NEO Signal (Iteration ', num2str(i), ')']);
        xlabel('Samples');
        ylabel('Amplitude');
        grid on;
    end

    maf_result(:, col) = neo_filtered(:);

    neo_seg = neo_filtered(26:end-25);
    if(shift == 0)
        segment = neo_seg(1:buffer_size);
        neo_signal = segment(:);  % Ensure column vector
    else
        segment = neo_seg((buffer_size - buffer_offset) + 1 : buffer_size);
        neo_signal = [neo_signal; segment(:)];
    end

    % ------------------------------------------------------------------- %
    %                          EDGE DETECTION                             %
    % ------------------------------------------------------------------- %

    % In this stage, the detection signal is created. To do this, the 
    % original signal is convolved with an edge detect kernel which gives a 
    % positive response when the input signal drops. The detection signal 
    % is then the element-wise product of the processed signal and the 
    % output of the edge detect kernel.

    kernel = -1 * [1 0 -1];
    convolved_sig = conv(artifacts_removed, kernel, 'same');
    % visualise_signal(convolved_sig, sampling_f, "Edge Detection");

    % convolve edge detection and neo processed signals
    f_t_sig = convolved_sig(2:end) .* neo_filtered;

    % remove negative values
    for z = 1:length(f_t_sig)
        if(f_t_sig(z) < 0)
            f_t_sig(z) = 0;
        end
    end

    % square values of detection signal to pronounce higher peaks
    f_t_sig = f_t_sig.^2;

    % remove noise with mean of the signal * scalar value as baseline 
    f_t_sig_base = mean(f_t_sig) * ed_scalar_val;
    
    % f_t_sig_base = mean(f_t_sig);
    for z = 1:length(f_t_sig)
        if(f_t_sig(z) < f_t_sig_base)
            f_t_sig(z) = 0;
        end
    end
    
    if exist('buffer_plot', 'var') && buffer_plot > 0 && exist('ed_plot', 'var') && ed_plot > 0 && (shift+1 >= buffer_num(1) && shift+1 <= buffer_num(2))
        subplot(7,1,7); % 4th subplot for Filtered NEO Signal
        plot(f_t_sig, 'r'); % Red color to differentiate
        title(['Edge Detection Signal (Iteration ', num2str(i), ')']);
        xlabel('Samples');
        ylabel('Amplitude');
        grid on;
    end

    ed_result(:,col) = f_t_sig(:);

    ed_seg = f_t_sig(26:end-25);
    if(shift == 0)
        segment = ed_seg(1:buffer_size);
        ed_signal = segment(:);  % Ensure column vector
    else
        segment = ed_seg((buffer_size - buffer_offset) + 1 : buffer_size);
        ed_signal = [ed_signal; segment(:)];
    end

    
    
% visualise_signal(f_t_sig, sampling_f, "Detection signal");

    % ------------------------------------------------------------------- %
    %                        ACTIVATION DETECTION                         %
    % ------------------------------------------------------------------- %

    % The activations are detected using a variable threshold method 
    % performed on the output of the two previous processing stages. The 
    % times at which the detection signal exceeds the threshold are marked 
    % as potential activation times. 

    % threshold calculation * scale
    
    mad = mean(f_t_sig) * actd_scalar_val;

    % finding where detection signal > threshold

    offset = buffer_offset * shift;
    find_result = find(f_t_sig > mad);
    loc = find_result + offset;

    if ~isempty(loc)
        fprintf("Appending %d activation(s): %s\n", length(loc), mat2str(loc));
    end
    
    % append to locs array
    
    locs = [locs loc];
    
    % calculate window shift
    shift = shift + 1;  
    i = i + buffer_offset;
    j = i + buffer_size;
    col = col + 1;
    
    lap = lap+1;
    time = time + toc;

end

if exist('lpf_data_save', 'var') && lpf_data_save > 0
    % Save low pass filtered results to CSV file
    ver_lpf_filename = strjoin({'ver_lpf_', dataname, '_', num2str(initial_sampling_f),'_ch', num2str(channel_select), '.csv'}, '');
    writematrix(lpf_result, ver_lpf_filename);
    disp(['Lowpass filterd signal saved as: ', ver_lpf_filename]);
end
if exist('hpf_data_save', 'var') && hpf_data_save > 0
    % Save high pass filtered results to CSV file
    ver_hpf_filename = strjoin({'ver_hpf_', dataname, '_', num2str(initial_sampling_f),'_ch', num2str(channel_select), '.csv'}, '');
    writematrix(hpf_result, ver_hpf_filename);
    disp(['High pass filtered signal saved as: ', ver_hpf_filename]);
end
if exist('ad_data_save', 'var') && ad_data_save > 0
    % Save artifact removed results to CSV file
    ver_ar_filename = strjoin({'ver_ad_', dataname, '_', num2str(initial_sampling_f), '_ch', num2str(channel_select), '.csv'}, '');
    writematrix(ar_result, ver_ar_filename);
    disp(['artifact removed signal saved as: ', ver_ar_filename]);
end
if exist('neo_data_save', 'var') && neo_data_save > 0
    % Save neo transformed results to CSV file
    ver_neo_filename = strjoin({'ver_neo_', dataname, '_', num2str(initial_sampling_f), '_ch', num2str(channel_select), '.csv'}, '');
    writematrix(neo_result, ver_neo_filename);
    disp(['neo transformed signal saved as: ', ver_neo_filename]);
end 
if exist('maf_data_save', 'var') && maf_data_save > 0
    % Save moving average filtered results to CSV file
    ver_maf_filename = strjoin({'ver_maf_', dataname, '_', num2str(initial_sampling_f), '_ch', num2str(channel_select), '.csv'}, '');
    writematrix(maf_result, ver_maf_filename);
    disp(['moving average filtered signal saved as: ', ver_maf_filename]);
end
if exist('ed_data_save', 'var') && ed_data_save > 0
    % Save edge detection results to CSV file
    ver_ed_filename = strjoin({'ver_ed_', dataname, '_', num2str(initial_sampling_f), '_ch', num2str(channel_select), '.csv'}, '');
    writematrix(ed_result, ver_ed_filename);
    disp(['edge detection signal saved as: ', ver_ed_filename]);
end
if exist('actdpre_data_save', 'var') && actd_data_save > 0
    % Save activation detection results to CSV file
    actdpre_result(:, 1) = locs(:);
    ver_actdpre_filename = strjoin({'ver_actdpre_', dataname, '_', num2str(initial_sampling_f), '_ch', num2str(channel_select), '.csv'}, '');
    writematrix(actdpre_result, ver_actdpre_filename);
    disp(['pre activation detections saved as: ', ver_actdpre_filename]);
end

% remove points in close proximity of each other (max 1 activation per window)
for i = 1:length(locs)
    activation_prox_arr = diff(locs);
    close_act_idx = find(activation_prox_arr < close_prox_act_remove_threshold, 1);
    if ~isempty(close_act_idx)
        fprintf("Removing close prox activation : %d\n", locs(close_act_idx + 1));
        locs(close_act_idx + 1) = [];
    end

end


% remove points past the window length
for i = 1:length(locs)
    fprintf("%d\n", locs(i));
    if locs(i) > length(signal)
        % fprintf("%d\n", locs(i));
        % fprintf("%d\n", length(signal));
        % locs(i:numel(locs)) = [];
        error('activation detection passed the window length');    
        break
    end
end

if exist('actd_data_save', 'var') && actd_data_save > 0
    % Save activation detection results to CSV file
    actd_result(:, 1) = locs(:);
    ver_actd_filename = strjoin({'ver_actd_', dataname, '_', num2str(initial_sampling_f), '_ch', num2str(channel_select), '.csv'}, '');
    writematrix(actd_result, ver_actd_filename);
    disp(['final activation detections saved as: ', ver_actd_filename]);
end

% ----------------------------------------------------------------------- %
%                               PLOTTING                                  %
% ----------------------------------------------------------------------- %

if exist('whole_signal_plot', 'var')
    figure;
    set(gcf, 'Position', [100 100 1200 750]);

    curr_count = 1;

    subplot(plot_count,1,curr_count);
    plot(signal, 'r');
    title('Original Signal');
    for idx = 1:length(locs)
        xline(locs(idx), 'b--');
    end
    
    if(lowpass_plot > 0)
        curr_count = curr_count +1;
        subplot(plot_count,1,curr_count);
        plot(lpf_signal, 'r');
        title('Low-Pass Filtered Signal');
        for idx = 1:length(locs)
            xline(locs(idx), 'b--');
        end
    end

    if exist("highpass_plot","var")
        curr_count = curr_count +1;
        subplot(plot_count,1,curr_count);
        plot(hpf_signal, 'r');
        title('High-Pass Filtered Signal');
        for idx = 1:length(locs)
            xline(locs(idx), 'b--');
        end
    end

    if exist("ar_plot", "var" )
        curr_count = curr_count +1;
        subplot(plot_count,1,curr_count);
        plot(ar_signal, 'r');
        title('Artifact Removed Signal');
        for idx = 1:length(locs)
            xline(locs(idx), 'b--');
        end
    end
    
    if exist("neo_maf_plot", "var")
        curr_count = curr_count +1;
        subplot(plot_count,1,curr_count);
        plot(neo_signal, 'r');
        title('NEO Transformed Signal');
        for idx = 1:length(locs)
            xline(locs(idx), 'b--');
        end
    end

    if exist("ed_plot", "var")
        curr_count = curr_count + 1;
        subplot(plot_count,1,curr_count);
        plot(ed_signal, 'r');
        title('Edge Detection Signal');
        for idx = 1:length(locs)
            xline(locs(idx), 'b--');
        end
    end
end


padded_signal = ones(1, length(signal) + 120) * signal(1);
padded_signal(61:end-60) = signal; % padded signal
padded_signal(end-60:end) = signal(end);
lowpass_signal = lowpass(padded_signal,3,sampling_f); % filtering
lowpass_signal = lowpass_signal(61:end-60); % removing padding

T = 1 / sampling_f;
N = length(lowpass_signal);
t = 0:T:T * (N - 1);
figure
plot(t,lowpass_signal);
hold on
plot(t(locs),lowpass_signal(locs),'r*');
set(gca,'XMinorTick','on');
title("Activation points of signal");
xlabel("Time (s)");
ylabel("Signal, mV") 

disp(time/lap);

