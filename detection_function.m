function locs = detection_function(buffer, threshold, sampling_f)

% ------------------------- Low Pass Filter ----------------------------- %

padded_buffer = ones(1, length(buffer) + 120) * buffer(1);
padded_buffer(61:end-60) = buffer; % padded signal
padded_buffer(end-60:end) = buffer(end);

lowpass_signal = lowpass(padded_buffer,3,sampling_f); % filtering
lowpass_signal = lowpass_signal(61:end-60); % removing padding

% visualise_signal(lowpass_signal, sampling_f, "Low pass filtered signal butter");
% hfvt = fvtool(b, a, 'FrequencyScale', 'log'); % filter analysis tool

% ------------------------- High Pass Filter ---------------------------- %

importfile("fir_51.mat"); % Highpass kernel
load fir_51.mat fir_high_pass

% Filter info: w = 0.32Hz, windowed sinc FIR
padded_lowpass = ones(1, length(lowpass_signal) + 100) * lowpass_signal(1);
padded_lowpass(51:end - 50) = lowpass_signal;
padded_lowpass(end - 50:end) = lowpass_signal(end); % padded signal

highpass_signal = conv(padded_lowpass, fir_high_pass); % filtering
highpass_signal = highpass_signal(51:end - 50); % removing padding

% visualise_signal(highpass_signal, sampling_f, "High pass filtered signal");

% ------------------------- Artifact Removal ---------------------------- %

artifacts_removed = highpass_signal;
window_width = 100; 
ii = 1;
jj = ii + window_width;

art_size = 60;

while jj < length(artifacts_removed)
    window = artifacts_removed(ii:jj);
    loc = artifact_detect(window, threshold);

    if (isnan(loc) ~= 1)
    % artifact detected
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

% ----------------- Non Linear Energy (NEO) Transform ------------------- %

test_neo = NEO_transform(artifacts_removed);
% visualise_signal(test_neo, sampling_f, "NEO Signal");

% ---------------------- Moving Average Filter -------------------------- %

neo_filtered = moving_average_1s_window(test_neo, sampling_f);
% visualise_signal(neo_filtered, sampling_f, "Filtered NEO Signal");

% ---------------------- Edge Detection Signal -------------------------- %

kernel = -1 * [1 0 -1];
convolved_sig = conv(artifacts_removed, kernel, 'same');
% visualise_signal(convolved_sig, sampling_f, "Edge Detection Signal");

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
f_t_sig_base = mean(f_t_sig) * 59;
for z = 1:length(f_t_sig)
    if(f_t_sig(z) < f_t_sig_base)
        f_t_sig(z) = 0;
    end
end

% visualise_signal(f_t_sig, sampling_f, "Filtered Detection signal");

% ----------------------- Activation Detection -------------------------- %

% threshold calculation * scale
mad = mean(f_t_sig) * 5.9; 

% finding where detection signal > threshold
locs = find(f_t_sig > mad); 

% remove points in close proximity of each other (max 1 activation per window)
for i = 1:length(locs)
    locs(find(diff(locs) < 500, 1) + 1) = [];
end

% remove points past the window length
for i = 1:length(locs)
    if locs(i) > length(buffer)
        locs(i:numel(locs)) = [];
        break
    end
end

