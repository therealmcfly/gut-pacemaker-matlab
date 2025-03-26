function output = moving_average_1s_window(input_signal, sample_rate)
% MOVING_AVERAGE_1S_WINDOW applies a moving average filter to the input
% signal, using a window of width 1 second.
%
% INPUT:    input_signal, the input signal as a 1xN row vector
%           sample_rate, the sample rate of the input signal in Hz
% OUTPUT:   output_signal, the filtered signal as a 1xN row vector

time = 1; % second
N = length(input_signal);
M = sample_rate * time; % size of window as number of points

% prepare vectors
sum = zeros(1, N);
output = zeros(1, N);
    
% calculate moving average
for i = 1:(N - M)
    for j = 1:M
        
        sum(i) = sum(i) + input_signal(i + j);
    end
    output(i) = sum(i) / M;
end
end

