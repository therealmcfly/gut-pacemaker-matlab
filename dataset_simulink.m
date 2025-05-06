% Prepare workspace
clear
clc
close all

filename = 'exp_16_output_32_ch1.csv';
Fs = 32;
Ts = 1 / Fs;

% Read signal data
data = readmatrix(filename);
num_samples = length(data);

% Connect to TCP server
t = tcpclient("127.0.0.1", 8080);

% Send filename (as uint8 char array)
write(t, uint8(filename), "uint8");
pause(0.1);  % Allow server to receive filename

% Send data at ~32 Hz
% for i = 1:length(data)
for i = 1:20  % or use num_samples
    tic;
    write(t, data(i), "double");
    while toc < Ts
        % active wait (spins CPU but tighter timing)
    end
end

% Clean up
clear t;
