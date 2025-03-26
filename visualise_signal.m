function visualise_signal(data,sampling_f,plot_title)
% VISUALISE_SIGNAL plots the data on a graph with the given sample_f and
% applies all required formatting such as title, labels, etc.
% INPUT: data, the signal to plot
%        sample_f, the sampling frequency in Hz
%        title, a string containing the desired title for the figure
% Author: Larissa Marr

T = 1 / sampling_f;
N = length(data);
    
% Create time vector
% t = 0:T:T * (N - 1);
    
% Plot and format
figure;
% plot(t, data);
plot(data);
title(plot_title);
% xlabel("Time (seconds)");
xlabel("Samples");
ylabel("Signal, mV") % TO DO, check unit is correct

end

