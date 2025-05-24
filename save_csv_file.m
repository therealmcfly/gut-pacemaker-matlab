% Prepare workspace
clear
clc
close all

load('ver_chdata_exp_16_output_512_ch1.csv');
data = readmatrix('ver_chdata_exp_16_output_512_ch1.csv');
clear
data = readmatrix('ver_chdata_exp_16_output_512_ch1.csv');
signal = downsample(data,16);
writematrix(signal, 'exp_16_output_32_ch1.csv');