%% SCRIPT TO READ AND DISPLAY ONE-DIMENSIONAL WAVEFORM FROM FILE

clear all
close all
clc

% Figure print parameters
szAxScale = 20;
szAxLabel = 20;

%% FILE NAME
filePath = '/home/soundscape/Downloads'; %'/export/SoundScape/reduction/volumetric/tx';
fileName = 'x30';%'x61';
fileExt = 're';
slash = '/';
dot = '.';

arg = [filePath slash fileName dot fileExt];

%% READ MANIFOLD DATA
% Read data as a column vector

fid = fopen(arg);
data = fread(fid, inf, 'double', 0, 'b');
fclose(fid);

%% PLOT MANIFOLD
n = length(data);
fs = 1e6;

figure('name', 'Waveform from file')
plot((0:n-1)/fs,data);
title(['Waveform from file: ', fileName, dot, fileExt], 'fontsize', szAxLabel)
xlabel('Samples, n', 'fontsize', szAxLabel)
ylabel('Real amplitude', 'fontsize', szAxLabel)
axis 'square';
grid on
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);


figure
plot(1/(n/fs)*(-n/2:n/2-1), abs(fftshift(fft(data))))
axis 'square';
grid on
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
