%% SCRIPT TO READ AND DISPLAY MANIFOLD DATA FROM FILE

clear all
close all
clc

% Figure print parameters
szAxScale = 20;
szAxLabel = 20;

%% FILE NAME
filePath = '/export/SoundScape/probe/planar/rx';
fileName = 'x';
fileExt = 're';

arg = [filePath '/' fileName '.' fileExt];

%% NUMBER OF RECEIVE ELEMENTS

N = 64;

%% READ MANIFOLD DATA
% Read data as a column vector

fid = fopen(arg);
data = fread(fid, inf, 'double');
fclose(fid);

%% RESHAPE DATA
% Reshape data to form a matrix with a given dimension

rowMajor = true;

if rowMajor == true;
    data = reshape(data, [], N);
    data = data.';
else
    data = reshape(data, N, []);
end;

%% PLOT MANIFOLD

figure('name', 'Manifold from file')
colormap(jet(256));
imagesc(data);
title('Manifold', 'fontsize', szAxLabel)
xlabel('Samples, n', 'fontsize', szAxLabel)
ylabel('Receiver elements, m', 'fontsize', szAxLabel)
axis 'square';
axis 'xy';
colorbar
grid on
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
