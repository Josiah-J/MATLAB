%% MATCHED FILTERING
% This script performs matched filtering on real ultrasonic data.

%% OPTIONS
clear all
close all
clc

% Figure print parameters
szAxScale = 20;
szAxLabel = 20;
figRes = '-r300';

% Output path
outpath = '';

%% IMPORT THE TRANSMIT SIGNAL DATA
importfile('Tx 270kHz 0deg.log','tx')
tx = tx.';
ntx = length(tx);

%% IMPORT THE RECEIVE SIGNAL DATA
importfile('Rx 270kHz 0deg.log','rx');
rx = rx.';
nrx = length(rx);

if ntx < nrx
    tx = [tx zeros(1,nrx-ntx)];
    ntx = length(tx);
end;

%% SAMPLING/SIGNAL PARAMETERS

c = 1500;   % Speed of wave propagation
fc = 270e3; % Centre frequency
fs = 4e6;   % Sampling rate
dt = 1/fs;  % Sampling interval
n = nrx;    % Number of samples

t = dt*(0:n-1); % Time vector
df = 1/(n*dt);  % Frequency interval
f = df*(0:n-1); % Frequency vector

%% Plot transmit and receive signals

figure('name','Transmitted signal')
plot(t,tx);
title('Transmitted pulse', 'fontsize',szAxLabel)
xlabel('Time, t [sec]', 'fontsize',szAxLabel)
ylabel('Amplitude, V_t_x [V]', 'fontsize', szAxLabel)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'
print('-dpng', figRes, [outpath 'Transmitted signal'])

figure('name','Received signal')
plot(t,rx);
title('Received pulse', 'fontsize',szAxLabel)
xlabel('Time, t [sec]', 'fontsize',szAxLabel)
ylabel('Amplitude, V_r_x [V]', 'fontsize', szAxLabel)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'
print('-dpng', figRes, [outpath 'Received signal'])

%% MATCHED FILTERING

fs0b = fft(tx); % FFT of transmitted signal
fsb = fft(rx);  % FFT of received signal

fsmb = fsb.*conj(fs0b);                % Matched filtering
fsmb = fsmb.*(f <= f(n/2-1));          % Delete the negative frequency component
smb = ifft(fsmb).*exp(-1j*2*pi*fc*t);  % Basebanding

G = abs(smb);
figure('name', 'Signal after matched filtering')
plot(t,G);
title('Signal after matched filtering', 'fontsize',szAxLabel)
xlabel('Time, t [sec]', 'fontsize',szAxLabel)
ylabel('Amplitude, V_m_f [V]', 'fontsize', szAxLabel)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'
print('-dpng', figRes, [outpath 'Signal after matched filtering'])
