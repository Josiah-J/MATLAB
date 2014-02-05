%% PULSE TRAIN PICKER SIMULATOR

%% OPTIONS
close all
clear all
clc

% Figure print parameters
szAxScale = 20;
szAxLabel = 20;

%% WAVEFORM CONFIGURATION

fc = 420e3; % Centre frequency
Tp = 50e-6; % Pulse width
Np = 4; % Number of pulses

dTp = 2.0*Tp; % Interpulse period
PRF = 1/dTp; % Pulse repetition frequency

%% TEMPORAL SAMPLING CONFIGURATION

fs = 20*fc; % Temporal sampling frequency
dt = 1/fs; % Temporal sampling interval

%% TEMPORAL SAMPLING VECTOR

Ts = 0;
Tf = (Np*dTp);

n = 2*ceil((Tf - Ts)/(2*dt));

t = Ts+dt*(0:n-1);

df = 1/(n*dt);
f = df*(-n/2:n/2-1);

%% THE PULSE TRAIN

pt = zeros(1,n);

for i = 1:Np
    tdi = t - (i-1)*dTp;
    pha = 2*pi*fc*tdi;
    temp = exp(1j*pha).*(tdi >= 0 & tdi <= Tp);
    
    pt = pt + temp;
end;

figure('name', 'The pulse train')
plot(t - Ts, real(pt))
title('The pulse train', 'fontsize', szAxLabel)
xlabel('Time, t [s]', 'fontsize', szAxLabel)
ylabel('Amplitude (Real), p(t)', 'fontsize', szAxLabel)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'; axis tight;

%% FOURIER SPECTRUM

fpt = fft(pt);
fpt = fftshift(fpt);

figure('name', 'Pulse train spectrum')
plot(f, abs(fpt));
xlabel('Frequency, f [Hz]', 'fontsize', szAxLabel)
ylabel('Magnitude, |P(f)|', 'fontsize', szAxLabel)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'; axis tight;


%% PULSE PICKER

% Multiply pulse train by a time-dependent phase

for i = 0:Np-1
    
    tdi = dTp*i;
    
    pha = 2*pi*f*(tdi - 0.25*dTp);
    ftemp = fpt.*exp(1j*pha);
    temp = ifft(fftshift(ftemp));
    
    pti = temp(t-Ts >= 0 & t-Ts <= dTp);
    
    figure('name', ['Pulse ', num2str(i+1)])
    plot(t(t-Ts >= 0 & t-Ts <= dTp) - Ts, real(pti))
    title(['Pulse ', num2str(i+1)], 'fontsize', szAxLabel)
    xlabel('Time, t [s]', 'fontsize', szAxLabel)
    ylabel('Amplitude (Real), p(t)', 'fontsize', szAxLabel)
    h_fig=get(gcf,'CurrentAxes');
    set(h_fig, 'fontsize', szAxScale);
    axis 'square'; axis tight;
end;


