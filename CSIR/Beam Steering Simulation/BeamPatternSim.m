%% PHASED-TRANSDUCER ARRAY BEAM PATTERN SIMULATION
% This code simulates the beam pattern of a phased ultrasound array given
% the number of elements, the required steering angle.
% A effect of phase and amplitude errors on the steering angle is
% demonstrated.

%% OPTIONS
close all
clear all
clc

% Figure print parameters
szAxScale = 20;
szAxLabel = 20;
szAxLabel2 = 14;
figRes = '-r300';

%% PULSE PARAMETERS
c = 1500; % Speed of sound in water
fc = 420e3; % Centre frequency in Hz
B = 20e3; % Signal bandwidth
B0 = B/2; % Baseband bandwidth
Tp = 2e-3; % Chirp pulse width
K = B0/Tp; % Chirp rate

beta = fc - B0;
lambda_c = c/fc; % Wavelength at the centre frequency
lambda_max = c/(fc - B0); % Maximum wavelength
lambda_min = c/(fc + B0); % Minimum wavelength

%% SAMPLING PARAMETERS
% fs = 4*(fc+B0); % Temporal sampling rate
fs = 10*B0;
dt = 1/fs; % Temporal sampling interval

%% ARRAY PARAMETERS

N = 100; % Number of transmit elements in the array
d = 12.5e-3;%25e-3;%10.6e-3;%0.493*lambda_c; % Inter-element spacing
a = 3.8e-3; % Element width
D = d*(N-1)+a; % Aperture size

theta_s = 0; % Required steering angle
delta_tau = (d/c)*sind(theta_s); % Inter-element time delay for given steering angle
tau = delta_tau*(0:N-1); % Time delay vector

%% PHASE ERRORS

pha_error_type = 'selected-random';
phaErrorInt = [-pi/18 pi/18];
nFaultyElms = ceil(1*N); % Number of elements with phase errors

a = phaErrorInt(1);
b = phaErrorInt(2);

switch pha_error_type
    case 'none'
        pha_error = zeros(1,N); % No phase error on array
    case 'random'
        rng('shuffle');
        pha_error = a + (b-a)*rand(1,N); % Random phase errors within the specified interval
    case 'selected-random'
        rng('shuffle');
        faulty_elm = randi([1,N], nFaultyElms, 1);
        pherror = a + (b-a)*rand(1,nFaultyElms); % Random phase errors within the specified interval
        pha_error = zeros(1,N);
        pha_error(faulty_elm) = pherror;
    case 'equal'
        rng('shuffle');
        pha_error = (a + (b-a)*rand(1,1))*ones(1,N); % Random phase errors within the specified interval        
    otherwise
        pha_error = zeros(1,N); % No phase error on array
end;

figure('name','Phase error')
plot(1:N, (180/pi)*(pha_error));
xlabel('Element', 'fontsize', szAxLabel)
ylabel('Phase, \theta [deg]', 'fontsize', szAxLabel)
ylim([-180 180]); xlim([1 N]);
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'

%% AMPLITUDE ERRORS
amp_error_type = 'random';
ampErrorInt = [0.7 1];

switch amp_error_type 
    case 'none'
        amp_error = ones(1,N); % No amplitude error on each transmit element
    case 'random'
        rng('shuffle');
        a = ampErrorInt(1);
        b = ampErrorInt(2);
        amp_error = a+(b-a)*rand(1,N); % Random amplitude error on transmit elements
    otherwise
        amp_error = ones(1,N); % No amplitude error on each transmit element
end;

figure('name','Amplitude error')
plot(1:N, amp_error);
xlabel('Element', 'fontsize', szAxLabel)
ylabel('Amplitude', 'fontsize', szAxLabel)
ylim([0 2]); xlim([1 N]);
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'

figure('name','Amplitude error in dB')
plot(1:N, 20*log10(amp_error));
xlabel('Element', 'fontsize', szAxLabel)
ylabel('Amplitude', 'fontsize', szAxLabel)
xlim([1 N]);
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'


%% SHADING FUNCTION
shading_function = 'rect';

switch shading_function
    case 'rect'
        w = ones(1,N);
    case 'hamming';
        w = hamming(N);
    case 'hann'
        w = hann(N);
    case 'blackman-harris'
        w = window(@blackmanharris,N);
    case 'triangle'
        w = window(@triang, N);
    otherwise
        w = ones(1,N);
end;

figure('name','Shading function')
plot(1:N, w);
xlabel('Element', 'fontsize', szAxLabel)
ylabel('Amplitude', 'fontsize', szAxLabel)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square', axis 'tight'

%% SIMULATION PARAMETERS/VECTORS
rFar = D.^2/lambda_c;
r = 400; % Simulation range
dtheta = 0.04; % Angular acquisition interval in degrees
theta_start = -3;
theta_end = 3;
theta_scan = theta_end - theta_start;
m = ceil(theta_scan/dtheta); % Number of angular bins
vecTheta = theta_start + dtheta*(0:m-1); % Angular vector

Ts = (r/1.1)/c; % Sampling start time
Te = r/c + Tp*1.3; % Sampling end time
n = 2*ceil((Te - Ts)/(2*dt));
t = Ts + dt*(0:n-1); % Time vector

df = 1/(n*dt);
f = df*(0:n-1); % Frequency vector

%% SIMULATION

%% The transmitted signal
td = t-Ts;
pha = 2*pi*(beta*td + K*td.^2);
pt = exp(1j*pha).*(td>=0 & td <= Tp);

% figure('name','The transmitted pulse in the time domain')
% plot(td,real(pt))
% title('Transmitted pulse', 'fontsize', szAxLabel)
% xlabel('Time, t [s]', 'fontsize', szAxLabel)
% ylabel('Amplitude (Real), p(t)', 'fontsize', szAxLabel)
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);
% axis 'square'
% 
% figure('name','The magnitude spectrum of transmitted pulse')
% plot(f, abs((fft(pt))))
% title('Magnitude spectrum of transmitted pulse', 'fontsize', szAxLabel)
% xlabel('Frequency, f [Hz]', 'fontsize', szAxLabel)
% ylabel('Magnitude, |P(f)|', 'fontsize', szAxLabel)
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);
% axis 'square'; axis tight
% 
pt = pt.*exp(-1j*2*pi*fc*t);

% figure('name','The baseband transmitted pulse in the time domain')
% plot(td,real(pt))
% title('Basebanded transmitted pulse', 'fontsize', szAxLabel)
% xlabel('Time, t [s]', 'fontsize', szAxLabel)
% ylabel('Amplitude (Real), p(t)', 'fontsize', szAxLabel)
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);
% axis 'square'
% 
% figure('name','The magnitude spectrum of transmitted pulse')
% plot(df*(-n/2:n/2-1), abs(fftshift((fft(pt)))))
% title('Magnitude spectrum of basebanded transmitted pulse', 'fontsize', szAxLabel)
% xlabel('Frequency, f [Hz]', 'fontsize', szAxLabel)
% ylabel('Magnitude, |P(f)|', 'fontsize', szAxLabel)
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);
% axis 'square'; axis tight

%% The received signal

st = zeros(m,n);

for i = 1:m
    i
    for j = -N/2:N/2-1
        R = (r.^2+(abs(j)*d).^2 - 2*r*(j*d)*cosd(90 - vecTheta(i))).^0.5;
        
        td = t - tau(j+N/2+1) - R/c;
        pha = 2*pi*(beta*td + K*td.^2);
        
        st(i,:) = st(i,:) + (amp_error(j+N/2+1) * exp(1j*pha_error(j+N/2+1)) * w(j+N/2+1) * exp(1j*pha) .*(td >= 0 & td <= Tp) );
    end;
    
    st(i,:) = st(i,:).*exp(-1j*2*pi*fc*t);
end;

theta_slice = theta_s;
theta_bin = ceil((theta_slice - theta_start)/dtheta)+1;

% figure
% plot( t, real(st(theta_bin,:)) );
% title(['Received signal at \theta = ',num2str(theta_slice),'^o'], 'fontsize', szAxLabel)
% xlabel('Time, t [s]', 'fontsize', szAxLabel)
% ylabel('Amplitude', 'fontsize', szAxLabel)
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);
% axis 'square'

% figure('name','The magnitude spectrum of received signal at theta = 0')
% plot(f, abs((fft(st(theta_bin,:)))))
% title(['Magnitude spectrum of received signal at \theta = ', num2str(theta_slice), '^o'], 'fontsize', szAxLabel)
% xlabel('Frequency, f [Hz]', 'fontsize', szAxLabel)
% ylabel('Magnitude, |S(f)|', 'fontsize', szAxLabel)
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);
% axis 'square'; axis tight

G = real(st);
figure('name','Acquired ultrasound data')
colormap(hot)
imagesc(c*t,vecTheta,G);
axis 'square'; axis 'xy'; colorbar;
xlabel('Time, t [s]', 'fontsize', szAxLabel)
ylabel('Angle, \theta [deg]', 'fontsize', szAxLabel)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);

%% BEAM PATTERN PROCESSING USING MATCHED FILTER

f = df*(-n/2:n/2-1);

s0b = pt;
fs0b = fty(s0b).*exp(-1j*2*pi*f*(r/c - Ts));

fsb = fty(st);
fsmb = fsb.*(ones(m,1)*conj(fs0b)); % Matched filtering and windowing
fsmb = fsmb.*(ones(m,1)*(abs(f) <= B0));
smb = ifty(fsmb);

tm=(r/c)+dt*(-n/2:n/2-1);   % fast-time array after matched filtering

G=abs(smb);
xg=max(max(G)); ng=min(min(G)); cg=255/(xg-ng);

figure('name', 'SAS signal after range compression')
colormap(hot(256))
imagesc(c*tm, vecTheta, G);
title('Matched filtered signal variation with angle', 'fontsize', szAxLabel)
xlabel('Range, ct [m]', 'fontsize', szAxLabel)
ylabel('Angle, \theta [deg]','fontsize', szAxLabel)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis('square'); axis('xy'); colorbar

figure('name','Received signal at boresight after matched filtering');
plot(c*tm, G(theta_bin,:))
title(['Magnitude of match filtered signal at \theta = ', num2str(theta_slice), '^o'], 'fontsize', szAxLabel)
xlabel('Range, ct [m]', 'fontsize', szAxLabel)
ylabel('Amplitude, s_m(t)', 'fontsize', szAxLabel)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'; axis tight

%% DIG OUT VALUES OF THE PROCESSED DATA FROM AT THE IMAGING RANGE

act_Beam = zeros(1,m);
for i = 1:m
%     temp = abs(smb(i,:));
%     temp1 = (temp == max(temp));
%     idx = find(temp1,1,'first');
%     act_Beam(i) = abs(smb(i,idx));
    act_Beam(i) = abs(smb(i,n/2+1));
end;
act_Beam_normalised = act_Beam./max(act_Beam);

% figure('name','Measured radiation pattern (extracted from processed data)')
% plot(vecTheta, act_Beam_normalised)
% title('Measured radiation pattern on linear scale', 'fontsize', szAxLabel)
% xlabel('Angle, \theta [deg]', 'fontsize', szAxLabel)
% ylabel('Normalized level', 'fontsize', szAxLabel)
% axis 'square'; axis 'tight'
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);

figure('name','Measured radiation pattern - dB scale')
plot(vecTheta, 20*log10(act_Beam_normalised))
title(['Radiation pattern with r = ',num2str(r),'m (',num2str(rFar,4),'m); amp-error = ', num2str(20*log10(ampErrorInt(1)),1), 'dB;', char(10),'phase-error = +-', num2str(180/pi*phaErrorInt(2)), '^o; window = ', shading_function ], 'fontsize', szAxLabel2)
xlabel('Angle, \theta [deg]', 'fontsize', szAxLabel)
ylabel('Normalized level, [dB]', 'fontsize', szAxLabel)
axis 'square'; axis 'tight'
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
print('-dpdf', figRes, './Radiation pattern 1')
print('-dpng', figRes, './Radiation pattern 1')

figure('name','Measured radiation pattern - Polar plot')
polar(vecTheta*pi/180, act_Beam_normalised)
title(['Measured radiation pattern'], 'fontsize', szAxLabel)
axis 'square'; %axis tight
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);