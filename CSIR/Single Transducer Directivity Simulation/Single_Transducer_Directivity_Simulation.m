%% TRANSDUCER DIRECTIVITY SIMULATION
% The script simulates the acquisition and processing of ultrasonic signals
% to calibrate the directivity of a single sonar transducer. A sonar
% transducer consists of 6 ceramic elements. This simulation investigates
% the effect of varying the amplitude weighting of individual elements on
% the directivity characteristics namely, the beamwidth and sidelobe levels


%% GENERAL OPTIONS
close all
clear all
clc

% Figure print options
szAxScale = 20;
szAxLabel = 20;
szAxLabel2 = 14;
figRes = '-r300';

% File output options
outpath = './';
filename = 'Radiation pattern';

% Beam pattern interpolation
beamInterp = 1;
interpFactor = 100;

% Near-field compensation option
nFieldComp = 0;

%% PULSE CONFIGURATION
c = 1500; % Speed of sound in water
fc = 420e3; %200e3;% Centre frequency in Hz
B = 5e3; % Signal bandwidth
B0 = B/2; % Baseband bandwidth
Tp = 4e-3; % Chirp pulse width
K = B0/Tp; % Chirp rate

beta = fc - B0;
lambda_c = c/fc; % Wavelength at the centre frequency
lambda_max = c/(fc - B0); % Maximum wavelength
lambda_min = c/(fc + B0); % Minimum wavelength

%% SAMPLING CONFIGURATION

% fs = 4*(fc+B0); % Temporal sampling rate
fs = 5*B;
dt = 1/fs; % Temporal sampling interval

%% ARRAY CONFIGURATION

N = 6; %100; % Number of elements in the array
d = 3.0e-3 + 0.1e-3;  %12.5e-3;  %20e-3;%10.6e-3;% 2*lambda_c; %Inter-element spacing
a = 3.0e-3;%d; % Element width
D = d*(N-1)+a; % Aperture size
phi_d = asind(0.88*lambda_c/D); % Divergence angle of a rectangular aperture

thetaScan_min = 2*asind(lambda_min/(2*d)); % Maximum scanning angle at highest frequency
thetaScan_c = 2*asind(lambda_c/(2*d)); % Maximum scanning angle at centre frequency
thetaScan_max = 2*asind(lambda_max/(2*d)); % Maximum scanning angle at lowest frequency
theta_s = 0; % Required steering angle

delta_tau = (d/c)*sind(theta_s); % Inter-element time delay for given steering angle
tau = delta_tau*(0:N-1); % Time delay vector

%% Array phase errors

pha_error_type = 'none';
phaErrorInt = [-pi/12 pi/12];

if strcmp(pha_error_type, 'none')
       phaErrorInt = [0 0];
end;
nFaultyElms = ceil(0.15*N); % Number of elements with phase errors

b1 = phaErrorInt(1);
b2 = phaErrorInt(2);

switch pha_error_type
    case 'none'
        pha_error = zeros(1,N); % No phase error on array
    case 'random'
        rng('shuffle');
        pha_error = b1 + (b2-b1)*rand(1,N); % Random phase errors within the specified interval
    case 'selected-random'
        rng('shuffle');
        faulty_elm = randi([1,N], nFaultyElms, 1);
        pherror = b1 + (b2-b1)*rand(1,nFaultyElms); % Random phase errors within the specified interval
        pha_error = zeros(1,N);
        pha_error(faulty_elm) = pherror;
    case 'equal'
        rng('shuffle');
        pha_error = (b1 + (b2-b1)*rand(1,1))*ones(1,N); % Random phase errors within the specified interval        
    otherwise
        pha_error = zeros(1,N); % No phase error on array
end;

figure('name','Phase error')
stem(1:N, (180/pi)*(pha_error));
title('Phase error of transducer array elements', 'fontsize', szAxLabel2)
xlabel('Elements, n', 'fontsize', szAxLabel)
ylabel('Phase, \theta_e [deg]', 'fontsize', szAxLabel)
ylim([-180 180]); xlim([1 N]);
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'

%% Array amplitude errors

amp_error_type = 'none';
ampErrorInt = [0.5 1];
b1 = ampErrorInt(1);
b2 = ampErrorInt(2);
nFaultyElms = ceil(0.15*N); % Number of elements with amplitude errors

if strcmp(amp_error_type, 'none')
       ampErrorInt = [1, 1];
end;
switch amp_error_type 
    case 'none'
        amp_error = ones(1,N); % No amplitude error on each transmit element
    case 'random'
        rng('shuffle');        
        amp_error = b1+(b2-b1)*rand(1,N); % Random amplitude error on transmit elements
    case 'selected-random'
        rng('shuffle');
        faulty_elm = randi([1,N], nFaultyElms, 1);
        aperror = b1 + (b2-b1)*rand(1,nFaultyElms); % Random amplitude errors within the specified interval
        amp_error = ones(1,N);
        amp_error(faulty_elm) = aperror;
    otherwise
        amp_error = ones(1,N); % No amplitude error on each transmit element
end;

% figure('name','Amplitude error')
% stem(1:N, amp_error);
% title('Amplitude error of transducer array elements', 'fontsize', szAxLabel2)
% xlabel('Elements, n', 'fontsize', szAxLabel)
% ylabel('Amplitude', 'fontsize', szAxLabel)
% ylim([0 2]); xlim([1 N]);
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);
% axis 'square'

figure('name','Amplitude error in dB')
stem(1:N, 20*log10(amp_error));
title('Amplitude error of transducer array elements', 'fontsize', szAxLabel2)
xlabel('Elements, n', 'fontsize', szAxLabel)
ylabel('Amplitude, [dB]', 'fontsize', szAxLabel)
xlim([1 N]);
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'

%% Array shading function
shading_function = 'rect';

% Optimum combination: [0.3, 0.68, 1.0, 1.0, 0.68, 0.3]
elm1 = 0.3;
elm2 = 0.68;
elm3 = 1.0;
elm4 = elm3;
elm5 = elm2; 
elm6 = elm1;

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
    case 'chebyshev'
        w = window(@chebwin, N, 100);
    case 'custom'
        w = [elm1, elm2, elm3, elm4, elm5, elm6];
    otherwise
        w = ones(1,N);
end;

figure('name','Shading function')
stem(1:N, w);
title(['Shading function: ', shading_function], 'fontsize', szAxLabel)
xlabel('Elements, n', 'fontsize', szAxLabel)
ylabel('Sensitivity', 'fontsize', szAxLabel)
ylim([0 1.1]); xlim([1 N]);
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square';
print('-dpng', '-r300', './Shading function')


%% SIMULATION PARAMETERS AND VECTORS

rFarField = D.^2/(lambda_c); % Far field distance
rFarField_min = D.^2/(lambda_max);
rFarField_max = D.^2/(lambda_min);

rSensing = 4; % Simulation range
rFocus = rSensing; % Focusing range for near field distance
dtheta = 0.05; % Angular acquisition interval in degrees
theta_start = -45;
theta_end = 45;
theta_scan = theta_end - theta_start;
m = ceil(theta_scan/dtheta); % Number of angular bins
vecTheta = theta_start + dtheta*(0:m-1); % Angular vector

Ts = rSensing/c - Tp; % Sampling start time
Te = rSensing/c + Tp*1.3; % Sampling end time
n = 2*ceil((Te - Ts)/(2*dt)); % Number of time samples
t = Ts + dt*(0:n-1); % Time vector

df = 1/(n*dt); % Frequency interval
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
pt = pt.*exp(-1j*2*pi*fc*t); % Basebanding the transmitted signal

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
fprintf('Simulating the received signal at a distance of %0.2f m from array...\n', rSensing);
for i = 1:m
    i
    for k = -N/2:N/2-1
        R = (rSensing.^2+(k*d).^2 - 2*rSensing*(k*d)*cosd(90 - vecTheta(i))).^0.5;
        
        % The total time delay consists of the wave propagation delay and 
        % the beam-steering delay for the given element
        td = t - tau(k+N/2+1) - R/c;
        
        % Extra time delay for near-field compensation
        if nFieldComp == 1
            td = td + (((k*d).^2/(2*rFocus))/c);
        end;
        
        pha = 2*pi*(beta*td + K*td.^2);
        
        sti = (amp_error(k+N/2+1) * exp(1j*pha_error(k+N/2+1)) * w(k+N/2+1) * exp(1j*pha) .*(td >= 0 & td <= Tp) );
        
        st(i,:) = st(i,:) + sti;
    end;
    
    st(i,:) = st(i,:).*exp(-1j*2*pi*fc*t);
end;
fprintf('Simulation complete...\n');

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
colormap(hot(256))
imagesc(t,vecTheta,G);
title('Acquired signal at various angles', 'fontsize', szAxLabel2)
axis 'square'; axis 'xy'; colorbar;
xlabel('Time, t [s]', 'fontsize', szAxLabel)
ylabel('Angle, \theta [deg]', 'fontsize', szAxLabel)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);

%% RECEIVED SIGNAL PROCESSING
% A matched filter is used to compress the received signal at the sensing
% range.

f = df*(-n/2:n/2-1);

s0b = pt;
fs0b = fty(s0b).*exp(-1j*2*pi*f*(rSensing/c - Ts));

fsb = fty(st);
fsmb = fsb.*(ones(m,1)*conj(fs0b)); % Matched filtering and windowing
fsmb = fsmb.*(ones(m,1)*(abs(f) <= B0));
smb = ifty(fsmb);

tm=(rSensing/c)+dt*(-n/2:n/2-1); % Fast-time array after matched filtering

G=abs(smb);
xg=max(max(G)); ng=min(min(G)); cg=255/(xg-ng);

figure('name', 'SAS signal after range compression')
colormap(hot(256))
imagesc(c*tm, vecTheta, G);
title('Pressure intensity variation with angle at the sensing range', 'fontsize', szAxLabel2)
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

%% BEAM PATTERN EXTRACTION
% The values of the processed signal over the scan angle at the sensing
% range are extracted.

act_Beam = abs(smb(:,n/2+1)).';
act_Beam_normalised = act_Beam./max(act_Beam);
act_Beam_dB = 20*log10(act_Beam_normalised);
threedB_level = -3;

if beamInterp == 1
    interpBeam = real(interp(act_Beam,interpFactor));
    interpBeamNorm = interpBeam./max(interpBeam);
    interpBeamdB = 20*log10(interpBeamNorm);
    newVecTheta = theta_start + (dtheta)/interpFactor*(0:length(interpBeam)-1); % Angular vector

    figure('name','Measured radiation pattern - dB scale')
    plot(newVecTheta, interpBeamdB)
    title(['Radiation pattern at R = ',num2str(rSensing),'m (Far-field distance = ',num2str(rFarField,4),'m)',char(10),' with amp-error = ', num2str(20*log10(ampErrorInt(1)),2), 'dB; phase-error = +-', num2str(180/pi*phaErrorInt(2)), '^o; window = ', shading_function ], 'fontsize', szAxLabel2)
    xlabel('Angle, \theta [deg]', 'fontsize', szAxLabel)
    ylabel('Normalized level, [dB]', 'fontsize', szAxLabel)
    axis 'square'; axis 'tight';
    h_fig=get(gcf,'CurrentAxes');
    set(h_fig, 'fontsize', szAxScale);
    print('-dpdf', figRes, [outpath filename])
    print('-dpng', figRes, [outpath filename])

    figure('name','Measured radiation pattern - Polar plot')
    polar(newVecTheta*pi/180, interpBeamNorm)
    title(['Radiation pattern at ', num2str(rSensing,3), ' m'], 'fontsize', szAxLabel)
    axis 'square'; %axis tight
    h_fig=get(gcf,'CurrentAxes');
    set(h_fig, 'fontsize', szAxScale);
    
    threedB_points = [ newVecTheta(find(real(interpBeamdB) >= threedB_level, 1, 'first')),...
        newVecTheta(find(real(interpBeamdB) >= threedB_level, 1, 'last')) ];
    beamwidth = threedB_points(2) - threedB_points(1);

    fprintf('Theoretical beamwidth = %.2f deg\n', phi_d);
    fprintf('Measured beamwidth = %.2f deg\n', beamwidth);
else
    figure('name','Measured radiation pattern without interpolation')
    plot(vecTheta, act_Beam_dB)
    title(['Radiation pattern at R = ',num2str(rSensing),'m (Far-field distance = ',num2str(rFarField,4),'m) with amp-error = ', num2str(20*log10(ampErrorInt(1)),2), 'dB;', char(10),'phase-error = +-', num2str(180/pi*phaErrorInt(2)), '^o; window = ', shading_function ], 'fontsize', szAxLabel2)
    xlabel('Angle, \theta [deg]', 'fontsize', szAxLabel)
    ylabel('Normalized level, [dB]', 'fontsize', szAxLabel)
    axis 'square'; axis 'tight'
    h_fig=get(gcf,'CurrentAxes');
    set(h_fig, 'fontsize', szAxScale);
    print('-dpdf', figRes, [outpath filename])
    print('-dpng', figRes, [outpath filename])

    figure('name','Measured radiation pattern - Polar plot')
    polar(vecTheta*pi/180, act_Beam_normalised)
    title(['Radiation pattern at ', num2str(rSensing,3), ' m'], 'fontsize', szAxLabel)
    axis 'square'; %axis tight
    h_fig=get(gcf,'CurrentAxes');
    set(h_fig, 'fontsize', szAxScale);
    
    threedB_points = [ vecTheta(find(act_Beam_dB >= threedB_level, 1, 'first')),...
        vecTheta(find(act_Beam_dB >= threedB_level, 1, 'last')) ];
    beamwidth = threedB_points(2) - threedB_points(1);

    fprintf('Theoretical beamwidth = %.2f deg\n', phi_d);
    fprintf('Measured beamwidth = %.2f deg\n', beamwidth);    
end;

if nFieldComp == 1
    focusingDepth = [1/(1/rFocus+1/rFarField) 1/(1/rFocus - 1/rFarField)];
    fprintf('Depth of focus: %.2f m to %.2f m\n', focusingDepth(1), focusingDepth(2));
end;