%% MANIFOLD SIMULATION
% This code simulates the manifold obtained when the transmitter and
% receiver are co-located and a single ideal target is placed at some
% distance away from the transceiver.

%% OPTIONS
clear all
close all
clc

% Figure print parameters
szAxScale = 20;
szAxLabel = 20;

%% PULSE CONFIGURATION

c = 1500; % Propagation velocity in water
fc = 420e3; % Centre frequency
Tp = 50e-6; % Pulse width

lambda_c = c/fc;

%% TEMPORAL SAMPLING CONFIGURATION
fs = 20*fc; % Fast-time sampling rate
dt = 1/fs; % Fast-time sampling interval

%% RECEIVE-ARRAY CONFIGURATION

N = 64; % Number of elements in receive-array
d = 12.5e-3; % Element width
a = d; % Inter-element spacing
D = d*(N-1)+a; % Receive-array length (Synthetic aperture half-length)

du = D/N; % Element spacing (Synthetic aperture sampling interval)

%% RECEIVER BEAM PROPERTIES

rFarField = D.^2/(lambda_c); % Far field distance of receive array

%% RECEIVER POSITION VECTOR

xrx = 0.0; % Receiver array x-coordinate
yrx = (du*(-N/2:N/2-1)).'; % Receiver array y-coordinate position vector
zrx = 0; % Receiver array z-coordinate

rxCoords = [ones(N,1)*xrx yrx ones(N,1)*zrx];

%% TRANSMIT ELEMENT PROPERTIES

xtx = 0.0; % x-coordinate of the transmit element
ytx = 0.0; % y-coordinate of the transmit element
ztx = 0.1;  % z-coordinate of the transmit element

txCoord = [0, 0, ztx]; % Transmitter coordinate.

%% TARGET COORDINATES
% Target coordinates are specified in rectangular coordinates.

Rtarget = 4; % Line of sight range to target

ntargets = 1;
targets = [Rtarget, 0, 0];

%% IMAGING RANGE INTERVAL

Rmin = 0.95*Rtarget; % Minimum one-way range
Rmax = 1.05*Rtarget; % Maximum one-way range

Rc = Rtarget; % Mean range
R0 = (Rmax-Rmin)/2;

%% TEMPORAL SAMPLING VECTOR

Ts = 2*Rmin/c; % Fast-time sampling start time
Tf = 2*Rmax/c + Tp*1.1; % Fast-time sampling end time
n = 2*ceil((Tf - Ts)/(2*dt)); % Number of fast-time samples

df = 1/(n*dt);

t = Ts+dt*(0:n-1);
f = df*(-n/2:n/2-1);

%% THE TRANSMITTED PULSE

td = t-Ts;
pha = 2*pi*fc*td;
pt = exp(1j*pha).*(td>=0 & td <= Tp);

figure('name','Transmitted pulse')
plot(td, real(pt));
title('Transmitted pulse', 'fontsize', szAxLabel)
xlabel('Time, t [s]', 'fontsize', szAxLabel)
ylabel('Amplitude (Real), p(t)', 'fontsize', szAxLabel)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'; axis tight;

%% THE RECEIVED SIGNAL

% 2D array to hold the signal received by each receiving element:
stu = zeros(N, n); 

for i = 1:N % For each receiver element...
    for j = 1:ntargets % ... for each target ...
        
        % Vector connecting target to transmitter:
        temp1.x = targets(j,1) - txCoord(1);
        temp1.y = targets(j,2) - txCoord(2);
        temp1.z = targets(j,3) - txCoord(3);        
        temp1.abs = (temp1.x^2 + temp1.y^2 + temp1.z^2).^0.5;
        
        % Vector connecting target to ith receiver:
        temp2.x = targets(j,1) - rxCoords(i,1);
        temp2.y = targets(j,2) - rxCoords(i,2);
        temp2.z = targets(j,3) - rxCoords(i,3);
        temp2.abs = (temp2.x^2 + temp2.y^2 + temp2.z^2).^0.5;
        
        % Two-way delay from transmitter to target to ith receiver:
        R = (temp1.abs + temp2.abs);
        tau = R/c;
        td = t - tau;
        
        % The phase contribution of nth target:
        pha = 2*pi*fc*td;
        
        % Received signal:
        stu(i,:) = stu(i,:) + (exp(1j*pha).*(td <= Tp & td >= 0));
    end;
end;

G = real(stu);
figure('name', 'Manifold')
colormap(hot(256));
imagesc(c*t/2, 1:N, G);
title('Manifold', 'fontsize', szAxLabel)
xlabel('Range, [m]', 'fontsize', szAxLabel)
ylabel('Receiver elements, m', 'fontsize', szAxLabel)
axis 'square'; axis 'xy'; colorbar
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);

% figure('name', 'Received signal from reference element')
% plot(t,real(stu(Nelm/2+1,:)))
% title('Received signal from reference element', 'fontsize', szAxLabel2)
% xlabel('Time, t [s]', 'fontsize', szAxLabel1)
% ylabel('Amplitude (real), s_i(t)', 'fontsize', szAxLabel1)
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);
% axis 'square'; axis tight;


%% WRITE MANIFOLD TO BINARY FILE

filePath = '.';
fileName = 'TwoWayManifold';
fileExt = 're';

arg = [filePath '/' fileName '.' fileExt];

fid = fopen(arg, 'w');
fwrite(fid, G, 'double');
fclose(fid);
