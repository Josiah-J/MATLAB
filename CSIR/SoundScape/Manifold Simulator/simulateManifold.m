%% MANIFOLD SIMULATION
% This code simulates the manifold obtained from a transmitter placed at a
% given range away from the receiver array.

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
zrx = 0.0; % Receiver array z-coordinate

rxCoords = [ones(N,1)*xrx yrx ones(N,1)*zrx];

%% TRANSMIT ELEMENT PROPERTIES

rtx = 4; % Line of sight range to transmitter
ytx = 0.0; % y-coordinate of the transmit element
ztx = 0.0;  % z-coordinate of the transmit element

txCoord = [rtx, ytx, ztx]; % Transmitter coordinate.

%% IMAGING RANGE INTERVAL

Rmin = 0.95*rtx; % Minimum one-way range
Rmax = rtx*1.05; % Maximum one-way range

Rc = Rmin + (Rmax - Rmin)/2; % Mean range
R0 = (Rmax - Rmin)/2; 

%% TEMPORAL SAMPLING VECTOR

Ts = Rmin/c; % Fast-time sampling start time
Tf = Rmax/c + Tp*1.1; % Fast-time sampling end time
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
         
        % Vector connecting target to transmitter:
        temp1.x = txCoord(1) - rxCoords(i,1);
        temp1.y = txCoord(2) - rxCoords(i,2);
        temp1.z = txCoord(3) - rxCoords(i,3);        
        temp1.abs = (temp1.x^2 + temp1.y^2 + temp1.z^2).^0.5;
                
        % One-way delay from transmitter to ith receiver:
        R = temp1.abs;
        tau = R/c;
        td = t - tau;
        
        % The phase contribution of nth target:
        pha = 2*pi*fc*td;
        
        % Received signal:
        stu(i,:) = stu(i,:) + (exp(1j*pha).*(td <= Tp & td >= 0));
end;

G = real(stu);
figure('name', 'Manifold')
colormap(jet(256));
imagesc(c*t, 1:N, G);
title('Manifold', 'fontsize', szAxLabel)
xlabel('Range, R [m]', 'fontsize', szAxLabel)
ylabel('Receiver elements, m', 'fontsize', szAxLabel)
axis 'square'; axis 'xy'; colorbar
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
% 
% figure('name', 'Received signal from reference element')
% plot(t,real(stu(N/2+1,:)))
% title('Received signal from reference element', 'fontsize', szAxLabel2)
% xlabel('Time, t [s]', 'fontsize', szAxLabel1)
% ylabel('Amplitude (real), s_i(t)', 'fontsize', szAxLabel1)
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);
% axis 'square'; axis tight;

%% WRITE DATA TO BINARY FILE

filePath = '.';
fileName = 'OneWayManifold';
fileExt = 're';

arg = [filePath '/' fileName '.' fileExt];

fid = fopen(arg, 'w');
fwrite(fid, G.', 'double');
fclose(fid);

