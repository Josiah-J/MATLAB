%% FAN BEAM IMAGING: TEMPORAL-SPATIAL
% This code simulates fan beam imaging on a given plane using a single
% transmit element and an array of receive elements.

%% OPTIONS
clear all
close all
clc

% Figure print parameters
szAxScale = 20;
szAxLabel1 = 20;
szAxLabel2 = 14;

% Automatic near field compensation
autoNFComp = 1; 

%% PULSE CONFIGURATION

c = 1500; % Propagation velocity in water
fc = 350e3; % Centre frequency
B = 10e3; % Bandwidth
B0 = B/2; % Baseband bandwidth
Tp = 3e-3; % Chirp pulse width
K = B0/Tp; % Chirp rate

beta = fc - B0;
lambda_c = c/fc;
lambda_max = c/(fc-B0);
lambda_min = c/(fc+B0);

%% TEMPORAL SAMPLING CONFIGURATION
fs = 5*fc; % Fast-time sampling rate
dt = 1/fs; % Fast-time sampling interval

%% RECEIVE-ARRAY CONFIGURATION

N = 64; % Number of elements in receive-array
d = 12.5e-3; % Element width
a = d;%3.8e-3; % Inter-element spacing
D = d*(N-1)+a; % Receive-array length (Synthetic aperture half-length)

du = D/N; % Element spacing (Synthetic aperture sampling interval)
zrx = 0; % Receiver array z-coordinate

rxCoords = [zeros(N,1) (du*(-N/2:N/2-1)).' ones(N,1)*zrx];
u = du*(-N/2:N/2-1);

elShading = 'rectangular';
% Divergence angle of a single receive element
switch elShading
    case 'rectangular'
        phi_D = asind(0.88*lambda_c/d);
    case 'triangular'
        phi_D = asind(1.27*lambda_c/d);
    otherwise
        phi_D = asind(0.88*lambda_c/d);
end;

arrayShading = 'hamming';
switch arrayShading
    case 'rect'
        w = ones(N,1);
    case 'hamming'
        w = hamming(N);
    case 'hann'
        w = hann(N);
    case 'blackman-harris'
        w = window(@blackmanharris, N);
    case 'triangle'
        w = window(@triang, N);
    case 'chebyshev'
        w = window(@chebwin, N, 100);
    otherwise
        w = ones(N, 1);
end;

%% BEAM PROPERTIES

rFarField = D.^2/(lambda_c); % Far field distance of receive array
rFarField_min = D.^2/(lambda_max);
rFarField_max = D.^2/(lambda_min);

% Beam resolution
theta_3dB = asind(lambda_c/D); % Without aperture weighting

% Maximum view angle
theta_max = 2*asind(lambda_c/(2*du));

Nbeams = N*2-1; % Number of beams to be formed

%% IMAGING RANGE INTERVAL

Rmin = 10; % Minimum one-way range
Rmax = 20; % Maximum one-way range

Rc = Rmin + (Rmax - Rmin)/2; % Mean range
R0 = (Rmax-Rmin)/2;

if Rc < rFarField && autoNFComp == 1
    nFieldComp = 1;
else
    nFieldComp = 0;
end;

%% TEMPORAL SAMPLING VECTOR

Ts = 2*Rmin/c; % Fast-time sampling start time
Tf = 2*Rmax/c + Tp*1.1; % Fast-time sampling end time
n = 2*ceil((Tf - Ts)/(2*dt)); % Number of fast-time samples

df = 1/(n*dt);

t = Ts+dt*(0:n-1);
f = df*(-n/2:n/2-1);

%% TRANSMIT ELEMENT PROPERTIES

ztx = 0.1;  % z-coordinate of the transmit element
txCoord = [0, 0, ztx]; % Transmitter coordinate.

%% TARGET REGION PROPERTIES

Bmax = Rmax*sind(phi_D); % Maximum beamwidth of single element at the maximum range

Xc = Rc; %10; % Centre of the target region in the range direction
X0 = R0; %1; % Half swath-width of the target region

Yc = 0; % Centre of the target region in the cross-range direction
Y0 = D + Bmax/2; % Half-width of the target region in the cross-range direction

%% TARGET COORDINATES
% Target coordinates are specified in rectangular coordinates.

% ntargets = 1;
% targets = [Rc, 0, 0];

% ntargets = 1;
% targets = [Rc, 0.2*Y0, 0];

% ntargets = 2;
% targets = [Rc, 0.2*Y0, 0; Rc, -0.2*Y0, 0];

% ntargets = 3;
% targets = [Xc, 0, 0; Xc, 0.2*Y0, 0; Xc, -0.2*Y0, 0];

% ntargets = 3;
% targets = [Xc, 0, 0; Xc-15e-3, 0, 0; Xc+15e-3, 0, 0];

ntargets = 9;
targets = [Xc, 0, 0; Xc, 0.3*Y0, 0; Xc, -0.3*Y0, 0; ...
    Xc-0.5*X0, 0, 0; Xc-0.5*X0, 0.3*Y0, 0; Xc-0.5*X0, -0.3*Y0, 0; ...
    Xc+0.5*X0, 0, 0; Xc+0.5*X0, 0.3*Y0, 0; Xc+0.5*X0, -0.3*Y0, 0];

% ntargets = 8;
% targets = [10, 0, 0; 10, 0.5, 0;...
%     20, 0, 0; 20, 0.5, 0;...
%     30, 0, 0; 30, 0.5, 0;...
%     50, 0, 0; 50, 0.5, 0];

% ntargets = 5;
% targets = [Rc, 0, 0; Rc, 0.5*Y0, 0; Rc, -0.5*Y0, 0; Rc-0.5*R0, 0, 0; Rc+0.5*R0, 0, 0];

%% THE TRANSMITTED PULSE

td = t-Ts;
pha = 2*pi*(beta*(td) + K*td.^2);
pt = exp(1j*pha).*(td>=0 & td <= Tp);

pt = pt.*(exp(-1j*2*pi*fc*td)); % Basebanding
Pf = fftshift(fft(pt));

figure('name','Transmitted pulse (basebanded)')
plot(td, real(pt));
title('Transmitted pulse', 'fontsize', szAxLabel2)
xlabel('Time, t [s]', 'fontsize', szAxLabel1)
ylabel('Amplitude (Real), p(t)', 'fontsize', szAxLabel1)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'; axis tight;

figure('name', 'Magnitude spectrum of basebanded transmitted signal')
plot(f, abs(Pf));
title('Magnitude spectrum of transmitted pulse', 'fontsize', szAxLabel2)
xlabel('Frequency, f [Hz]', 'fontsize', szAxLabel1)
ylabel('Magnitude, |P(f)|', 'fontsize', szAxLabel1)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'; axis tight

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
        pha = 2*pi*(beta*td + K*td.^2);
        
        % Received signal:
        stu(i,:) = stu(i,:) + (exp(1j*pha).*(td <= Tp & td >= 0));
    end;
end;

G = real(stu);
figure('name', 'The received signal')
colormap(hot(256));
imagesc(t, u, G);
title('The received signal', 'fontsize', szAxLabel2)
xlabel('Time, t [sec]', 'fontsize', szAxLabel1)
ylabel('Receive aperture, u [m]', 'fontsize', szAxLabel1)
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

%% PREPROCESSING - RECEIVED SIGNAL BASEBANDING

% Baseband conversion
stu = stu .* (ones(N,1)*exp(-1j*2*pi*fc*t));

G = real(stu);
figure('name', 'The received signal after basebanding')
colormap(hot(256));
imagesc(t, u, G);
title('The received signal after basebanding', 'fontsize', szAxLabel2)
xlabel('Time, t [sec]', 'fontsize', szAxLabel1)
ylabel('Receive aperture, u [m]', 'fontsize', szAxLabel1)
axis 'square'; axis 'xy'; colorbar
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);

%% PREPROCESSING - ARRAY SHADING FUNCTION

% Applying the array shading function to the received data
w = w*ones(1,n);
stu = stu.*w;

G = real(stu);
figure('name', 'The received signal after applying artificial shading')
colormap(hot(256));
imagesc(t, u, G);
title(['The received signal after applying artificial shading - ', arrayShading], 'fontsize', szAxLabel2)
xlabel('Time, t [sec]', 'fontsize', szAxLabel1)
ylabel('Receive aperture, u [m]', 'fontsize', szAxLabel1)
axis 'square'; axis 'xy'; colorbar
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);

%% RANGE COMPRESSION

Sfu = fftshift((fft(stu.')).',2);

Sfu = Sfu.*(ones(N,1)*(conj(Pf)));

stu = (ifft(fftshift(Sfu,2).')).';

G = abs(stu);
figure('name','Range-compressed data')
colormap(hot(256));
% imagesc(G);
imagesc(c*t, u, G);
title('The received signal after range compression', 'fontsize', szAxLabel2)
xlabel('Two-way range, R = ct [m]', 'fontsize', szAxLabel1)
ylabel('Receive aperture, u [m]', 'fontsize', szAxLabel1)
axis 'square'; axis 'xy'; colorbar
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);

%% NEAR-FIELD COMPENSATION

Sfu = fftshift((fft(stu.')).',2);

% Compensation for near-field effects
if nFieldComp == 1
    compFilter = zeros(N,n);
    for k = -N/2:N/2-1
        td_comp = ((k*du).^2/(2*Xc))/c;
        compFilter(k+N/2+1,:) = exp(+1j*2*pi*(f+fc)*td_comp);
    end;
    
    Sfu = Sfu.*compFilter;
end;

stu = (ifft(fftshift(Sfu,2).')).';

%% DELAY-SUM BEAMFORMER

% deltaR = 0.25*(c/(2*B));
% Nr = 2*ceil((Rmax - Rmin)/(2*deltaR));
% r = Rmin + deltaR*(0:Nr-1);
% 
% fBeam = zeros(Ntheta, Nr); % The reconstruction matrix with N-beams and n range bins
% 
% % Beam angle vector
% theta_s = dtheta*(-Ntheta/2:Ntheta/2-1);
% 
% for i = -Ntheta/2:Ntheta/2-1 % For each beam angle ...
%     i+Ntheta/2+1
% 
%     for j = 1:Nr % ... at each range ...
%         Sfoc = 0.0;
%         
%         % The reference coordinate in rectangular coordinates:
%         % (x =  Rcos(theta), y = Rsin(theta), z)
%         refCoord = [r(j)*cosd(theta_s(i+Ntheta/2+1)), r(j)*sind(theta_s(i+Ntheta/2+1)), zrx];        
%         
%         % Vector connecting transmitter to the reference coordinate:
%         temp1.x = refCoord(1) - txCoord(1);
%         temp1.y = refCoord(2) - txCoord(2);
%         temp1.z = refCoord(3) - txCoord(3);
%         temp1.abs = (temp1.x^2+temp1.y^2+temp1.z^2).^0.5;        
%         
%         % Vector connecting reference (middle) receiver to the reference
%         % coordinate:
%         temp2.x = refCoord(1) - rxCoords(Nelm/2+1,1);
%         temp2.y = refCoord(2) - rxCoords(Nelm/2+1,2);
%         temp2.z = refCoord(3) - rxCoords(Nelm/2+1,3);
%         temp2.abs = (temp2.x^2+temp2.y^2+temp2.z^2).^0.5;
%         
%         % Two-way range to the reference coordinate
%         rRef = temp1.abs + temp2.abs;
%         
%         % Two-way delay to the reference coordinate
%         tRef = rRef/c;
%         
%         for k = -Nelm/2:Nelm/2-1%0% % ... for each receiving element...
%             % Vector connecting the kth receiver to the reference
%             % coordinate:
%             temp2.x = refCoord(1) - rxCoords(k+Nelm/2+1,1);
%             temp2.y = refCoord(2) - rxCoords(k+Nelm/2+1,2);
%             temp2.z = refCoord(3) - rxCoords(k+Nelm/2+1,3);
%             temp2.abs = (temp2.x^2+temp2.y^2+temp2.z^2).^0.5;
%             
%             % Two-way delay from transmitter to reference coordinate and
%             % back to the kth receiver:
%             td = (temp1.abs+temp2.abs)/c;
%             
%             % Time index of two-way delay:
%             td_idx = 2*ceil((td-Ts)/(2*dt));
%             
%             if (td_idx > 0 && td_idx <= n)
%                 Sfoc = Sfoc + ( stu(k+Nelm/2+1,td_idx) * exp(1j*2*pi*fc*(td - tRef)) );
%             end;
%         end;
%         
%         fBeam(i+Ntheta/2+1,j) = Sfoc;
%     end;
% end;
% 
% G = abs(fBeam);
% figure('name','Fan-beam image in polar coordinates')
% colormap(hot(256));
% imagesc(r, theta_s, G);
% title('The fan-beam image in polar coordinates', 'fontsize', szAxLabel2)
% xlabel('Range, R [m]', 'fontsize', szAxLabel1)
% ylabel('Angle, \theta [deg]', 'fontsize', szAxLabel1)
% axis 'square'; axis 'xy'; colorbar
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);

%% RESAMPLING TO RECTANGULAR COORDINATES

% Xmin = Xc-X0;
% Xmax = Xc+X0;
% dX = deltaR*0.5;
% nX = 2*ceil((Xmax - Xmin)/(2*dX));
% X = Xmin+dX*(0:nX-1);
% 
% dY = du*0.5;
% nY = 2*ceil(Y0/dY);
% Y = dY*(-nY/2:nY/2-1);
% 
% fBeamRspd = max(abs(fBeam(:)))*ones(nY, nX);
% 
% for i = -nY/2:nY/2-1%0%
%     for j = -nX/2:nX/2-1%0%
%         
%         R = ( X(j+nX/2+1)^2 + Y(i+nY/2+1)^2  ).^0.5;
%         R_idx = ceil((R - Rmin)/deltaR);
%         if R_idx <= 0 || R_idx > Nr
%             continue
%         end;
%         
%         my_theta = atand(Y(i+nY/2+1)/X(j+nX/2+1));
%         my_theta_idx = ceil(my_theta/dtheta)+Ntheta/2+1;
%         if my_theta_idx <= 0 || my_theta_idx > Ntheta
%             continue
%         end; 
%         
%         fBeamRspd(i+nY/2+1, j+nX/2+1) = fBeam(my_theta_idx, R_idx);
%     end;
% end;
% 
% G = abs(fBeamRspd).';
% figure('name','Fan-beam image in rectangular coordinates')
% colormap(hot(256));
% imagesc(Y, X, G);
% title('The fan-beam image in rectangular coordinates', 'fontsize', szAxLabel2)
% xlabel('Azimuth, y [m]', 'fontsize', szAxLabel1)
% ylabel('Range, x [m]', 'fontsize', szAxLabel1)
% axis 'square'; colorbar
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);

%% FREQUENCY DOMAIN BEAMFORMING

% Zero pad in the spatial domain
m = Nbeams;
newD = d*(m - 1) + a;
newdu = newD/m;

u = newdu*(-m/2:m/2-1);

stu = [ zeros(floor((m - N)/2),n); stu; zeros(ceil((m - N)/2),n) ];

G = abs(stu);
figure('name','Range-compressed data after zero-padding')
colormap(hot(256));
imagesc(c*t/2, u, G);
title('Range-compressed data after zero-padding', 'fontsize', szAxLabel2)
xlabel('Two-way range, R = ct [m]', 'fontsize', szAxLabel1)
ylabel('Receive aperture, u [m]', 'fontsize', szAxLabel1)
axis 'square'; axis 'xy'; colorbar
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);

Sfu = fftshift((fft(stu.')).');

Sfku = fftshift(fft(Sfu),1);

Stku = ifft(fftshift(Sfku.',1)).';

dtheta = asind(lambda_c/newD);%*N/(m));
newVecTheta = dtheta*(-m/2:m/2-1);

G = abs(Stku);
figure('name','Planar image')
colormap(hot(256));
% imagesc(G);
imagesc(c*t/2, newVecTheta, G);
title('Planar image', 'fontsize', szAxLabel2)
xlabel('Range, ct [m]', 'fontsize', szAxLabel1)
ylabel('Angle, \theta [deg]', 'fontsize', szAxLabel1)
axis 'square'; axis 'xy'; colorbar
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);

rSlice = Rc;
rBin = 2*ceil((2*rSlice/c - Ts)/(2*dt));

% figure('name','Beam profile at a given range')
% plot(newVecTheta, abs(Stku(:,rBin)));
% title('Beam profile at a given range', 'fontsize', szAxLabel2)
% xlabel('Angle, \theta [deg]', 'fontsize', szAxLabel1)
% ylabel('Magnitude', 'fontsize', szAxLabel1)
% axis 'square'; axis 'tight'
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);

figure('name','Beam profile at a given range - dB scale')
plot(newVecTheta, 20*log10(abs(Stku(:,rBin))./max(abs(Stku(:,rBin)))));
title('Beam profile at a given range', 'fontsize', szAxLabel2)
xlabel('Angle, \theta [deg]', 'fontsize', szAxLabel1)
ylabel('Magnitude [dB]', 'fontsize', szAxLabel1)
axis 'square'; axis 'tight'
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);

if nFieldComp == 1
    focusingDepth = [1/(1/Xc+1/rFarField) 1/(1/Xc - 1/rFarField)];
    fprintf('Depth of focus: %.2f m to %.2f m\n', focusingDepth(1), focusingDepth(2));
end;

%% WRITE DATA TO BINARY FILE

filePath = '.';
fileName = 'Fanbeam';
fileExt = 'bin';
fileData = G.';

arg = [filePath '/' fileName '.' fileExt];

fid = fopen(arg, 'w');
fwrite(fid, fileData, 'double');
fclose(fid);
