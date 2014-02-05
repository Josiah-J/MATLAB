%% 3D VOLUMETRIC SONAR IMAGING USING THE SINGLE PING METHODOLOGY
% Each elevation beam is allocated a given frequency.  The transmitted
% waveform is generated such that all the beams are encoded in one pulse.
% The 3D sonar volume is built by processing the returns from each beam
% using conventional 2D beamforming techniques.
% Note:
% No basebanding of the received data is performed.
% Spatial beamforming precedes pulse compression.

%% GENERAL OPTIONS
close all
clear all
clc

% Figure print parameters
szAxScale = 20;
szAxLabel1 = 20;
szAxLabel2 = 14;

% Automatic near field compensation
autoNFComp = 1; 

% Debug mode
debug = 1;

%% GENERAL SONAR OPTIONS

c = 1500; % - speed of propagation
f0 = 130e3; % - start frequency
f1 = 330e3; % - stop frequency
Tp = 1.5e-3; % - Pulse duration

% Derived options:
fc = 0.5*(f0+f1); % Centre frequency
B = f1 - f0; % - Total bandwidth
B0 = B/2; % - Baseband bandwidth
lambda_c = c/fc; % Wavelength at centre frequency
lambda_max = c/f0; % Wavelength at start frequency
lambda_min = c/f1; % Wavelength at stop frequency

%% TRANSMIT ARRAY PROPERTIES

tx.M = 32; % - number of elements
tx.d = 25e-3; % - inter-element spacing
tx.a = tx.d; % - element width
tx.N = 31; % - number of beams
tx.Pos.xOffset = 0;
tx.Pos.yOffset = 0;
tx.Pos.zOffset = 0;
tx.Ref = tx.M/2 + 1;

% Derived properties:
tx.D = tx.d*tx.M;
if( tx.N > (2*tx.M - 1) || ~mod(tx.N,2) ) tx.N = 2*tx.M - 1; end; % - number of beams (ensure it is within the theoretical limit)
tx.Pos.x = tx.Pos.xOffset + zeros(1, tx.M); % - position vector of each element
tx.Pos.y = tx.Pos.yOffset + zeros(1, tx.M);
tx.Pos.z = tx.Pos.zOffset + ((0:tx.M-1) - tx.M/2).*tx.d;

tx.Beam.BFraction = 0.7;
tx.Beam.B = tx.Beam.BFraction*(f1 - f0)/tx.N;
tx.Beam.K = 0.5*tx.Beam.B/Tp;
tx.Beam.f0 = f0 + ( 0:(tx.N - 1)) * ((f1 - f0)/tx.N); % - beam bandwidth (start frequency and stop frequency) vector
tx.Beam.f1 = f0 + ( 0:(tx.N - 1)) * ((f1 - f0)/tx.N) + tx.Beam.B;
tx.Beam.fc = 0.5*(tx.Beam.f1 + tx.Beam.f0);
tx.Beam.lambda =  c./tx.Beam.f1;
tx.Beam.thetaMax = 2 * asin( tx.Beam.lambda/(2*tx.d) ); % - divergence angle
tx.Beam.theta3dB = asin( tx.Beam.lambda/tx.D );
tx.Beam.dtheta = min(tx.Beam.thetaMax)/tx.N;
tx.Beam.angle = ( -(tx.N-1)/2:(tx.N-1)/2 ) * tx.Beam.dtheta;
tx.Beam.Fresnel = (tx.D)^2/min(tx.Beam.lambda);

%% RECEIVE ARRAY PROPERTIES

rx.M = 64; % - number of elements
rx.d = 12.5e-3; % - inter-element spacing
rx.a = rx.d; % - element width
rx.N = 2*rx.M - 1; % - number of beams
rx.Shading.Type = 'hann'; % - shading function
rx.AmpError.Type = 'none'; % - array amplitude error vector
rx.PhaseError.Type = 'none'; % - array phase error vector
rx.Pos.xOffset = 0;
rx.Pos.yOffset = 0;
rx.Pos.zOffset = -(tx.D*1.05)/2;
rx.Ref = rx.M/2 + 1;

% Derived properties:
rx.D = rx.M*rx.d;
if(rx.N > (2*rx.M - 1) || ~mod(rx.N,2) ) rx.N = (2*rx.M - 1); end; % - number of beams (ensure it is within the theoretical limit)
rx.Pos.x = rx.Pos.xOffset + zeros(1, rx.M); % - position vector of each element
rx.Pos.y = rx.Pos.yOffset + ((0:rx.M-1) - rx.M/2).*rx.d;
rx.Pos.z = rx.Pos.zOffset + zeros(1, rx.M);

rx.Beam.thetaMax = 2*asin( lambda_min/(2*rx.d) ); % - array divergence angle
rx.Beam.theta3dB = asin(lambda_min/rx.D);
rx.Beam.dtheta = rx.Beam.thetaMax/rx.N;
rx.Beam.angle = ( -(rx.N - 1)/2:(rx.N - 1)/2 ) * rx.Beam.dtheta;
rx.Beam.Fresnel = (rx.D)^2/lambda_min;

%% IMAGING VOLUME PROPERTIES

R0 = 3.0; % - Start range
R1 = 6.0; % - stop range

Rc = 0.5*(R0 + R1); % Focal range
Y0 = 0.5*Rc*rx.Beam.thetaMax; % - Width (within the receive array max. steering angle at the given range)
Z0 = 0.5*Rc*min(tx.Beam.thetaMax); % - Height (within the transmit array max steering angle at the given range)

NearTarget.Pos.x = R0; NearTarget.Pos.y = 0; NearTarget.Pos.z = 0;
Tx2Target.x = NearTarget.Pos.x - tx.Pos.x(tx.Ref);
Tx2Target.y = NearTarget.Pos.y - tx.Pos.y(tx.Ref);
Tx2Target.z = NearTarget.Pos.z - tx.Pos.z(tx.Ref);
Tx2Target.range = sqrt(Tx2Target.x^2 + Tx2Target.y^2 + Tx2Target.z^2);

Rx2Target.x = NearTarget.Pos.x - rx.Pos.x(rx.Ref);
Rx2Target.y = NearTarget.Pos.y - rx.Pos.y(rx.Ref);
Rx2Target.z = NearTarget.Pos.z - rx.Pos.z(rx.Ref);
Rx2Target.range = sqrt(Rx2Target.x^2 + Rx2Target.y^2 + Rx2Target.z^2);

X0 = (Tx2Target.range + Rx2Target.range);

FarTarget.Pos.x = R1; FarTarget.Pos.y = Y0; FarTarget.Pos.z = Z0;
Tx2Target.x = FarTarget.Pos.x - tx.Pos.x(tx.Ref);
Tx2Target.y = FarTarget.Pos.y - tx.Pos.y(tx.Ref);
Tx2Target.z = FarTarget.Pos.z - tx.Pos.z(tx.Ref);
Tx2Target.range = sqrt(Tx2Target.x^2 + Tx2Target.y^2 + Tx2Target.z^2);

Rx2Target.x = FarTarget.Pos.x - rx.Pos.x(rx.Ref);
Rx2Target.y = FarTarget.Pos.y - rx.Pos.y(rx.Ref);
Rx2Target.z = FarTarget.Pos.z - rx.Pos.z(rx.Ref);
Rx2Target.range = sqrt(Rx2Target.x^2 + Rx2Target.y^2 + Rx2Target.z^2);

X1 = (Tx2Target.range + Rx2Target.range);

Xc = 0.5*(X0+X1);

%% TEMPORAL SAMPLING PROPERTIES

fs = 5*f1; % - Temporal sampling frequency
dt = 1/fs; % - Temporal Sampling interval

Ts = X0/c; % - sampling start time
Te = X1/c + 1.1*Tp; % - Sampling end time

n = 2*ceil((Te - Ts)/(2*dt)); % - number of temporal samples

df = 1/(n*dt); % - frequency rate

t = Ts + (0:n-1)*dt; % - time vector
f = df*(-n/2:n/2-1); % - frequency vector

%% TARGET COORDINATES

targets.N = 3;
targets.Pos.x = [ Rc, Rc, Rc];
targets.Pos.y = [ 0, 0.5*Y0, 0];
targets.Pos.z = [ 0, 0.35*Z0, -0.35*Z0];
targets.Pos.r = sqrt(targets.Pos.x.^2 + targets.Pos.y.^2 + targets.Pos.z.^2);
targets.Pos.phi = pi/2 - acos(targets.Pos.z./targets.Pos.r);
targets.Pos.theta = atan(targets.Pos.y./targets.Pos.r);

%% THE TRANSMIT PULSE

% K = 0.5*tx.Beam.B/Tp;

pt = zeros(tx.M, n);

for i = 1:tx.M
    cValue = zeros(1, n);
    for j = 1:tx.N
        beamDelay = (tx.M - i)*tx.d/c * sin(tx.Beam.angle(j));
        
        if(Rc < tx.Beam.Fresnel) 
            fresnelComp = ((i - tx.M/2 - 1)*tx.d)/(2*Rc)/c; 
        else
            fresnelComp = 0.0;
        end;
        
        td = t - Ts + fresnelComp - beamDelay - (n*dt - Tp)/2;
        
        beta = tx.Beam.f0(j);
        
        pha = 2*pi*(beta*td + tx.Beam.K*td.^2);
        
        cValue = cValue + exp(1j*pha).*(td >= 0 & td <= Tp);%.*(exp(-1j*2*pi*tx.Beam.fc(j)*(t-Ts)));        
    end;
    
    pt(i, :) = cValue;
    
    plot(t-Ts, real(cValue));
    title(['Transmitted pulse: Element ', num2str(i)], 'fontsize', szAxLabel2)
    xlabel('Time, t [s]', 'fontsize', szAxLabel1)
    ylabel('Amplitude (Real), p(t)', 'fontsize', szAxLabel1)
    h_fig=get(gcf,'CurrentAxes');
    set(h_fig, 'fontsize', szAxScale);
    axis 'square'; axis tight;
    pause(0.2);
end;

Pf = fftshift(fft(pt(tx.Ref,:)));
figure('name', 'Magnitude spectrum of transmitted signal')
plot(f/1000, abs(Pf));
title('Magnitude spectrum of transmitted pulse', 'fontsize', szAxLabel2)
xlabel('Frequency, f [kHz]', 'fontsize', szAxLabel1)
ylabel('Magnitude, |P(f)|', 'fontsize', szAxLabel1)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'; axis tight

%% THE RECEIVED SIGNAL

% 2D array to hold the signal received by the receiver array
stu = zeros(rx.M, n);

for i = 1:tx.N
    fprintf('Simulating received signal from transmit beam at %.2f deg\n',tx.Beam.angle(i)*180/pi);
    for j = 1:targets.N
        withinTxBeam = (abs(targets.Pos.phi(j) - tx.Beam.angle(i)) <= tx.Beam.theta3dB(i)/2);
        
        for k = 1:rx.M
            % Vector connecting targets to transmitter
            Tx2Target.x = targets.Pos.x(j) - tx.Pos.x(tx.Ref);
            Tx2Target.y = targets.Pos.y(j) - tx.Pos.y(tx.Ref);
            Tx2Target.z = targets.Pos.z(j) - tx.Pos.z(tx.Ref);
            Tx2Target.range = sqrt(Tx2Target.x^2 + Tx2Target.y^2 + Tx2Target.z^2);

            % Vector connecting target to kth receiver
            Rx2Target.x = targets.Pos.x(j) - rx.Pos.x(k);
            Rx2Target.y = targets.Pos.y(j) - rx.Pos.y(k);
            Rx2Target.z = targets.Pos.z(j) - rx.Pos.z(k);
            Rx2Target.range = sqrt(Rx2Target.x^2 + Rx2Target.y^2 + Rx2Target.z^2);

            % Two-way delay from transmitter to target to kth element
            R = Tx2Target.range + Rx2Target.range;
            tau = R/c;
            td = t - tau;

            % The phase contribution of the jth target
            beta = tx.Beam.f0(i);
            pha = 2*pi*(beta*td + tx.Beam.K*td.^2);

            stu(k,:) = stu(k,:) + exp(1j*pha).*(td >= 0 & td <= Tp) * withinTxBeam;
        end;
    end;
end;

G = real(stu);
figure('name', 'The received signal')
colormap(hot(256));
imagesc(t, 1:rx.M, G);
title('The received signal', 'fontsize', szAxLabel2)
xlabel('Time, t [sec]', 'fontsize', szAxLabel1)
ylabel('Receiver element index', 'fontsize', szAxLabel1)
axis 'square'; axis 'xy'; colorbar
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);

%% ARRAY SHADING FUNCTION

% Applying the array function to the received data
switch rx.Shading.Type
    case 'rect'
        w = ones(rx.M, 1);
    case 'hamming'
        w = hamming(rx.M);
    case 'hann'
        w = hann(rx.M);
    case 'blackman-harris'
        w = window(@blackmanharris, rx.M);
    case 'triangle'
        w = window(@triang, rx.M);
    case 'chebyshev'
        w = window(@chebwin, rx.M, 100);
    otherwise
        w = ones(rx.M, 1);
end;

w = w*ones(1, n);
stu = stu.*w;

%% NEAR-FIELD COMPENSATION

RxToFocalRange.x = Rc - rx.Pos.x(rx.Ref);
RxToFocalRange.y = -rx.Pos.y(rx.Ref);
RxToFocalRange.z = -rx.Pos.z(rx.Ref);
RxToFocalRange.r = sqrt(RxToFocalRange.x^2 + RxToFocalRange.y^2 + RxToFocalRange.z^2);

if RxToFocalRange.r < rx.Beam.Fresnel && autoNFComp == 1
    
    Sfu = fftshift((fft(stu.')).', 2);

    % Compensation for near-field effects
    compFilter = zeros(rx.M, n);
    for k = -rx.M/2:rx.M/2-1
        td_comp = ((k*rx.d).^2/(2*RxToFocalRange.r))/c;
        compFilter(k+rx.M/2+1, :) = exp(+1j*2*pi*(f)*td_comp);
    end;

    Sfu = Sfu.*compFilter;

    stu = (ifft(fftshift(Sfu,2).')).';

    if (debug == 1)
        G = real(stu);
        figure('name', ['Near field compensation'])
        colormap(hot(256));
        imagesc(c*t/2, 1:rx.M, G);
        title(['Near field compensation'], 'fontsize', szAxLabel2)
        xlabel('Time, t [sec]', 'fontsize', szAxLabel1)
        ylabel('Receiver element index', 'fontsize', szAxLabel1)
        axis 'square'; axis 'xy'; colorbar
        h_fig=get(gcf,'CurrentAxes');
        set(h_fig, 'fontsize', szAxScale);
    end;    
end;

%% FREQUENCY DOMAIN BEAMFORMING

% Zero pad in the spatial domain

stuphi = [ zeros(floor((rx.N - rx.M)/2),n); stu; zeros(ceil((rx.N - rx.M)/2),n) ];

Sfu = fftshift((fft(stuphi.')).');
Sfku = fftshift(fft(Sfu),1);
Stku = ifft(fftshift(Sfku.',1)).';

stuphi = Stku;

if (debug == 1)
    G = abs(stuphi);
    figure('name',['After beamforming'])
    colormap(hot(256));
    % imagesc(G);
    imagesc(c*t/2, 1:rx.N, G);
    title(['Beam-time domain'], 'fontsize', szAxLabel2)
    xlabel('Range, ct [m]', 'fontsize', szAxLabel1)
    ylabel('Angle, \theta [deg]', 'fontsize', szAxLabel1)
    axis 'square'; axis 'xy'; colorbar
    h_fig=get(gcf,'CurrentAxes');
    set(h_fig, 'fontsize', szAxScale);
end;

%% RANGE COMPRESSION

stkuphi = zeros(rx.N, n, tx.N);

td = t-Ts;
for i = 1:tx.N  
    beta = tx.Beam.f0(i);
    
    pha = 2*pi*(beta*td + tx.Beam.K*td.^2);
    pt = exp(1j*pha).*(td >= 0 & td <= Tp);

%     pt = pt.*(exp(-1j*2*pi*tx.Beam.fc(i)*td)); % Basebanding
    Pf = fftshift(fft(pt));
    
    Sfkuphi = fftshift((fft(stuphi.')).',2);
    temp = Sfkuphi.*(ones(rx.N, 1)*conj(Pf));
    
    stkuphi(:, :, i) = (ifft(fftshift(temp,2).')).';
    
    if (i < 4 && debug == 1)
        G = abs(stkuphi(:,:,i));
        figure('name',['Range compression image for beam ', num2str(tx.Beam.angle(i)*180/pi), ' deg'])
        colormap(hot(256));
        imagesc(c*t/2, 1:rx.N, G);
        title(['2D image of range compressed data for beam ',num2str(tx.Beam.angle(i)*180/pi), ' deg'], 'fontsize', szAxLabel2)
        xlabel('Range, ct [m]', 'fontsize', szAxLabel1)
        ylabel('Angle, \theta [deg]', 'fontsize', szAxLabel1)
        axis 'square'; axis 'xy'; colorbar
        h_fig=get(gcf,'CurrentAxes');
        set(h_fig, 'fontsize', szAxScale);
        pause(2)
    end;
end;

%% BASEBANDING

% for i = 1:tx.N
%     %fprintf('Beam band: %.2f kHz - %.2f kHz\n', tx.Beam.f0(i)/1000.0, tx.Beam.f1(i)/1000);
%     
%     bbKernel = ones(rx.N, 1)*exp(-1j*2*pi*tx.Beam.fc(i)*t);
%     temp = stkuphi(:,:,i).*bbKernel;
%     
%     ftemp = fftshift((fft(temp.')).',2).*(ones(rx.N,1)*(abs(f) <= tx.Beam.B/2));
%     
%     stkuphi(:, :, i) = (ifft(fftshift(ftemp,2).')).';
%     
%     if (i == (tx.N-1)/2+1 && debug == 1)     
% %         figure('name', 'The received signal')
% %         plot(t, real(stkuphi(rx.Ref,:,i)));
% %         title(['Received signal at ', num2str(tx.Beam.angle(i)*180/pi), ' deg'], 'fontsize', szAxLabel2)
% %         xlabel('Time, t [s]', 'fontsize', szAxLabel1)
% %         ylabel('Amplitude (Real), p(t)', 'fontsize', szAxLabel1)
% %         h_fig=get(gcf,'CurrentAxes');
% %         set(h_fig, 'fontsize', szAxScale);
% %         axis 'square'; axis tight;
%         
%         figure('name', 'The received signal after basebanding')
%         plot(t, abs(stkuphi(rx.Ref,:,i)));
%         title(['Received signal after basebanding at ', num2str(tx.Beam.angle(i)*180/pi), ' deg'], 'fontsize', szAxLabel2)
%         xlabel('Time, t [s]', 'fontsize', szAxLabel1)
%         ylabel('Amplitude (Real), p(t)', 'fontsize', szAxLabel1)
%         h_fig=get(gcf,'CurrentAxes');
%         set(h_fig, 'fontsize', szAxScale);
%         axis 'square'; axis tight;
%         
%         figure('name', 'Magnitude spectrum of basebanded received signal')
%         plot(f/1000, abs(ftemp(rx.Ref,:)));
%         title('Magnitude spectrum of basebanded received signal', 'fontsize', szAxLabel2)
%         xlabel('Frequency, f [kHz]', 'fontsize', szAxLabel1)
%         ylabel('Magnitude, |P(f)|', 'fontsize', szAxLabel1)
%         h_fig=get(gcf,'CurrentAxes');
%         set(h_fig, 'fontsize', szAxScale);
%         axis 'square'; axis tight
%     end;
%     
% %     G = real(stuphi(:,:,i));
% %     colormap(hot(256));
% %     imagesc(t, 1:rx.M, G);
% %     title('The received signal after basebanding', 'fontsize', szAxLabel2)
% %     xlabel('Time, t [sec]', 'fontsize', szAxLabel1)
% %     ylabel('Receive aperture, u [m]', 'fontsize', szAxLabel1)
% %     axis 'square'; axis 'xy'; colorbar
% %     h_fig=get(gcf,'CurrentAxes');
% %     set(h_fig, 'fontsize', szAxScale);
% %     pause(5);
% end;

%% TEMPORAL DECIMATION

decimate.Ratio = 1;
decimate.Fs = fs/decimate.Ratio;
decimate.n = 2*floor(n/(2*decimate.Ratio));
decimate.dt = 1/decimate.Fs;

if(decimate.Ratio > 1)
    
    stuphiDecimated = zeros(size(stuphi,1), decimate.n, size(stuphi,3));

    v = floor( linspace(1, n, decimate.n) );

    for i = 1:tx.N
        for j = 1:rx.N
            for k = 1:decimate.n
                stuphiDecimated(j, k, i) = stkuphi( j, v(k), i );
            end;
        end;
    end;

else
    
    stuphiDecimated = stkuphi;
    
end;


%% RESAMPLE TO TARGET REGION

nR = 2*ceil(((X1 - X0)/c)/(2*decimate.dt));

stuphiTargetRegion = zeros(size(G,1), nR, size(G,3));

for i = 1:tx.N
    for j = 1:rx.N
        for k = 1:nR
            stuphiTargetRegion(j, k, i) = stuphiDecimated(j, k, i);
        end;
    end;
end;

newD = rx.d*(rx.N);
dtheta = asind(lambda_min/newD);

if(tx.N > 1)
    figure
    H = vol3d('cdata', abs(stuphiTargetRegion), 'texture', '3D', 'YData', (-rx.N/2:rx.N/2-1)*dtheta, 'ZData', tx.Beam.angle*180/pi, 'XData', (Ts+(0:nR-1)*decimate.dt)*c/2);
    view(3);
    xlabel('Range, ct [m]', 'fontsize', szAxLabel1)
    ylabel('Angle, \theta [deg]', 'fontsize', szAxLabel1)
    zlabel('Elevation, \phi [deg]', 'fontsize', szAxLabel1)
    axis tight; axis 'square';  colorbar; grid on
    h_fig=get(gcf,'CurrentAxes');
    set(h_fig, 'fontsize', szAxScale); 
    alphamap('rampup');
    alphamap(0.07.*alphamap);
else
    G = abs(stuphiTargetRegion);
    figure('name', 'Image')
    colormap(hot(256));
    % imagesc(G);
    imagesc((Ts+(0:nR-1)*decimate.dt)*c/2, (-rx.N/2:rx.N/2-1)*dtheta, G);
    title('Beam-time domain', 'fontsize', szAxLabel2)
    xlabel('Range, ct [m]', 'fontsize', szAxLabel1)
    ylabel('Angle, \theta [deg]', 'fontsize', szAxLabel1)
    axis 'square'; axis 'xy'; colorbar;
    h_fig=get(gcf,'CurrentAxes');
    set(h_fig, 'fontsize', szAxScale);
end;

