%% PHASED-TRANSDUCER ARRAY SIMULATION
% This code simulates beam-steering using a phased ultrasound transducer
% array. 

%% OPTIONS
close all
clear all
clc

%% PULSE PARAMETERS

c = 1500; % Speed of sound in water
fc = 200e3; % Centre frequency in Hz
wc = 2*pi*fc; % Centre frequency in rad/s
k = wc/c; % Wavenumber at the centre frequency
lambda_c = c/fc; % Wavelength at the centre frequency

%% ARRAY PARAMETERS

N = 32; % Number of transmit elements in the array
d = 25e-3; % Inter-element spacing
a = 25e-3; % Element width
D = d*(N-1)+a; % Aperture size

theta_s = 30; % Steering angle in degrees

%% SIMULATION REGION PARAMETERS

X0 = 5;
Y0 = 10;

dX = 2*X0/300;
dY = Y0/300;

nx = 2*ceil(2*X0/(2*dX));
ny = 2*ceil(Y0/(2*dY));

x = dX*(-nx/2:nx/2-1);
y = dY*(0:ny-1);


%% PRESSURE DISTRIBUTION

p = zeros(ny, nx);
t = 0.05;
for i = -nx/2:nx/2-1
    for j = 1:ny
        r = (x(i+nx/2+1).^2 + y(j).^2).^0.5;
%         t = r/c;
        
        theta = atand(x(i+nx/2+1)/y(j));
        
        temp1 = ( sind(pi*a*sind(theta)/lambda_c) )/( pi*sind(theta)/lambda_c );
        temp2 = sind(pi*d*(sind(theta_s) - sind(theta))*N/lambda_c);
        temp3 = sind(pi*d*(sind(theta_s) - sind(theta))/lambda_c);
        
        arg1 = pi*a*sind(theta)/lambda_c;
        arg2 = pi*d*(sind(theta_s) - sind(theta))*(N-1)/lambda_c;
        arg3 = 2*pi*fc*t - k*r;        
        
        p(j, i+nx/2+1) = temp1*(temp2/temp3)*exp(-1j*arg1)*exp(-1j*arg2)*exp(1j*arg3);        
    end;
end;

G=abs(p);
xg=max(max(G)); ng=min(min(G)); cg=255/(xg-ng);

figure('name', 'Pressure distribution')
colormap(hot(256))
imagesc(x, y, G);
% title('Matched filtered signal variation with angle', 'fontsize', szAxLabel)
% xlabel('Range, ct [m]', 'fontsize', szAxLabel)
% ylabel('Angle, \theta [deg]','fontsize', szAxLabel)
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);
axis('square'); colorbar
