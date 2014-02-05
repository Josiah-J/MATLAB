
close all
clear all
clc;

% Nsd = 44; % What are these values? Probably the turns calculated using Phil's method
% N1d = 49; % Shading requirements
% N2d = 38;
% N3d = 29;

Sp1 = 0.30;
Sp2 = 0.68;
Sp3 = 1.0;

R1 = 20e3;%17880;%17000; % What are these values and where are they from? % Impedances of the transducers
R2 = 20e3;%16700;%16000;
R3 = 20e3;%12470;%13000;

Al = 1760e-9;%4420d-9; % Al-value

F = 420e3; % Frequency

R1_X1_ratio = 0.75;%2.1;%2.3;  % Where does 2.3 come from? - magnetic impedance to ceramic impedance ratio to avoid loading
% ceramics: empirical values which started at 10 and can go as low as 0.5 with good results.

X1 = R1*R1_X1_ratio;

L1 = X1/(2*pi*F); % Inductance that gives the given reactance

N1 = round(sqrt(L1/Al)); % What are these values? Assuming these are number of turns values
% N2 = round(N2d*N1/N1d);
% N3 = round(N3d*N1/N1d);

% X2 = R2*R1_X1_ratio;
% X3 = R3*R1_X1_ratio;
% L2 = X2/(2*pi*F);
% L3 = X3/(2*pi*F);
% N2 = round(sqrt(L2/Al));
% N3 = round(sqrt(L3/Al));


N2 = round(N1*(Sp2/Sp1)-N1);
N3 = round(N1*(Sp3/Sp1)-N1-N2);

L2 = Al*N2^2;
X2 = 2*pi*F*L2;
R2_X1_X2_ratio = (X1+X2)/R2;

L3 = Al*N3^2;
X3 = 2*pi*F*L3;
R3_X1_X2_X3_ratio = (X1+X2+X3)/R3;


Ns = 12;%10%round(Nsd*N1/N1d)

Rt3 = R3*Ns^2/(N1+N2+N3)^2;
Rt2 = R2*Ns^2/(N1+N2)^2;
Rt1 = R1*Ns^2/(N1)^2;

Rs = ((Rt1*Rt2)/(Rt1+Rt2))*Rt3/(((Rt1*Rt2)/(Rt1+Rt2))+Rt3);

Ls = Al*Ns^2;
Xs = 2*pi*F*Ls;

Lp1 = Al*(N1)^2;
Xp1 = 2*pi*F*Lp1;
nr1 = Xp1/Xs;

Lp2 = Al*(N1+N2)^2;
Xp2 = 2*pi*F*Lp2;
nr2 = Xp2/Xs;

Lp = Al*(N1+N2+N3)^2;
Xp = 2*pi*F*Lp;
nrp = Xp/Xs;
