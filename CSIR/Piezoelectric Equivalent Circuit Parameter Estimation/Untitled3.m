
clear all
close all
clc

% Plot properties
szAxScale = 18;
szAxLabel = 18;
szAxLabel2 = 14;

% File path
impedanceFileDir= 'I:\M&Mtek\SST\086MG\HACS07X - Pakistan\3. Project Outputs\THABISO';
impedanceFileName = 'Module 21 element 4 3-4';

%% IMPORT DATA, SET VECTORS AND PLOT IMPEDANCE AND ADMITTANCE VALUES

impedanceFilePath = [impedanceFileDir '\' impedanceFileName '.csv'];
dataset = importImpedanceFile(impedanceFilePath);

f = dataset(:,1);
z = dataset(:,2).*exp(1j*pi*dataset(:,4)/180);
y = dataset(:,3).*exp(1j*pi*dataset(:,4)/180);

% Impedance Plots
figure('name','Impedance Plots')
subplot(121); plot(f, real(z),'LineWidth',2);
hold on
plot(f, imag(z), '--r','LineWidth',2);
xlabel('Frequency, f [Hz]', 'fontsize', szAxLabel)
ylabel('Impedance, Z [Ohms]', 'fontsize', szAxLabel)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
legend('Re\{Z\}','Im\{Z\}');
axis 'square'
grid on
set(gca, 'GridLineStyle', '-');
grid(gca,'minor')

subplot(122); plot(f,abs(z), 'LineWidth',2)
xlabel('Frequency, f [Hz]', 'fontsize', szAxLabel)
ylabel('Magnitude Impedance, |Z| [Ohms]', 'fontsize', szAxLabel)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'
grid on
set(gca, 'GridLineStyle', '-');
grid(gca,'minor')

figure('name','Measured magnitude impedance of piezoelectric')
plot(f,abs(z), 'LineWidth',2)
title('Magnitude impedance (measured)', 'fontsize', szAxLabel2)
xlabel('Frequency, f [Hz]', 'fontsize', szAxLabel)
ylabel('Magnitude impedance, |Z| [Ohms]', 'fontsize', szAxLabel)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'
grid on
set(gca, 'GridLineStyle', '-');
grid(gca,'minor')
print('-dpng', '-r300', './Measured impedance')

% Admittance Plots
% figure('name','Admittance Plots')
% subplot(121); plot(f, real(y),'LineWidth',2);
% hold on
% plot(f, imag(y), '--r','LineWidth',2);
% xlabel('Frequency, f [Hz]', 'fontsize', szAxLabel)
% ylabel('Admittance, Y [S]', 'fontsize', szAxLabel)
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);
% legend('Re\{Y\}','Im\{Y\}');
% axis 'square'
% grid on
% set(gca, 'GridLineStyle', '-');
% grid(gca,'minor')
% 
% subplot(122); plot(f,abs(y), 'LineWidth',2)
% xlabel('Frequency, f [Hz]', 'fontsize', szAxLabel)
% ylabel('Magnitude Admittance, |Y| [S]', 'fontsize', szAxLabel)
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);
% axis 'square'
% grid on
% set(gca, 'GridLineStyle', '-');
% grid(gca,'minor')

% Impedance Hodograph
% figure('name','Impedance Hodograph')
% plot(real(z), imag(z),'LineWidth',2);
% xlabel('Resistance, Re\{Z\} [Ohm]', 'fontsize', szAxLabel)
% ylabel('Reactance, Im\{Z\} [Ohm]', 'fontsize', szAxLabel)
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);
% axis 'square'
% grid on
% set(gca, 'GridLineStyle', '-');
% grid(gca,'minor')

%% EXTRACT IMPEDANCE AND ADMITTANCE VALUES AT A PARTICULAR FREQUENCY

fs = f(abs(y) ==  max(abs(y)));
ws = 2*pi*fs;

fp = f(abs(y) == min(abs(y)));
wp = 2*pi*fp;

fprintf('Impedance at 420kHz = %.2f %.2fj Ohms\n', real(z(f==420e3)), imag(z(f==420e3)));
fprintf('Magnitude impedance at 420kHz = %.2f Ohms\n\n', abs(z(f==420e3)));
fprintf('Series resonance frequency = %.2f kHz\nParallel resonance frequency = %.2f kHz\n', fs/1000, fp/1000);

%% COMPUTE TRANSDUCER EQUIVALENT CIRCUIT VALUES AT RESONANCE

Yfs = y(f==fs); % Admittance at series resonance

Gs = real(Yfs);
Bs = imag(Yfs);

Rm = 1/Gs; % Motional resistance at series resonance
C0 = Bs/ws; % Clamp capacitance at series resonance

Cm = C0*(wp^2/ws^2 - 1); % Motional capacitance
Lm = 1/(ws^2 * Cm); % Motional Inductance

Qm = (1/Rm)*(sqrt(Lm/Cm));

fp_computed = sqrt(fs^2*(Cm/C0+1));

fprintf('Lm = %0.2f mH\nRm = %.2f Ohms\nCm = %.2f pF\nC0 = %.2f pF\n',Lm*1000, Rm, Cm*1e12, C0*1e12);

%% COMPUTE TRANSDUCER ADMITTANCE OVER A RANGE OF FREQUENCY FROM THE CIRCUIT VALUES AT RESONANCE

w = 2*pi*f;

temp1 = w*Cm - (w*C0).*(w.^2*Lm*Cm - 1) + 1j*w.^2*Cm*Rm*C0;
temp2 = Rm*w*Cm + 1j*(w.^2*Lm*Cm - 1);

Y0_computed = temp1./temp2;
% figure('name','Admittance Plots using calculated equivalent circuit values')
% subplot(121); plot(f, real(Y0_computed),'LineWidth',2);
% hold on
% plot(f, imag(Y0_computed), '--r','LineWidth',2);
% xlabel('Frequency, f [Hz]', 'fontsize', szAxLabel)
% ylabel('Computed admittance, Y [S]', 'fontsize', szAxLabel)
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);
% legend('Re\{Y\}','Im\{Y\}');
% axis 'square'
% grid on
% set(gca, 'GridLineStyle', '-');
% grid(gca,'minor')
% 
% subplot(122); plot(f,abs(Y0_computed), 'LineWidth',2)
% xlabel('Frequency, f [Hz]', 'fontsize', szAxLabel)
% ylabel('Computed magnitude Admittance, |Y| [S]', 'fontsize', szAxLabel)
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);
% axis 'square'
% grid on
% set(gca, 'GridLineStyle', '-');
% grid(gca,'minor')

%% NON-LINEAR REGRESSION TO ESTIMATE ACTUAL EQUIVALENT CIRCUIT PARAMETERS

[beta, Resid, Jaco, MSC, ErrorMod] = nlinfit(f, abs(y), @admittanceModel, [Rm Lm Cm C0]);

Rm_hat = beta(1);
Lm_hat = beta(2);
Cm_hat = beta(3);
C0_hat = beta(4);

fs_hat = 1/(2*pi*sqrt(Lm_hat*Cm_hat));
fp_hat = sqrt(fs_hat^2*(Cm_hat/C0_hat+1));

fprintf('After nonlinear regression:\nEstimated series frequency = %.2f kHz\nEstimated parallel frequency = %.2f kHz\n',fs_hat/1000, fp_hat/1000);
fprintf('Lm = %.3e H\nRm = %.3e Ohms\nCm = %.3e F\nC0 = %.3e F\n',Lm_hat, Rm_hat, Cm_hat, C0_hat);


%% COMPUTE TRANSDUCER ADMITTANCE OVER A RANGE OF FREQUENCY FROM THE ESTIMATED CIRCUIT VALUES

temp1 = w*Cm_hat - (w*C0_hat).*(w.^2*Lm_hat*Cm_hat - 1) + 1j*w.^2*Cm_hat*Rm_hat*C0_hat;
temp2 = w*Rm_hat*Cm_hat + 1j*(w.^2*Lm_hat*Cm_hat - 1);

Y0_estimated = temp1./temp2;
% figure('name','Admittance Plots using regressed equivalent circuit values')
% subplot(121); plot(f, real(Y0_estimated),'LineWidth',2);
% hold on
% plot(f, imag(Y0_estimated), '--r','LineWidth',2);
% xlabel('Frequency, f [Hz]', 'fontsize', szAxLabel)
% ylabel('Computed admittance, Y [S]', 'fontsize', szAxLabel)
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);
% legend('Re\{Y\}','Im\{Y\}');
% axis 'square'
% grid on
% set(gca, 'GridLineStyle', '-');
% grid(gca,'minor')
% 
% subplot(122); plot(f,abs(Y0_estimated), 'LineWidth',2)
% xlabel('Frequency, f [Hz]', 'fontsize', szAxLabel)
% ylabel('Computed magnitude Admittance, |Y| [S]', 'fontsize', szAxLabel)
% h_fig=get(gcf,'CurrentAxes');
% set(h_fig, 'fontsize', szAxScale);
% axis 'square'
% grid on
% set(gca, 'GridLineStyle', '-');
% grid(gca,'minor')

Z0_estimated = 1./Y0_estimated;
figure('name','Impedance Plots using regressed values')
subplot(121); plot(f, real(Z0_estimated),'LineWidth',2);
hold on
plot(f, imag(z), '--r','LineWidth',2);
xlabel('Frequency, f [Hz]', 'fontsize', szAxLabel)
ylabel('Impedance, Z [Ohms]', 'fontsize', szAxLabel)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
legend('Re\{Z\}','Im\{Z\}');
axis 'square'
grid on
set(gca, 'GridLineStyle', '-');
grid(gca,'minor')

subplot(122); plot(f, abs(Z0_estimated), 'LineWidth',2)
xlabel('Frequency, f [Hz]', 'fontsize', szAxLabel)
ylabel('Magnitude Impedance, |Z| [Ohms]', 'fontsize', szAxLabel)
h_fig=get(gcf,'CurrentAxes');
set(h_fig, 'fontsize', szAxScale);
axis 'square'
grid on
set(gca, 'GridLineStyle', '-');
grid(gca,'minor')
