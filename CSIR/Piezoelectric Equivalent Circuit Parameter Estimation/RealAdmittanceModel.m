function yHatReal = RealAdmittanceModel(beta, f)

Rm = beta(1);
Lm = beta(2);
Cm = beta(3);
C0 = beta(4);

w = 2*pi*f;

temp1 = w*Cm - (w*C0).*(w.^2*Lm*Cm - 1) + 1j*w.^2*Cm*Rm*C0;
temp2 = Rm*w*Cm + 1j*(w.^2*Lm*Cm - 1);

yHatReal = real(temp1./temp2);
