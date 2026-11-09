T = readtable('HRFlib_500Hz.csv');

T.Properties.VariableNames

figure;
plot(T.time, T.afni_spmg1)


Tbasis = readtable('HRFbasis_500Hz.csv');

figure; plot(Tbasis.time, [Tbasis.pc1, Tbasis.pc2, Tbasis.pc3])


V = [Tbasis.pc1(:), Tbasis.pc2(:), Tbasis.pc3(:)];
G = V' * V

Vn = V / sqrt(G(1,1));
Vn' * Vn   % should now be the identity to ~1e-16


figure;

figure; plot(Tbasis.time, [0.2*Tbasis.pc1 +  0.3*Tbasis.pc2  - 0.1* Tbasis.pc3])
