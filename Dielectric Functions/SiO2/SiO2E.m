function [epsComplex] = SiO2E(wavenumber)

SiO2ExpData = importdata('Kischkat_SiO2.txt', '\t'); % From Kischkat Simulation
SiO2Wavenumber = 10000 ./ SiO2ExpData.data(:, 1);

nSiO2=interp1(SiO2Wavenumber, SiO2ExpData.data(:, 2), wavenumber);
kSiO2=interp1(SiO2Wavenumber, SiO2ExpData.data(:, 3), wavenumber);

eps1 = nSiO2.^2 - kSiO2.^2;
eps2 = 2 .* nSiO2 .* kSiO2;

epsComplex = eps1 + 1i.*eps2;

end