function [epsComplex] = PMMA_E(wavenumber)

wl = (1./wavenumber).*1e4; % convert from cm-1 to um

PMMAExpData = csvread('PMMA_Zhang_Tomson.csv'); % From Zhang 2020
PMMAWavenumber = 10000 ./ PMMAExpData(:, 1);

nPMMA=interp1(PMMAWavenumber, PMMAExpData(:, 2), wavenumber);
kPMMA=interp1(PMMAWavenumber, PMMAExpData(:, 3), wavenumber);

eps1 = nPMMA.^2 - kPMMA.^2;
eps2 = 2 .* nPMMA .* kPMMA;

epsComplex = eps1 + 1i.*eps2;

end