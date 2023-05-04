function [epsComplex] = Glass_soda_E(wavenumber)

wl = (1./wavenumber).*1e4; % convert from cm-1 to um

GlassExpData = csvread('Rubin_1985.csv'); % From Zhang 2020
GlassWavenumber = 10000 ./ GlassExpData(:, 1);

nGlass=interp1(GlassWavenumber, GlassExpData(:, 2), wavenumber);
kGlass=interp1(GlassWavenumber, GlassExpData(:, 3), wavenumber);

eps1 = nGlass.^2 - kGlass.^2;
eps2 = 2 .* nGlass .* kGlass;

epsComplex = eps1 + 1i.*eps2;

end