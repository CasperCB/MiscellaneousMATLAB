function [epsComplex] = STOE(wavenumber)

STOExpData = importdata('STO.xlsx'); % from McArdle 2020
STOWavenumber = STOExpData(:,1);

nSTO=interp1(STOWavenumber, STOExpData(:, 2), wavenumber);
kSTO=interp1(STOWavenumber, STOExpData(:, 3), wavenumber);

eps1 = nSTO.^2 - kSTO.^2;
eps2 = 2 .* nSTO .* kSTO;

epsComplex = eps1 + 1i.*eps2;

end