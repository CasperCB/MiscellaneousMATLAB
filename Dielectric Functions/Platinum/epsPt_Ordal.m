function [epsComplex] = epsPt_Ordal(wave)

% parameters taken from Ordal 1985, Appl. Opt.

wt = 5.58e2; % cm^-1
wp = 4.15e4; % cm^-1

eps1 = 1 - (wp.^2)./((wave.^2)+(wt.^2)); 
eps2 = (wt.*(wp.^2))./(wave.*((wave.^2)+(wt.^2))); 

epsComplex = complex(eps1, eps2); 

end