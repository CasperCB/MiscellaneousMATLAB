% Two layer extention to the Finite Dipole Model, adapted from Hauer's
% Thesis (2015).
%
% Top layer: air (eps=1). epsSample1: top layer, thickness d. epsSample2:
% bottom layer, thickness undefined.
%
% Earl Ritchie  10/2019.
function [Alpha_eff] = ML_FDM(Trad, L, height, epsSamp1, d, epsSamp2, wavenumber)

beta12 = (1-epsSamp1) ./ (1+epsSamp1);                  % Layer 1 Reflection Coefficient
beta23 = (epsSamp1-epsSamp2) ./ (epsSamp1+epsSamp2);    % Layer 2 Reflection Coefficient

q = 1/Trad;                                % Momentum of tip, m^-1
angWavenum = 2 .* pi .* wavenumber .* 100; % Angular wavenumber of illumination, m^-1
k = sqrt(angWavenum.^2 - q.^2);
B_eff = -(beta12 + beta23.*exp(1i.*2.*k.*d)) ./ (1 + beta12.*beta23.*exp(1i.*2.*k.*d));

% induced charges
W0 = 1.31*Trad;
W1 = 0.5*Trad; % position of charge [Q0, Q1] with respect to radius, m
g = 0.9*exp(0.06*1i); % contributing portion of near-field induced charge
% g = (1 - 0.7) * exp(-1 * height / (Trad/2)) + 0.7*exp(0.06*1i);
% g = 0.98*exp(0.08*1i);

% finite-dipole geometric factors f0 and f1, formula broken into two parts
f0_1 = (g-((Trad+2*height+W0)./(2*L)));
f0_2 = (log((4*L)./(Trad+(4*height)+(2*W0))))./(log((4*L)./Trad)); 
f0 = f0_1.*f0_2;

f1_1 = (g-((Trad+2*height+W1)./(2*L)));
f1_2 = (log((4*L)./(Trad+(4*height)+(2*W1))))./(log((4*L)./Trad)); 
f1 = f1_1.*f1_2;

Alpha_eff = (0.5.*((B_eff.*f0)./(1-(B_eff.*f1))))+1; % complex polarizability of tip-sample system
% Alpha_eff = (1+ B_eff) .* Alpha_eff;

end