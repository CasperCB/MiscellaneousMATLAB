function [n] = Altermatt_ii(Nd)
% function computes incomplete ionization of P-doped Si based on Altermatt
% et. al (2006). It uses the curve plotted in Fig. 6, which was then
% digitized.

% data range is from 6e17 to 1e21 cm-3

% both n & Nd are in cm-3

% Clayton B. Casper (03/05/2020)

ref = csvread('Altermatt_ii.csv'); 

n = interp1(ref(:,1), ref(:,2), Nd); 

end

