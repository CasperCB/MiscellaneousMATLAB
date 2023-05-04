function [cond_AC] = cond_AC(cond_DC, wave, tau)

c = 2.99792 * 10^8; % m/s

freq = c*(wave*100);

cond_AC = cond_DC ./ (1 - 1i.*freq.*tau); 

end

