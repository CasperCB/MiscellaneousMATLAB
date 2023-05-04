function [n k] = eps2nk(eps_1, eps_2)

n = sqrt((sqrt((eps_1.^2) + (eps_2.^2)) + eps_1)./2);
k = sqrt((sqrt((eps_1.^2) + (eps_2.^2)) - eps_1)./2); 

% n_complex = complex(n, k); 

end