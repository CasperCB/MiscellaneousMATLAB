function [eps] = nk2eps(n, k)
    eps1 = n.^2 - k.^2;
    eps2 = 2*n.*k; 
    
    eps = complex(eps1, eps2);

end