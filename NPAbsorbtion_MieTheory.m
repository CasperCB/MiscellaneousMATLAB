% Clayton Casper (08/26/16) - code for calculating the absorption cross
% section of a nanoparticle using Mie Theory. This code is an
% adaptation of one provided by Guillaume Baffou (2016)

clear; clc; 
%% Define parameters & optical constants
wavelength = 10.63*(10^-6); 
wavenumber = (1./wavelength)./(100); 

% define material dielectric function
eps = epsSiC_Taubner2004(wavenumber, 0); 

ref_index = sqrt ( ( sqrt( (real(eps).^2) + (imag(eps).^2) ) + real(eps)) ./ 2 ); 
ext_coeff = sqrt ( ( sqrt( (real(eps).^2) + (imag(eps).^2) ) - real(eps)) ./ 2 ); 

n_mat = complex(ref_index, ext_coeff); 

n_m = 1; % refractive index of the medium 
r0 = 100*(10^-9); % particle radius in m

%% calculation

% create empty arrays
numPoints = length(wavelength);
Qext_mat= zeros(numPoints(1), 1);
Qsca_mat = zeros(numPoints(1), 1);


for s = 1:numPoints
    
    m = n_mat(1,s) / n_m; % dimensionless parameter: ratio of indicies of refraction
    k = (2*pi*n_m) / wavelength(1,s); % wavevector of light in the medium

    x = k*r0; % dimensionless parameter related to propagation of light outside the sphere
    z = m*x;  % dimensionless parameter related to propagation of light inside the sphere
    
    N = round(2 + x + 4 *x^(1/3)); % estimate of maximum N-pole contribution as suggested by Bohren & Huffman (1983)
    

    j = (1:N); % define index for Bessel functions

    % define spherical Bessel functions of 1st & 2nd kind + Hankel function
    sqr = sqrt(pi*x/(2)); 
    sqr_m = sqrt(pi*x/(2));
    phi = sqr.*besselj(j+0.5, x); 
    xi = sqr.*(besselj(j+0.5, x) + 1i*bessely(j+0.5, x)); 
    phi_m = sqr_m.*besselj(j+0.5, z);
    phi_l = [sin(x), phi(1:N-1)]; 
    phi_lm = [sin(z), phi_m(1:N-1)]; 
    y = sqr*bessely(j+0.5, x); 
    yl = [-cos(x), y(1:N-1)]; 

    % define derivatives
    phi_prime = (phi_l-(j/x).*phi); 
    phi_mprime = (phi_lm-(j/z).*phi_m);
    xi_prime = (phi_l + 1i*yl)-(j/x).*xi; 

    % calculate extinction and scattering cross sections
    a_j = (m*phi_m.*phi_prime-phi.*phi_mprime) ./ (m*phi_m.*xi_prime - xi.*phi_mprime); 
    b_j = (phi_m.*phi_prime - m*phi.*phi_mprime) ./ (phi_m.*xi_prime - m*xi.*phi_mprime); 
    Qext_mat(s) = (2*pi/(k^2))*sum((2*j+1).*real(a_j + b_j)); 
    Qsca_mat(s) = (2*pi/(k^2))*sum((2*j+1).*((abs(a_j).*abs(a_j))+(abs(b_j).*abs(b_j))));
    
end

% calculate absorption cross section
Qabs_mat1 = Qext_mat - Qsca_mat; 

% %% plot normalized results
% 
% plot(wavelength.*(10^9), Qsca_mat, '-o',  'LineWidth', 2);
% 
% hold on
% 
% % plot(wavelength, Qsca_mat1, '-o', 'LineWidth', 2);
% 
% % ax = gca;
% % ax.XLim = [0.4*10^-6, 2.0*10^-6]; 
% 
% xlabel('\bf Wavelength \rm (\mum)'); 
% ylabel('\bf Cross section \rm (m^2)'); 
% title('Absorbance cross section for 200 nm Au NP'); 
% % legend('1x10^1^9 cm^3', '2x10^1^9 cm^3', '3x10^1^9 cm^3', '4x10^1^9 cm^3', '5x10^1^9 cm^3', '6x10^1^9 cm^3', '7x10^1^9 cm^3', '8x10^1^9 cm^3', '9x10^1^9 cm^3');


   










