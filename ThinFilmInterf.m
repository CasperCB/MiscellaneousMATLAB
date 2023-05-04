% The goal of this script is to calculate the optical contrast for thin
% flakes on a single or multilayered substrate. In addition, it combines
% the reflectance spectrum with color-matching functions in order to
% predict the observed color of a thin flake given its thickness. 

% This program uses the addon jreftran_rt (v 2.1.0.0) by Shawn Divitt
% (2016) to do the transfer-matrix calculations: 
% https://tinyurl.com/rw77gqb

% The color matching functions used are based off of the cone fundamentals
% measured by Stockman & Sharpe (2000). And can be found at:
% http://www.cvrl.org/

% Clayton B. Casper (3/4/2020)

clear; clc;
%% add path dependencies 
addpath(genpath('D:\Atkin Research\Code & Simulation\MATLAB\Addons\')); 
addpath(genpath('D:\Atkin Research\Code & Simulation\MATLAB\Dielectric Functions\'));

%% define modeling parameters
thick = linspace(1, 20, 100); % thickness vector for flakes, nm
d_pts = length(thick); % no. of thickness points

wave = linspace(400, 800, 500); % wavelength vector, nm
wave_pts = length(wave); % no. of wavelength points

t0 = 0; % angle of incidence, rad
polarization = 1; % 0 for TE polarization, 1 for anything else

Norm = 1; % normalize reflectance to substrate? 0 for no, 1 for yes

% intialize matrices
R_mat = zeros(wave_pts, d_pts); % reflectance of flake at given wavelength and thickness
C_mat = R_mat; % optical contrast of flake at given wavelength and thickness
R_ref_mat = zeros(wave_pts, 1); % substrate reflectance spectrum
R_diff = R_mat; % difference spectrum of flakes from substrate

%% define sample parameters
% the sample system is defined as a stack, with 1 being the top layer. For
% bulk layers, define thickness as NaN. You can add as many layers as you
% want, but each need a thickness and refractive index.

% thickness of layers in nm
d1 = NaN;
d2 = thick; 
d3 = NaN; 
% d4 = NaN;

% complex refractive indices of layers
n1 = 1; % air
n2 = 1.55;
n3 = GoldE_vis(wave); 
% n4 = 4; 

%% calculate reflectance w/ transfer-matrix method

for i = 1:wave_pts % loop over all wavelengths
    
    % build refractive index vector
    n = [n1 n2 n3(i)];
    n_ref = [n1 n1 n3(i)]; % for substrate reflectance (treat flake as air)

    for k = 1:d_pts % loop over all flake thicknesses
        
    % build thickness vector
    d = [d1 d2(k) d3];    
    
        if k == 1 % calculate substrate reflectance (only need to do this once)
            [r_ref, t_ref, R_ref, T_ref, A_ref] = jreftran_rt(wave(i), d, n_ref, t0, 1); 
            R_ref_mat(i) = R_ref; 
        end
    
    % calculate reflectance of system at given wavelength and thickness
    [r, t, R, T, A] = jreftran_rt(wave(i), d, n, t0, 1);
    
    % calculate optial contrast from substrate
    C = ((R_ref_mat(i) - R)./ R_ref_mat(i)).*100; 

    % store values 
    R_mat(i,k) = R;
    C_mat(i,k) = abs(C); 

    % reflectance value after substrate subtraction
    R_diff(i,k) = R_mat(i,k) - R_ref_mat(i); 
        
    end
    
    % simulate white light optical contrast by averaging
    C_avg = mean(C_mat, 1); 
    
    % if normalizing to substrate, use difference spectrum and shift to
    % positive values
    if Norm == 1
        R_mat(i,:) = R_diff(i,:) - min(R_diff(i,:)); 
    end

    
end

clear i k 

%% predict flake colors w/ RGB color matching functions

% intialize matrix to hold RGB values for a given thickness
col_mat = zeros(d_pts, 3); 

% read CIE functions and define sensitivity vectors
CIE = csvread('CIE_2.csv');
wave_CIE = CIE(:, 1); 
R_CIE = CIE(:, 2); 
G_CIE = CIE(:, 3);
B_CIE = CIE(:, 4); 


for k = 1:d_pts % loop over all flake thicknesses
    
    for i = 1:wave_pts % loop over all wavelengths
    
    % grab sensitivity for a given wavelength
    R_i  = interp1(wave_CIE, R_CIE, wave(i)); 
    G_i = interp1(wave_CIE, G_CIE, wave(i)); 
    B_i = interp1(wave_CIE, B_CIE, wave(i)); 
    
    % multiply reflectance value of system by sensitivity
    X(i) = R_mat(i,k).*R_i;
    Y(i) = R_mat(i,k).*G_i;
    Z(i) = R_mat(i,k).*B_i; 
    
    % do the same for the reference substrate
        if k == 1 
            X_ref(i) = R_ref_mat(i).*R_i;
            Y_ref(i) = R_ref_mat(i).*G_i;
            Z_ref(i) = R_ref_mat(i).*B_i; 
        end
        
    end
    
    % get RGB by completing convolution via summing products and normalizing
    tot = sum(X)+sum(Y)+sum(Z); 
    col_mat(k, 1) = sum(X)/tot; 
    col_mat(k, 2) = sum(Y)/tot;
    col_mat(k, 3) = sum(Z)/tot; 
    
    % do the same for reference substrate
    tot_ref = sum(X_ref)+sum(Y_ref)+sum(Z_ref);
    col_ref(1,1) = sum(X_ref)/tot_ref;
    col_ref(1,2) = sum(Y_ref)/tot_ref; 
    col_ref(1,3) = sum(Z_ref)/tot_ref; 
    
end

%% plot
% average optical contrast w/ each point colored
hold on
% h = scatter(d2, C_avg, 'filled'); 
% c = h.CData;
% c = repmat(c, [d_pts, 1]); 
% c = col_mat.*1;
% h.CData = c.*1.0; % artificially brighten RGB values for clarity
xlabel('\bf flake thickness \rm (nm)'); 
ylabel('\bf optical contrast \rm (%)'); 
set(gca, 'fontsize', 16);
box on; 
