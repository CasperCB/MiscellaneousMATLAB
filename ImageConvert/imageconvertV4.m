% The purpose of this script is to convert near-field images/linecuts of a
% doped semiconductor into images of electrical properties (density,
% mobility, conductivity). These values are derived by fitting values of
% amplitude & phase on every pixel to a parameterized (density and
% mobility) Finite-Dipole Model (Clayton B. Casper, 10/30/19)

% This script calls on several other functions which will need to be
% located somewhere on the file path MATLAB uses:
% 
%   nrrdread = read near-field images prepared as .nrrd files by Gwyddion
%   Alpha_FDM = calculates Finite-Dipole Model with rp correction
%   ML_FDM = calculates Multilayer Finite-Dipole Model
%   demod_Sim = performs FFT on a signal and grabs harmonics of a frequency
%   nSiE_nomob/pSiE_nomob = custom eps functions for a parameterized 
%   density and mobility (does not need to be Si)

% Because many Fourier transforms are computed, this script uses 
% parallelization and requires the Parallel Computing Toolbox from 
% MathWorks. If you don't want to use it, replace parfor loops with for 
% loops. 

% Last updated (V4) 02/06/2020: added outputs: parameters.txt and .gwy file


%% system stuff
clear; clc;
workers = 4; % no. of compute nodes
tic 
f = waitbar(0, 'Importing data...'); 

%% paths for dependencies
% if your name is not Clayton, replace with your own paths
addpath(genpath('D:\Atkin Research\Code & Simulation\MATLAB\Tip Dipole Models\')); 
addpath(genpath('D:\Atkin Research\Code & Simulation\MATLAB\Dielectric Functions\'));
addpath(genpath('D:\Atkin Research\Code & Simulation\MATLAB\Data Workup\')); 

% inputpath = 'C:\Users\Earl\Desktop\Conductivity Workup\Nanowires\G152\12-12-19_Andrew\';
% inputpath = 'C:\Users\Earl\Desktop\Conductivity Workup\Nanowires\DHJ 223\2019-08-02\';
inputpath = 'D:\Atkin Research\Data\Si-NWs\Andrew Scans\G152\1-10_G152\stitch\';
addpath(genpath(inputpath));

ScanName = 'results.gwy';
data_mode = 2; % 1 for line cut, 2 for image
if data_mode == 2
    [mask, meta] = nrrdread('mask.nrrd');
    amp = nrrdread('amp.nrrd');
    phi = nrrdread('phi.nrrd'); 
    phi = phi  .* (180/pi); % Comment out if data already in degrees
elseif data_mode == 1
    amp_data = importdata('p-amp.txt');
    phi_data = importdata('p-phase.txt');
end 

%% constants
e = 1.602*10^-19; % C, elementary charge
m0 = 9.109*10^-31; % kg, electron rest mass
c = 2.99792 * 10^8; % m/s, speed of light in vaccuum 

%% modeling parameters & options
L_n = 500; % no. of carrier densities tested
L_mu = 500; % no. of carrier mobilities tested

n_param = logspace(18, 21, L_n); % cm^-3, active carrier density 
mu_param = linspace(1, 300, L_mu); % (cm^2)/(V*s), carrier mobility

FDM_mode = 1; % 0 = use regular FDM, 1 = use ML-FDM
Sinorm = 1; % 0 = data is normalized to Au, 1 = data is normalized to i-Si
type = 1; % 0 = p-type Si, 1 = n-type Si
gwy = 1; % 1 = save a .gwy file with results, 0 = no save

%% experimental parameters
harm_ex = 3; % tapping harmonic
wave = 940; % cm^-1, incident wavenumber
A_tap = 60*10^-9; % tapping amplitude, m
d = 3*10^-9; % nm, SiO2 thickness
angle = pi/3;
ax_start = str2num(meta.axismins); ax_end = str2num(meta.axismaxs);
x_range = (ax_end(1) - ax_start(1)) * 10^6; % x-range of input in microns (only needed for .gwy output)
y_range = (ax_end(2) - ax_start(2)) * 10^6; % y-range of input in microns (only needed for .gwy output)

if Sinorm == 0
    epsSub = GoldE(wave); % evaporated gold, Olmon
elseif Sinorm == 1
    epsSub = 11.68; % i-Si permittivity 
end

%% sampling parameters for Fourier transform
N = 12500; % number of points, this will drastically change compute time
N_FFT = 2^nextpow2(N); % change number of points to a power of 2 for faster FT
T = 6.9 * (10^-3); % sampling time, s 
t = (0:(N_FFT-1))/N_FFT; % time vector, s
t = t .* T; 

%% tip parameters
% geometric
r = 25*10^-9; % tip radius, m
L = 600*10^-9; % length of tip spheroid, m

% tapping
freq_tap = 250*10^3; % tapping frequency, Hz
H = A_tap*(1+cos(2*pi*freq_tap*t)); % time-dependent tip height, m

%% data import & export
% folder for .figs & .pngs (will overwrite!) 
savepath = [inputpath 'Processed\'];
% savepath = 'C:\Users\Earl\Desktop\Conductivity Workup\Infineon\Andrew\Processed\';

% need 3 image inputs: amplitude, phase, and mask
% make sure mask parity is correct (0-valued pixels are computed)
% if normalizing to i-Si, need an exp. value for i-Si (A_i.txt)

if data_mode == 2
%     mask = nrrdread('mask.nrrd');
%     amp = nrrdread('amp.nrrd');
%     phi = nrrdread('phi.nrrd'); 
    
    if Sinorm == 0
        A_i = 1;
    elseif Sinorm == 1
        A_i = 1;
%         A_i = importdata('A_i.txt');
    end

    x_dim = length(mask(1,:)); y_dim = length(mask(:,1)); 

elseif data_mode == 1
    
%     amp_data = importdata('amp.dat');
%     phi_data = importdata('phi.dat'); 
    
    x = amp_data.data(:,1).*(10^6); 
    amp = amp_data.data(:,2); 
    phi = phi_data.data(:,2) .* (180/pi);
    
%     if Sinorm == 0
        A_i = 1;
%     elseif Sinorm == 1
%         A_i = importdata('A_i.txt');
%     end

   x_dim = 1; 
   y_dim = length(amp); 

end

% intialize images to be calculated
n = zeros(y_dim, x_dim);
mu = zeros(y_dim, x_dim);
err = zeros(y_dim, x_dim); 

%% calculate FDM values for A & phi for parameter space
% calculate complex polarizability for substrate reference
if FDM_mode == 0
    Alpha_ref = Alpha_FDM(r, L, H, epsSub, epsSub, angle);
    
% if ML, account for oxide if normalizing to i-Si
elseif FDM_mode == 1
    if Sinorm == 0
        Alpha_ref = ML_FDM(r, L, H, SiO2E(wave), 0, epsSub, wave);
%         Alpha_ref = ML_FDM_refl(r, L, H, SiO2E(wave), 0, epsSub, wave); 
    elseif Sinorm == 1
        Alpha_ref = ML_FDM(r, L, H, SiO2E(wave), d, epsSub, wave); 
%         Alpha_ref = ML_FDM_refl(r, L, H, SiO2E(wave), d, epsSub, wave); 
    end
end

% Apply FFT & grab releveant harmonics for reference substrate
harmonics_ref = demod_Sim(N_FFT, T, freq_tap, Alpha_ref, 4); 

% calculate reference amplitude & phase
amp_ref = abs(harmonics_ref(harm_ex));
phi_ref = atan2(imag(harmonics_ref(harm_ex)), real(harmonics_ref(harm_ex))).*(180/pi); 

% initialize matrices for parameter-dependent FDM
amps_FDM = zeros(L_n, L_mu);
phis_FDM = amps_FDM;

% main calculation loop, calculate FDM values for n & mu parameter set
for i = 1:L_n
    
    waitbar((i/L_n)/2, f, 'Calculating FDM...');
    
    parfor (j = 1:L_mu, workers)
        
        % initialize matrices for polarizability calculation
        Alphaeff = zeros(N_FFT, 1); 
        harmonics = zeros(1, 4); 
       
        % calculate permittivity of semiconductor for a given n & mu
        if type == 1
            epsSamp = nSiE_noMob(wave, n_param(i)*(10^6), mu_param(j)*(10^-4)); 
        elseif type == 0
            epsSamp = pSiE_noMob(wave, n_param(i)*(10^6), mu_param(j)*(10^-4)); 
        end

        % complex polarizability & associated harmonics 
        if FDM_mode == 0
            Alphaeff(:) = Alpha_FDM(r, L, H, epsSamp, epsSub, angle);
        elseif FDM_mode == 1
            Alphaeff(:) = ML_FDM(r, L, H, SiO2E(wave), d, epsSamp, wave);
%             Alphaeff(:) = ML_FDM_refl(r, L, H, SiO2E(wave), d, epsSamp, wave); 
        end
        
        harmonics(:) = demod_Sim(N_FFT, T, freq_tap, Alphaeff(:), 4);  
        
        % expected amplitude by square modulus & normalize to reference
        expect_amp = abs(harmonics(harm_ex))./amp_ref; 
        expect_phi = (atan2(imag(harmonics(harm_ex)), real(harmonics(harm_ex))).* ...
            (180/pi)) - phi_ref; 
%         expect_phi = (atan2(imag(harmonics(harm_ex)), real(harmonics(harm_ex))))...
%             - phi_ref;
        
        amps_FDM(i,j) = expect_amp;
        phis_FDM(i,j) = expect_phi; 
        
    end
end

%% determine best-fit n & mu

% loop over all pixels
for k = 1:y_dim   
    
    waitbar(0.5 + (k/y_dim)/2, f, 'Building image...');
    
    for l = 1:x_dim
        
        if  data_mode == 1 || mask(k,l) == 0 % only calculate on relevant pixels
            
            % grab amplitude & phase of pixel
            A = amp(k,l);
            P = phi(k,l); 
            
            % scale experimental amplitude to reference
            A_s = A*(1/A_i); 
            
            % initialize error matrix
            err_param = 1000.*ones(L_n, L_mu); 
            
            % calculate % error between pixel and parameter space
            for s = 1:L_n 
                for t = 1:L_mu
                    
                    err_A = abs(A_s - amps_FDM(s,t)) / ((A_s + amps_FDM(s,t))/2); 
                    err_phi = abs(P - phis_FDM(s,t)) / ((abs(P) + abs(phis_FDM(s,t)))/2); 

%                     err_A = (A_s - amps_FDM(s,t)) / A_s; 
%                     err_phi = (P - phis_FDM(s,t)) / P;
                    
                    err_param(s,t) = abs(err_A) + abs(err_phi);
                    
                end
            end
            
            % find minimum of error matrix
            [min_cols, inxs] = min(err_param); % minimum values of each column
            [min_g, inx_col] = min(min_cols);  % global minimum by minimizing columns
            inx_row = inxs(inx_col); % row index of global minimum

            inx_g = [inx_row inx_col]; % location of global minimum (best-fit)

            % assign best fit values
            n(k,l) = n_param(inx_g(1));
            mu(k,l) = mu_param(inx_g(2));
            err(k,l) = err_param(inx_g(1), inx_g(2));
            
        else
            continue 
        end  
    end  
end

sigma = (1.602*10^-19).*n.*mu; % S/cm, conductivity 

%% plot
% images
if data_mode == 2
    
    % plot parameters
    xpos = 0; 
    ypos = 0; 
    x_im = x_dim*2;%4; 
    y_im = y_dim*2;%4; 
    font_sz = 24; 
    
    % input amplitude
    figure(1)
    imagesc(amp.*(1./A_i));
    colorbar; colormap('parula');
    caxis([0.2 2.5]);
    title('amplitude');
    set(gca, 'fontsize', font_sz);
    set(gcf, 'Position', [xpos ypos x_im y_im]);
    saveas(gcf, strcat(savepath, 'amp.fig'));
    saveas(gcf, strcat(savepath, 'amp.png')); 
    
    % input phase
    figure(2)
    imagesc(phi);
    colorbar; colormap('jet');
    title('phase');
    caxis([0 90]); 
    set(gca, 'fontsize', font_sz);
    set(gcf, 'Position', [xpos ypos x_im y_im]);    
    saveas(gcf, strcat(savepath, 'phi.fig'));
    saveas(gcf, strcat(savepath, 'phi.png')); 
    
    % fitted density
    figure(3) 
    imagesc(n); 
    colorbar; colormap('hot'); set(gca, 'colorscale', 'log'); 
    caxis([10^19 10^20]); 
    title('density');
    set(gca, 'fontsize', font_sz); 
    set(gcf, 'Position', [xpos ypos x_im y_im]);
    saveas(gcf, strcat(savepath, 'n.fig'));
    saveas(gcf, strcat(savepath, 'n.png')); 
    
    % fitted mobility
    figure(4)
    imagesc(mu);
    colorbar; colormap('copper');
    caxis([0 200]); 
    title('mobility');
    set(gca, 'fontsize', font_sz);
    set(gcf, 'Position', [xpos ypos x_im y_im]);
    saveas(gcf, strcat(savepath, 'mu.fig'));
    saveas(gcf, strcat(savepath, 'mu.png')); 

    % fitted conductivity
    figure(5)
    imagesc(sigma);
    colorbar; colormap('parula');
    title('conductivity \rm (S/cm)');
    caxis([0 1000]); 
    set(gca, 'fontsize', font_sz);
    set(gcf, 'Position', [xpos ypos x_im y_im]);
    saveas(gcf, strcat(savepath, 'cond.fig'));
    saveas(gcf, strcat(savepath, 'cond.png')); 
    
    % fitting error
    figure(6)
    imagesc(abs(err).*100);
    colorbar; colormap('bone');
    caxis([0 100]);
    title('% error');
    set(gca, 'fontsize', font_sz);
    set(gcf, 'Position', [xpos ypos x_im y_im]);
    saveas(gcf, strcat(savepath, 'err.fig'));
    saveas(gcf, strcat(savepath, 'err.png')); 
end    

% linecuts
if data_mode == 1
    figure(1)
    yyaxis left;
    plot(x, amp.*(1/A_i), '-r', 'LineWidth', 2); 
    set(gca, 'ycolor', 'r'); 
    ylabel('\bf amplitude');
    yyaxis right; 
    plot(x, phi, '-b', 'LineWidth', 2);
    set(gca, 'ycolor', 'b');
    ylabel('\bf phase'); 
    xlabel('\bf distance \rm (\mum)'); 
    set(gca, 'XLim', [0, x(end)]); 
    title('\bf amplitude & phase'); 
    set(gca, 'fontsize', 20); 
    saveas(gcf, strcat(savepath, 'inputs.fig'));
    saveas(gcf, strcat(savepath, 'inputs.png')); 
    
    figure(2)
    yyaxis left;
    plot(x, n, 'LineWidth', 2);
    ylabel('\bf density \rm (cm^-^3)')
    set(gca, 'YLim', [10^17 1.5*10^20]);
%     set(gca, 'YScale', 'log');
    yyaxis right; 
    plot(x, mu, 'LineWidth', 2);
    ylabel('\bf mobility \rm (cm^2/V*s)'); 
    set(gca, 'XLim', [0 x(end)]); 
    set(gca, 'YLim', [mu_param(1) mu_param(end)]); 
    xlabel('\bf distance \rm (\mum)');
    title('\bf density & mobility'); 
    set(gca, 'fontsize', 20); 
    saveas(gcf, strcat(savepath, 'outputs.fig'));
    saveas(gcf, strcat(savepath, 'outputs.png')); 
    
    figure(3)
    plot(x, abs(err).*100, '-k', 'LineWidth', 2);
    ylabel('\bf % error');
    set(gca, 'YLim', [0 25]); 
    set(gca, 'XLim', [0 x(end)]);
    title('\bf error'); 
    set(gca, 'fontsize', 20); 
    saveas(gcf, strcat(savepath, 'error.fig'));
    saveas(gcf, strcat(savepath, 'error.png')); 

    figure(4)
    plot(x, sigma, 'LineWidth', 2);
    ylabel('\bf conductivity \rm (S/cm)'); 
    set(gca, 'YLim', [0 1500]); 
    set(gca, 'XLim', [0 x(end)]);
    title('\bf conductivity'); 
    set(gca, 'fontsize', 20); 
    saveas(gcf, strcat(savepath, 'cond.fig'));
    saveas(gcf, strcat(savepath, 'cond.png'));
    
end
delete(f)

%% output .gwy results file
if gwy == 1 && data_mode == 2
    output.ScanName = ScanName;
    output.yPixelSize = y_dim/y_range; 
    output.yRes = y_dim; 
    output.yRange = y_range;

    output.xPixelSize = x_dim/x_range;
    output.xRes = x_dim; 
    output.xRange = x_range; 

    output.Data = cat(3, amp.*(1/A_i), phi, mask, log10(n+1), mu, sigma .* 100, abs(err)); 
    output.Labels = {'Input Amplitude', 'Input Phase', 'Input Mask', 'Output Density', 'Output Mobility', 'Conductivity', 'Output Error'}; 
    output.Units = {'', 'deg', '', '', 'cm^2/(V*s)', 'S/m', '%'}; 

    output.FileName = strcat(savepath, output.ScanName); 

    saveasgwy(output.FileName, output.Data, output.xRes, output.yRes, 0, output.xRange, 0, output.yRange, output.Labels, output.Units, now); 
end

%% write parameters used to text file
fid = fopen(strcat(savepath, '/','parameters.txt'), 'wt'); 

fprintf(fid, 'Modeling Parameters \n'); 
fprintf(fid, 'L_n = %u\n', L_n); 
fprintf(fid, 'L_mu = %u\n', L_mu); 
fprintf(fid, 'n tested from %e to %e\n', n_param(1), n_param(end)); 
fprintf(fid, 'mu tested from %u to %u\n', mu_param(1), mu_param(end)); 
fprintf(fid, 'FDM = %u\n', FDM_mode); 
fprintf(fid, 'Si Norm = %u\n', Sinorm); 
fprintf(fid, 'Si Type = %u\n\n', type); 
fprintf(fid, 'Experimental Parameters \n'); 
fprintf(fid, 'harmonic = %u\n', harm_ex); 
fprintf(fid, 'wavenumber = %f\n', wave); 
fprintf(fid, 'SiO2 thickness = %f\n', d); 
fprintf(fid, 'incdent angle = %f\n', angle); 
fprintf(fid, 'x-range (gwy) = %f\n', x_range);
fprintf(fid, 'y-range (gwy) = %f\n\n', y_range); 
fprintf(fid, 'FFT points = %u\n\n', N_FFT);
fprintf(fid,  'Tip Parameters \n'); 
fprintf(fid, 'tip radius = %f\n', r*(10^9)); 
fprintf(fid, 'tip length = %f\n', L*(10^9)); 
fprintf(fid, 'tapping frequency = %u\n', freq_tap); 
fprintf(fid, 'tapping amplitude = %f\n', A_tap*(10^9)); 

fclose(fid); 
%%
toc