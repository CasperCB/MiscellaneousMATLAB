function Harm = demod_Sim(N, T, freq_tap, alpha, k)
% Simulate demodulation through a discrete fourier transform
% N = number of points
% T = sampling time in seconds
% freq_tap = tapping frequency in Hz
% alpha = time-dependent complex permittivity
% k = no. of harmonics desired

FFT = fft(alpha); % take fast fourier transform
P = FFT/(N/2);
alpha_FFT = P(1:N/2); % fft is symmetric, discard negative half

freq = (0:N/2-1)/T; % frequency vector in Hz
freq_norm = freq./freq_tap; % normalize frequency vector

Harm = zeros(k, 1); % initialize matrix

    for j = 1:k
        harmInx = find((abs(freq_norm-j)) < (10^-4)); % locate harmonic in frequency vector within tolerance
        Harm(j) = alpha_FFT(harmInx); % store harmonic value
    end
Harm = Harm.';
end

