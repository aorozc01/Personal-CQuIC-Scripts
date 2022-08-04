function [X,f] = fftplot(data,dt)
X = fft(data)/length(data); % fourier transform of data
X_mag = abs(X); % magnitude for FT
X_shift_mag = abs(fftshift(X))/max(X_mag);
X_shifted_mag_sqrd = X_shift_mag.^2;
L = length(X); % number of bins
f = (0:(L-1))*(1/(L*dt));
end
