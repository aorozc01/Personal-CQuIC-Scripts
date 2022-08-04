function [freq,power_spectra]=PowerSpectra_PWC_DFT(f,dt)
% PowerSpectra_PWC_DFT(f,dt) calculates the power spectra (modulus squared 
% of the discrete Fourier transform) of the components of the column vector
% f. Note: This function is a modified version of the FTSetOfPoints function.
% UNDER CONSTRUCTION

freq = (-50:dt:50)'; % call this length p
n = (1:length(f))';   % call this length m
[Freq,N] = meshgrid(freq,n); % Freq has m replications of freq as columns 
                             % and N has p replications of n as rows 
basis_matrix = exp(-1i*Freq.*N*5); % this is a mxp matrix
F = f'*basis_matrix; % this is a row vector with length p
modulus_sqrd = (1/2*pi)*abs(F'.*(sin(freq*dt/2)./(freq/2))).^2;
power_spectra = modulus_sqrd/max(modulus_sqrd);
end