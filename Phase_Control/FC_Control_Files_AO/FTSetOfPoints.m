function [w,mod_sqrd_FT,FT]=FTSetOfPoints(f,dt)
% This function calculates the discrete Fourier transform of the components
% of the column vector f with sample point spacing dt.

w = linspace(-500,500,length(f));
F = zeros(length(f),1);
FT = zeros(length(f), 1);
for jj = 1:length(w); % evaluate the transform at each value of w (the frequency domain)
    F(jj) = 0;
for kk = 1:length(f) % sum through the sample points (i.e., points of f)
% F(jj) =(abs(f(kk)*exp(-1i*(kk-1)*w(jj)*dt))*dt*(sin (w(jj)*dt/2)/(w(jj)*dt/2)))^2 + F(jj);% Deutsch calculation
F(jj) = f(kk)*exp(-1i*(kk-1)*w(jj)*dt)+ F(jj);
end
FT(jj) = F(jj)*exp(-1i*w(jj)*dt/2)*dt*sinc(w(jj)*dt/2);
end
mod_sqrd_FT = (1/2*pi)*abs(FT).^2;
mod_sqrd_FT = mod_sqrd_FT/max(mod_sqrd_FT);
end
% figure;plot(w,H)