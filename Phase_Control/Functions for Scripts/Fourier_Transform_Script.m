y = fft(PWC_FC);
y_shift = fftshift(fft(PWC_FC));
N = length(PWC_FC);
delta_t = 2*pi/N; % here T = 2*pi/w which implicitly states dt is in units of w
ws = 1/delta_t;
dw = 1/(N*delta_t);
w = (-(N-1)/2:(N-1)/2)*dw;
figure;stem(w,abs(y_shift));
Y = ifft(y);
tt = (0:(N-1))*delta_t;
figure;plot(tt,Y,'r');
hold on;plot(tt,PWC_FC,'b')

wc = 3;dw = 0.01;
F = [wc,wc+dw]; % frequency intervals for lowpass filter. Starts at 0 
                % frequency and cuts off at wc where the transition
                % window is wc+dw wide.
A = [1,0]; %'1' passband, '0' is stopband
dev = [10^(0.01),0.0001]; % specs for attenuation and passband ripple
[M,Wn,beta,typ] = kaiserord(F,A,dev,ws);
b = fir1(M,Wn,typ,kaiser(M+1,beta),'noscale');
out = filter(b,1,PWC_FC);
figure;plot(tt,out)