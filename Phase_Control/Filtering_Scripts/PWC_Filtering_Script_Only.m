% The purpose for this function is to vary the dimension using optimization
% control function and find the dependence of the std on this parameter.

% This script is labeled archive because the following script was known to
% work and so has been archived while a newer version is modified.

clear

max_spin = 0.5;
spin = 7;%0.5:0.5:max_spin;
dim = 2*spin + 1;
K = length(spin);
plots = 1;% if plots = 1 then spectrum plots occur if anything else no.
% scalars 
std_FC = zeros(K,1);
std_PWC = zeros(K,1);
Avg_Fid_FC_per_spin = zeros(K,1);
Avg_Fid_PWC_per_spin = zeros(K,1);

for jj = 1:K;
    [std_FC(jj),std_PWC(jj),Avg_Fid_FC_per_spin(jj),Avg_Fid_PWC_per_spin(jj),... 
          freq_FC,FC_power_spectrum,freq_PWC,PWC_power_spectrum,PWC_FC,T_PWC_FC,...
          PWC_PWC,T_PWC_PWC,FC_dt,PWC_dt,FC_FT,PWC_FT,FC_phase,PWC_phase,aopt,...
          PWC_phiF_fft,PWC_phiG_fft,df,PWC_aopt,target_state,t_final,steps_FC,...
          initial_state,psi_f] =...
           Control_Optimization(spin(jj),1,1,1,10,1);
    %[FC_phase,T] = piecewise_repmat_leftmost_point(phiF,dt);
    %[psiG,TG] = piecewise_repmat_leftmost_point(phiG,dtG);
      
end
%[freq_pos1,freq_pos2,TSS1,TSS2] = TSS_to_SSS(freq_FC,FC_power_spectrum,...
    %freq_PWC,PWC_power_spectrum);
%figure;plot(freq_pos1,TSS1,freq_pos2,TSS2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Signal: PWC_PWC
%%% Signal Time: T_PWC_PWC
%%% Sampling Rate: dt_PWC_PWC 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constant for Filtering FC phase
dt_PWC_PWC  = median(diff(T_PWC_PWC)); %sampling interval fourier constrained phase
fs = 1/dt_PWC_PWC; %sampling frequency
w0 = 1;
PB = 60*(w0/(2*pi))/(fs/2); % normalized Passband
SB = PB + (3)*(w0/(2*pi))/(fs/2);
t = 0:dt_PWC_PWC:2*pi;

% Filter PWC Phase
tic 
psif = -Weak_Dressed_H(PWC_PWC,initial_state,target_state,spin,t_final,1,1);

df_PWC = fs/length(PWC_PWC);
freq_PWC = (0:(length(PWC_PWC)-1))*df_PWC;
dw = 2*pi*df_PWC;
ang_freq_PWC = (0:(length(PWC_PWC)-1))*dw;
ang_freq_norm_PWC = freq_PWC/(2*pi*fs/2); % normalized frequencies
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement designfilt through phase_filter function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[PWC_filtered,A,D] = phase_filter(PWC_PWC,PB(kk),SB(kk));
PWC_filtered = PWC_filtered(D+1:end);



figure;stem(ang_freq_PWC,abs(fft(PWC_PWC))/max(abs(fft(PWC_PWC))));
hold on;stem(ang_freq_PWC,abs(fft(PWC_filtered))/max(abs(fft(PWC_filtered))));
[H,W] = freqz(A,floor(length(PWC_PWC)/2)); % W is gvien in terms of rad/s
hold on; plot((W/pi)*(2*pi*fs/2),abs(H),'r');% W has been normalized (it ranges between 0 and pi)
ylabel('Normalized Power Spectrum')
xlabel('\omega_0')
axis([0,(2*pi*fs/2),0,1])
figure;plot(t(1:end-1),PWC_filtered,'k',t(1:end-1),PWC_PWC,'m')

d_PWC = abs(psi_f_filt'*psif);

toc

figure;plot(PB*(2*pi*fs/2),d_PWC,'om');title('Fourier Constrained Phase')
xlabel('lowpass Frequency [\Omega / 2\pi]');ylabel('abs(Fidelity - Filtered Phase Fidelity)')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement Butter Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PB_butter = (1:1:200)*(w0/(2*pi))/(fs/2); % normalized Passband
d_PWC_butter = zeros(length(PB_butter),1);
psif = Weak_Dressed_H(PWC_PWC,initial_state,target_state,spin,t_final,1,1);
Fid_PWC = d_PWC_butter;
for kk = 1:length(PB_butter)
    [B,A] = butter(5,PB_butter(kk)); % by default butter produces a lowpass filter
%   [H_butter,W_butter] = freqz(B,A,floor(length(PWC_PWC)/2));
% 
%     figure;plot((W_butter/pi)*(2*pi*fs/2),abs(H_butter),'r')
%     hold on;plot(ang_freq_PWC,abs(fftshift(fft(PWC_PWC)))/max(abs(fft(PWC_PWC))))
%     xlabel('\omega_0');ylabel('Power Spectrum Squared')
%     legend('Magnitude Response Function','Power Spectrum Squared')

    % filter signal using butter with delay introduced
%     PWC_butter = filter(B,A,PWC_PWC);
%     figure;plot(T_PWC_PWC,PWC_butter)
%     hold on;plot(T_PWC_PWC,PWC_PWC)
%     legend('Filtered','Not-Filtered')

    % filter signal using butter conpensating for delay
      PWC_butter = filtfilt(B,A,PWC_PWC);
%      figure;plot(T_PWC_PWC,PWC_butter)
%      hold on;plot(T_PWC_PWC,PWC_PWC)
%      legend('Filtered','Not-Filtered')
%      xlabel('1/\Omega');ylabel('Phase Amplitude')

    % Calculate Fidelity
    [psif_filt,Fid_butter] = Weak_Dressed_H(PWC_butter,initial_state,target_state,spin,t_final,1,1);
    d_PWC_butter(kk) = abs(psif_filt'*psif);
    Fid_PWC(kk) = abs(psif_filt'*target_state)^2;
end
figure;plot(PB_butter*(2*pi*fs/2),d_PWC_butter,'*r')
xlabel('Passband [\omega_0]')
ylabel('\langle \Psi_{evolved-filt}|\Psi_{evolved} \rangle')
str = sprintf('Dimension %i and spin %i',2*spin+1,spin);
title(str)
figure;plot(PB_butter*(2*pi*fs/2),Fid_PWC,'ob')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement Bessel Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[b,a] = besself(10,PB*(2*pi*fs/2));
freqs(b,a)
[H_bessel,W_bessel] = freqz(b,a,floor(length(PWC_PWC)/2));
figure;plot((W_bessel/pi)*(2*pi*fs/2),abs(H_bessel),'r')
hold on;stem(ang_freq_PWC,abs(fft(PWC_PWC))/max(abs(fft(PWC_PWC))))

PWC_bessel = filter(b,a,PWC_PWC);
figure;plot(PWC_bessel)
