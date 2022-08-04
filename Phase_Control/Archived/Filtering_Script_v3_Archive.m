% The purpose for this function is to vary the dimension using optimization
% control function and find the dependence of the std on this parameter.

% This script is labeled archive because the following script was known to
% work and so has been archived while a newer version is modified.

clear

max_spin = 0.5;
spin = 5;%0.5:0.5:max_spin;
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
           Control_Optimization(spin(jj),1,1,1,10,0);
    %[FC_phase,T] = piecewise_repmat_leftmost_point(phiF,dt);
    %[psiG,TG] = piecewise_repmat_leftmost_point(phiG,dtG);
      
end
%[freq_pos1,freq_pos2,TSS1,TSS2] = TSS_to_SSS(freq_FC,FC_power_spectrum,...
    %freq_PWC,PWC_power_spectrum);
%figure;plot(freq_pos1,TSS1,freq_pos2,TSS2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Signal: PWC_FC
%%% Signal Time: T_PWC_FC
%%% Sampling Rate: dt_PWC_FC 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Constant for Filtering FC phase
dt_PWC_FC  = median(diff(T_PWC_FC)); %sampling interval fourier constrained phase
fs = 1/dt_PWC_FC; %sampling frequency
w0 = 1;
PB = 60*(w0/(2*pi))/(fs/2); % normalized Passband
SB = PB + (3)*(w0/(2*pi))/(fs/2);
t = 0:dt_PWC_FC:2*pi;

% Filter FC Phase
tic 
psif = Weak_Dressed_H(PWC_FC,initial_state,target_state,spin,t_final,1,1);

freq_FC = (0:(length(PWC_FC)-1))*df;
dw = 2*pi*df;
ang_freq_FC = (0:(length(PWC_FC)-1))*dw;
ang_freq_norm_FC = freq_FC/(2*pi*fs/2); % normalized frequencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant for Filtering FC phase
dt_PWC_PWC  = median(diff(T_PWC_PWC)); %sampling interval fourier constrained phase
fs_PWC = 1/dt_PWC_PWC; %sampling frequency
w0 = 1;
PB = 60*(w0/(2*pi))/(fs_PWC/2); % normalized Passband
SB = PB + (3)*(w0/(2*pi))/(fs_PWC/2);
t = 0:dt_PWC_PWC:2*pi;

% Filter PWC Phase
tic 
psif_PWC = Weak_Dressed_H(PWC_PWC,initial_state,target_state,spin,t_final,1,1);

df_PWC = fs_PWC/length(PWC_PWC);
freq_PWC = (0:(length(PWC_PWC)-1))*df_PWC;
dw = 2*pi*df_PWC;
ang_freq_PWC = (0:(length(PWC_PWC)-1))*dw;
ang_freq_norm_PWC = freq_PWC/(2*pi*fs_PWC/2); % normalized frequencies
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement designfilt through phase_filter function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[FC_filtered,A,D] = phase_filter(PWC_FC,PB(kk),SB(kk));
FC_filtered = FC_filtered(D+1:end);

[PWC_filtered,A,D] = phase_filter(PWC_PWC,PB(kk),SB(kk));
PWC_filtered = PWC_filtered(D+1:end);


figure;stem(ang_freq_FC,abs(fft(PWC_FC))/max(abs(fft(PWC_FC))));
hold on;stem(ang_freq_FC,abs(fft(FC_filtered))/max(abs(fft(FC_filtered))));
[H,W] = freqz(A,floor(length(PWC_FC)/2)); % W is gvien in terms of rad/s
hold on; plot((W/pi)*(2*pi*fs/2),abs(H),'r');% W has been normalized (it ranges between 0 and pi)
ylabel('Normalized Power Spectrum')
xlabel('\omega_0')
axis([0,(2*pi*fs/2),0,1])
figure;plot(t(1:end-1),FC_filtered,'k',t(1:end-1),PWC_FC,'m')

d_FC = abs(psi_f_filt'*psif);

toc

figure;plot(PB*(2*pi*fs/2),d_FC,'om');title('Fourier Constrained Phase')
xlabel('lowpass Frequency [\Omega / 2\pi]');ylabel('abs(Fidelity - Filtered Phase Fidelity)')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement Butter Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PB_butter = (1:1:50)*(w0/(2*pi))/(fs/2); % normalized Passband
d_FC_butter = zeros(length(PB_butter),1);
psif = Weak_Dressed_H(PWC_FC,initial_state,target_state,spin,t_final,1,1);
Fid_FC = d_FC_butter;

PB_butter_PWC = (1:1:120)*(w0/(2*pi))/(fs_PWC/2); % normalized Passband
d_PWC_butter = zeros(length(PB_butter_PWC),1);
psif_PWC = Weak_Dressed_H(PWC_PWC,initial_state,target_state,spin,t_final,1,1);
Fid_PWC = d_PWC_butter;

for kk = 1:length(PB_butter)
    [B,A] = butter(5,PB_butter(kk)); % by default butter produces a lowpass filter
%     [H_butter,W_butter] = freqz(B,A,floor(length(PWC_FC)/2));
% 
%     figure;plot((W_butter/pi)*(2*pi*fs/2),abs(H_butter),'r')
%     hold on;stem(ang_freq_FC,(abs(fft(PWC_FC))/max(abs(fft(PWC_FC)))).^2)
%     xlabel('\omega_0');ylabel('Power Spectrum')
%     legend('Magnitude Response Function','Power Spectrum Squared')

    % filter signal using butter with delay introduced
%     FC_butter = filter(B,A,PWC_FC);
%     figure;plot(T_PWC_FC,FC_butter)
%     hold on;plot(T_PWC_FC,PWC_FC)
%     legend('Filtered','Not-Filtered')

    % filter signal using butter conpensating for delay
     FC_butter = filtfilt(B,A,PWC_FC);
%      figure;plot(T_PWC_FC,FC_butter)
%      hold on;plot(T_PWC_FC,PWC_FC)
%      legend('Filtered','Not-Filtered')
%      xlabel('1/\Omega');ylabel('Phase Amplitude')
%      figure;hold on;stem(ang_freq_FC,(abs(fft(FC_butter))/max(abs(fft(FC_butter)))).^2)

    % Calculate Fidelity
    psif_filt = Weak_Dressed_H(FC_butter,initial_state,target_state,spin,t_final,1,1);
    d_FC_butter(kk) = abs(psif_filt'*psif)^2;
%     Fid_FC(kk) = abs(psif_filt'*target_state)^2;
end
clear kk 
%%%%%%%%%%%%%%%%%%%%%%PWC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk = 1:length(PB_butter_PWC)
    [B,A] = butter(5,PB_butter_PWC(kk)); % by default butter produces a lowpass filter
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
    d_PWC_butter(kk) = abs(psif_filt'*psif_PWC)^2;
    Fid_PWC(kk) = abs(psif_filt'*target_state)^2;
end


figure;plot(PB_butter_PWC*(2*pi*fs_PWC/2),d_PWC_butter,'*r')
xlabel('Passband [\omega_0]')
ylabel('|\langle \Psi_{evolved-filt}|\Psi_{evolved} \rangle|^2')
str = sprintf('Dimension %i and spin %i',2*spin+1,spin);
title(str)
%figure;plot(PB_butter_PWC*(2*pi*fs/2),Fid_PWC,'ob')
grid on
hold on;plot(PB_butter*(2*pi*fs/2),d_FC_butter,'*b');
% figure;plot(PB_butter*(2*pi*fs/2),Fid_FC,'ob')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement Bessel Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[b,a] = besself(10,PB*(2*pi*fs/2));
freqs(b,a)
[H_bessel,W_bessel] = freqz(b,a,floor(length(PWC_FC)/2));
figure;plot((W_bessel/pi)*(2*pi*fs/2),abs(H_bessel),'r')
hold on;stem(ang_freq_FC,abs(fft(PWC_FC))/max(abs(fft(PWC_FC))))

FC_bessel = filter(b,a,PWC_FC);
figure;plot(FC_bessel)
