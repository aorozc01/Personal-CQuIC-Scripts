% The purpose for this function is to vary the dimension using optimization
% control function and find the dependence of the std on this parameter.

% This script is labeled archive because the following script was known to
% work and so has been archived while a newer version is modified.

clear

max_spin = 0.5;
spin = 1;%0.5:0.5:max_spin;
dim = 2*spin + 1;
K = length(spin);
plots = 1;% if plots = 1 then spectrum plots occur if anything else no.
% scalars
Avg_Fid_FC_per_spin = zeros(K,1);
Avg_Fid_PWC_per_spin = zeros(K,1);

for jj = 1:K;
    [~,~,Avg_Fid_FC_per_spin(jj),Avg_Fid_PWC_per_spin(jj),... 
          ~,~,~,~,PWC_FC,T_PWC_FC,...
          PWC_PWC,T_PWC_PWC,~,~,~,~,FC_phase,PWC_phase,~,...
          ~,~,df,~,target_state,t_final,steps_FC,...
          initial_state,psi_f] =...
           Control_Optimization(spin(jj),1,1,1,10,0);
    %[FC_phase,T] = piecewise_repmat_leftmost_point(phiF,dt);
    %[psiG,TG] = piecewise_repmat_leftmost_point(phiG,dtG);
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Signal: PWC_FC
%%% Signal Time: T_PWC_FC
%%% Sampling Rate: dt_PWC_FC 

% Constant for Filtering FC phase

dt_PWC_FC  = median(diff(T_PWC_FC)); %sampling interval fourier constrained phase
fs = 1/dt_PWC_FC; %sampling frequency
freq_FC = (0:(length(PWC_FC)-1))*df;
dw = 2*pi*df;
ang_freq_FC = (0:(length(PWC_FC)-1))*dw;
ang_freq_norm_FC = freq_FC/(2*pi*fs/2); % normalized frequencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant for Filtering PWC phase
dt_PWC_PWC  = median(diff(T_PWC_PWC)); %sampling interval fourier constrained phase
fs_PWC = 1/dt_PWC_PWC; %sampling frequency
w0 = 1;
PB = 60*(w0/(2*pi))/(fs_PWC/2); % normalized Passband
SB = PB + (3)*(w0/(2*pi))/(fs_PWC/2);
t = 0:dt_PWC_PWC:2*pi;

% Filter PWC Phase

df_PWC = fs_PWC/length(PWC_PWC);
freq_PWC = (0:(length(PWC_PWC)-1))*df_PWC;
dw = 2*pi*df_PWC;
ang_freq_PWC = (0:(length(PWC_PWC)-1))*dw;
ang_freq_norm_PWC = freq_PWC/(2*pi*fs_PWC/2); % normalized frequencies
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement Butter Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wmax_FC = 2*dim;
PB_butter_FC = (1:1:wmax_FC)*(w0/(2*pi))/(fs/2); % normalized Passband
d_FC_butter = zeros(length(PB_butter_FC),1);
psif = Weak_Dressed_H(PWC_FC,initial_state,target_state,spin,t_final,1,1);
Fid_FC = d_FC_butter;

wmax_PWC = 10*wmax_FC;
PB_PWC = 20*3 -30;
PB_PWC_0 = 1;%max(dim(ii),PB_PWC-80);
PB_PWC_f = PB_PWC+50;

PB_butter_PWC = (PB_PWC_0:PB_PWC_f)*(w0/(2*pi))/(fs_PWC/2);%(1:1:wmax_PWC)*(w0/(2*pi))/(fs_PWC/2); % normalized Passband
d_PWC_butter = zeros(length(PB_butter_PWC),1);
psif_PWC = Weak_Dressed_H(PWC_PWC,initial_state,target_state,spin,t_final,1,1);
Fid_PWC = d_PWC_butter;

for kk = 1:length(PB_butter_FC)
    [B,A] = butter(5,PB_butter_FC(kk)); % by default butter produces a lowpass filter
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
    
    Fid_FC(kk) = abs(psif_filt'*target_state)^2;
if d_FC_butter(kk) >= 0.99
    FC_w0 = kk;
    break
end
end

clear kk 
%%%%%%%%%%%%%%%%%%%%%%PWC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for jj = 1:length(PB_butter_PWC)
    [B,A] = butter(5,PB_butter_PWC(jj)); % by default butter produces a lowpass filter
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
    [psif_filt,Fid_butter] = Weak_Dressed_H(PWC_butter,initial_state,...
        target_state,spin,t_final,1,1);
    d_PWC_butter(jj) = abs(psif_filt'*psif_PWC)^2;
    Fid_PWC(jj) = abs(psif_filt'*target_state)^2;
if d_PWC_butter(jj) >= 0.99
    PWC_w0 = jj;
    break
end
end

%% Plot Results
figure;plot(PB_butter_PWC*(2*pi*fs_PWC/2),d_PWC_butter,'*r')
xlabel('Passband [\omega_0]')
ylabel('|\langle \Psi_{evolved-filt}|\Psi_{evolved} \rangle|^2')
str = sprintf('Dimension %i and spin %i',2*spin+1,spin);
title(str)
%figure;plot(PB_butter_PWC*(2*pi*fs/2),Fid_PWC,'ob')
grid on
hold on;plot(PB_butter_FC*(2*pi*fs/2),d_FC_butter,'*b');
% figure;plot(PB_butter*(2*pi*fs/2),Fid_FC,'ob')
