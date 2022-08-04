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
           Control_Optimization(spin(jj),1,1,1,10,1);
    %[FC_phase,T] = piecewise_repmat_leftmost_point(phiF,dt);
    %[psiG,TG] = piecewise_repmat_leftmost_point(phiG,dtG);
      
end
%[freq_pos1,freq_pos2,TSS1,TSS2] = TSS_to_SSS(freq_FC,FC_power_spectrum,...
    %freq_PWC,PWC_power_spectrum);
%figure;plot(freq_pos1,TSS1,freq_pos2,TSS2)

% Constant for Filtering FC phase
dt_PWC_FC  = median(diff(T_PWC_FC)); %sampling interval fourier constrained phase
fs = 1/dt_PWC_FC; %sampling frequency
w0 = 1;
PB = 60*(w0/(2*pi))/(fs/2); % normalized Passband
SB = PB + (3)*(w0/(2*pi))/(fs/2);
t = 0:dt_PWC_FC:2*pi;

% Filter FC Phase
tic 
psif = -Weak_Dressed_H(PWC_FC,initial_state,target_state,spin,t_final,1,1);

[FC_filtered,A,D] = phase_filter(PWC_FC,PB(kk),SB(kk));
FC_filtered = FC_filtered(D+1:end);
    
%[b,a] = besself(5,PB);
%[B,A] = butter(20,PB(kk),'low');
%FC_filtered = filter(B,A,PWC_FC);

freq_FC = (0:(length(PWC_FC)-1))*df;
dw = 2*pi*df;
ang_freq_FC = (0:(length(PWC_FC)-1))*dw;
ang_freq_norm_FC = freq_FC/(2*pi*fs/2); % normalized frequencies

figure;stem(ang_freq_FC,abs(fft(PWC_FC))/max(abs(fft(PWC_FC))));
hold on;stem(ang_freq_FC,abs(fft(FC_filtered))/max(abs(fft(FC_filtered))));
[H,W] = freqz(A,floor(length(PWC_FC)/2)); % W is gvien in terms of rad/s
hold on; plot((W/pi)*(2*pi*fs/2),abs(H),'r'); % W has been normalized (it ranges between 0 and pi)
figure;plot(t(1:end-1),FC_filtered,'k',t(1:end-1),PWC_FC,'m')
d_FC = abs(psi_f_filt'*psif);

toc

figure;plot(PB*(2*pi*fs/2),d_FC,'om');title('Fourier Constrained Phase')
xlabel('lowpass Frequency [\Omega / 2\pi]');ylabel('abs(Fidelity - Filtered Phase Fidelity)')
