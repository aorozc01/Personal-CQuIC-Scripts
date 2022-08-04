% The purpose for this function is to vary the dimension using optimization
% control function and find the dependence of the std on this parameter.

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
PB = (5:10)*(w0/(2*pi))/(fs/2); % normalized Passband
SB = PB + (1/2)*(w0/(2*pi))/(fs/2);
t = 0:dt_PWC_FC:2*pi;
d_FC = zeros(length(PB),1);
% Filter FC Phase
tic 
Fid = -Weak_Dressed_H(PWC_FC,initial_state,target_state,spin,t_final,1,1);

%this part is to analyse a sampled phase
%{
wmax = 20;
dt_samp = 2*pi/100;%(4*pi)/wmax;
[t_samples,samples] = samplegen(PWC_FC,T_PWC_FC,dt_samp,2,FC_dt);
dt_sample = median(diff(t_samples));
fs_samp = 1/dt_sample;
w0 = 1;
PB_samp = 1*(w0/(2*pi))/(fs_samp/2); % normalized Passband
SB_samp = PB_samp + (1/2)*(w0/(2*pi))/(fs_samp/2);
df_samp = 1/(length(samples)*dt_samp);
figure;plot(T_PWC_FC,PWC_FC,'or',t_samples,samples,'*b')
w_samp = (0:(length(samples)-1))*2*pi*df_samp;
figure;stem(w_samp,abs(fft(samples))/max(abs(fft(samples))))
[FC_samp_filtered,A_samp,D_samp] = phase_filter(samples,PB_samp,SB_samp);
figure;plot(FC_samp_filtered)
%}
for kk = 1:length(PB)
    [FC_filtered,A,D] = phase_filter(PWC_FC,PB(kk),SB(kk));
    FC_filtered = FC_filtered(D+1:end);
    
   %[b,a] = besself(5,PB);

%     [B,A] = butter(20,PB(kk),'low');
%     FC_filtered = filter(B,A,PWC_FC);
    
    %{
    freq_FC = (0:(length(PWC_FC)-1))*df;
    dw = 2*pi*df;
    ang_freq_FC = (0:(length(PWC_FC)-1))*dw;
    freq_norm_FC = freq_FC/(fs/2); % normalized frequencies
    figure;stem(ang_freq_FC,abs(fft(PWC_FC))/max(abs(fft(PWC_FC))));
    hold on;stem(freq_norm_FC*(fs/2)*(w0/2*pi),abs(fft(FC_filtered))/max(abs(fft(FC_filtered))));
    [H,W] = freqz(A,floor(length(PWC_FC)/2));
    hold on; plot(W/pi,abs(H),'r'); % W has been normalized (it ranges between 0 and pi)
    figure;plot(t(1:end-1),FC_filtered,'k',t(1:end-1),PWC_FC,'m')
    %}
    psi_f_filt = -Weak_Dressed_H(FC_filtered,initial_state,target_state,spin,t_final,1,1);
    d_FC(kk) = abs(psi_f_filt'*psi_f(:,:,end));

end

% Filter PWC Phase
%{
dt_PWC  = median(diff(T_PWC_PWC));
fs_PWC = 1/dt_PWC;
PB_PWC = (0.1:0.1:dim)/(fs_PWC/2); % normalized Passband
SB_PWC = PB_PWC + 0.25;
PWC_Fid = -PWC_Phase(PWC_PWC,initial_state,target_state,spin,t_final,1,1); 
d_PWC = zeros(length(PB),1);
%t = (0:length(PWC_PWC)-1)*dt_PWC;

clear ll
for ll=1:length(PB_PWC)
    [PWC_filtered,A,D] = phase_filter(PWC_PWC,PB_PWC(ll),SB_PWC(ll));
    PWC_filtered = PWC_filtered(D+1:end);
    %{
    freq_PWC = (0:(length(PWC_PWC)-1))*df;
    freq_norm_PWC = freq_PWC/(fs_PWC/2); % normalized frequencies
    figure;stem(freq_norm_PWC,abs(fft(PWC_PWC)));
    hold on;stem(freq_norm_PWC,abs(fft(PWC_filtered)));
    figure;plot(t,PWC_filtered,'k',t,PWC_PWC,'m')
    %}
    [~,~,~,~,psi_f_filt_PWC] = PWC_Phase(PWC_filtered,initial_state,...
        target_state,spin,t_final,1,1); 
    d_PWC(ll) = abs(psi_f_filt_PWC(:,:,end)'*psi_f(:,:,end));
end
%}
toc
figure;plot(PB*(2*pi*fs/2),d_FC,'om');title('Fourier Constrained Phase')
xlabel('lowpass Frequency [\Omega / 2\pi]');ylabel('abs(Fidelity - Filtered Phase Fidelity)')
%figure;plot(PB_PWC*(2*pi*fs_PWC/2),d_PWC,'ok');title('Unconstrained Phase')
%xlabel('lowpass Frequency [\Omega / 2\pi]');ylabel('abs(Fidelity - Filtered Phase Fidelity)')

%{
figure;grid on;
plot(dim,std_FC,'r',dim,std_PWC,'b')
hold on;plot(dim,std_FC,'or',dim,std_PWC,'ob')

title('Power Spectrum Standard Deviation')
xlabel('Hilbert Space Dimension')
ylabel('Standard Deviation')
legend('FC','PWC','location','NorthWest')
%}