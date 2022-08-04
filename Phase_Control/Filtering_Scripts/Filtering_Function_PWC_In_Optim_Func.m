% The purpose for this function is to vary the dimension using optimization
% control function and find the dependence of the std on this parameter.

% This script is labeled archive because the following script was known to
% work and so has been archived while a newer version is modified.

function [wmax_out,d_PWC_butter,PB_PWC_out,dim] = Filtering_Function_PWC_In_Optim_Func(...
    initial_state,target_state,PWC_PWC,T_PWC_PWC,t_final,spin,PB_PWC_0,PB_PWC_f)

dim = 2*spin +1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Signal: PWC_PWC
%%% Signal Time: T_PWC_PWC
%%% Sampling Rate: dt_PWC_PWC 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w0 = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant for Filtering PWC phase
dt_PWC_PWC  = median(diff(T_PWC_PWC)); %sampling interval fourier constrained phase
fs = 1/dt_PWC_PWC; %sampling frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement Butter Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PB_butter_PWC = (PB_PWC_0:PB_PWC_f)*(w0/(2*pi))/(fs/2); % normalized Passband with angular sample frequency conversion
% PB_butter is a vector of max low pass frequencies
psif_PWC = Weak_Dressed_H(PWC_PWC,initial_state,target_state,spin,t_final,1,1);
d_PWC_butter = zeros(length(PB_butter_PWC),1);

for jj = 1:length(PB_butter_PWC)
    [B,A] = butter(5,PB_butter_PWC(jj)); % by default butter produces a lowpass filter
    PWC_butter = filtfilt(B,A,PWC_PWC);
    psif_filt= Weak_Dressed_H(PWC_butter,initial_state,...
        target_state,spin,t_final,1,1);
    d_PWC_butter(jj) = abs(psif_filt'*psif_PWC)^2;
   
if d_PWC_butter(jj) >= 0.99
    wmax_PWC = PB_butter_PWC(jj);
    break
end
end

PB_PWC_out = PB_butter_PWC*(2*pi*fs/2);
wmax_out = wmax_PWC*(2*pi*fs/2);% unnormalize and convert to angular frequency
end