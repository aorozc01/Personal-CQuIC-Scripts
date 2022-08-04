% The purpose for this function is to vary the dimension using optimization
% control function and find the dependence of the std on this parameter.

% This script is labeled archive because the following script was known to
% work and so has been archived while a newer version is modified.

function [d_FC_butter,PB_FC_out,dim] = Filtering_Function_FC_Only(spin,PB_FC_0,PB_FC_f)

dim = 2*spin +1;       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Signal: PWC_FC
%%% Signal Time: T_PWC_FC
%%% Sampling Rate: dt_PWC_FC 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w0 = 1;
% Constant for Filtering FC phase
dt_PWC_FC  = median(diff(T_PWC_FC)); %sampling interval fourier constrained phase
fs = 1/dt_PWC_FC; %sampling frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement Butter Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PB_butter_FC = (PB_FC_0:1:PB_FC_f)*(w0/(2*pi))/(fs/2); % normalized Passband
d_FC_butter = zeros(length(PB_butter_FC),1);
psif = Weak_Dressed_H(PWC_FC,initial_state,target_state,spin,t_final,1,1);
Fid_FC = d_FC_butter;
Com_Fid_FC = Fid_FC;

for kk = 1:length(PB_butter_FC)
    [B,A] = butter(5,PB_butter_FC(kk)); % by default butter produces a lowpass filter
    FC_butter = filtfilt(B,A,PWC_FC);
    [psif_filt,Fid_FC(kk)] = Weak_Dressed_H(FC_butter,initial_state,target_state,spin,t_final,1,1);
    d_FC_butter(kk) = abs(psif_filt'*psif)^2;
    Com_Fid_FC(kk) = abs(psif_filt'*target_state)^2;
end
PB_FC_out = PB_butter_FC*(2*pi*fs/2);
end