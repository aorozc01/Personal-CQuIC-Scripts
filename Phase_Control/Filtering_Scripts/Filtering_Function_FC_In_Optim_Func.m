% The purpose for this function is to vary the dimension using optimization
% control function and find the dependence of the std on this parameter.

% This script is labeled archive because the following script was known to
% work and so has been archived while a newer version is modified.

function [wmax_out,d_FC_butter,PB_FC_out,dim] = Filtering_Function_FC_In_Optim_Func(...
    initial_state,target_state,PWC_FC,T_PWC_FC,t_final,spin,PB_FC_0,PB_FC_f)

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

PB_butter_FC = (PB_FC_0:0.5:PB_FC_f)*(w0/(2*pi))/(fs/2); % normalized Passband with angular sample frequency conversion
psif = Weak_Dressed_H(PWC_FC,initial_state,target_state,spin,t_final,1,1);


for kk = 1:length(PB_butter_FC)
    [B,A] = butter(5,PB_butter_FC(kk)); % by default butter produces a lowpass filter
    FC_butter = filtfilt(B,A,PWC_FC);
    psif_filt = Weak_Dressed_H(FC_butter,initial_state,target_state,spin,t_final,1,1);
    d_FC_butter = abs(psif_filt'*psif)^2;
        if d_FC_butter >= 0.99
            wmax_FC = PB_butter_FC(kk);
        break
        end
end
PB_FC_out = PB_butter_FC*(2*pi*fs/2);
wmax_out = wmax_FC*(2*pi*fs/2);% unnormalize and convert to angular frequency
end