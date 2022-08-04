% The purpose for this function is to vary the dimension using optimization
% control function and find the dependence of the std on this parameter.

% This script is labeled archive because the following script was known to
% work and so has been archived while a newer version is modified.

function [PB_PWC_out,d_PWC_butter,dim] = Filtering_Function_PWC_Only(spin,PB_PWC)

dim = 2*spin +1;

[~,PWC_PWC,T_PWC_PWC,~,~,~,target_state,t_final,initial_state,~] =...
           Control_Optimization_PWC_Only(spin,1,1,10,1);
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Signal: PWC_FC
%%% Signal Time: T_PWC_FC
%%% Sampling Rate: dt_PWC_FC 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w0 = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant for Filtering PWC phase
dt_PWC_PWC  = median(diff(T_PWC_PWC)); %sampling interval fourier constrained phase
fs_PWC = 1/dt_PWC_PWC; %sampling frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement Butter Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PB_butter_PWC = (1:1:PB_PWC)*(w0/(2*pi))/(fs_PWC/2); % normalized Passband
d_PWC_butter = zeros(length(PB_butter_PWC),1);
psif_PWC = Weak_Dressed_H(PWC_PWC,initial_state,target_state,spin,t_final,1,1);
Fid_PWC = d_PWC_butter;
Com_Fid_PWC = Fid_PWC;

%%%%%%%%%%%%%%%%%%%%%%PWC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk = 1:length(PB_butter_PWC)
    [B,A] = butter(5,PB_butter_PWC(kk)); % by default butter produces a lowpass filter
    PWC_butter = filtfilt(B,A,PWC_PWC);
    [psif_filt,Fid_PWC(kk)] = Weak_Dressed_H(PWC_butter,initial_state,target_state,spin,t_final,1,1);
    d_PWC_butter(kk) = abs(psif_filt'*psif_PWC)^2;
    Com_Fid_PWC(kk) = abs(psif_filt'*target_state)^2;
end
PB_PWC_out = PB_butter_PWC*(2*pi*fs_PWC/2);
%{
figure;plot(PB_butter_PWC*(2*pi*fs_PWC/2),d_PWC_butter,'*r')
xlabel('Passband [\omega_0]')
ylabel('|\langle \Psi_{evolved-filt}|\Psi_{evolved} \rangle|^2')
str = sprintf('Dimension %i and spin %i',2*spin+1,spin);
title(str)
%figure;plot(PB_butter_PWC*(2*pi*fs/2),Fid_PWC,'ob')
grid on
hold on;plot(PB_butter*(2*pi*fs/2),d_FC_butter,'*b');
% figure;plot(PB_butter*(2*pi*fs/2),Fid_FC,'ob')
%}