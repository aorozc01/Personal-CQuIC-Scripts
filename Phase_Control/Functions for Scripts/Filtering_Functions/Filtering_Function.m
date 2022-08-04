  % The purpose for this function is to vary the dimension using optimization
% control function and find the dependence of the std on this parameter.

% This script is labeled archive because the following script was known to
% work and so has been archived while a newer version is modified.

function [PB_FC_out,d_FC_butter,PB_PWC_out,d_PWC_butter,dim] = Filtering_Function(spin,PB_FC,PB_PWC)

dim = 2*spin +1;

plots = 0;% if plots = 1 then spectrum plots occur if anything else no.
[~,~,~,~,PWC_FC,T_PWC_FC,PWC_PWC,T_PWC_PWC,~,~,~,~,~,~,~,...
~,~,~,~,target_state,t_final,~,initial_state,~] =...
           Control_Optimization_for_filtering(spin,1,1,1,10,plots,1);
       
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
% Constant for Filtering PWC phase
dt_PWC_PWC  = median(diff(T_PWC_PWC)); %sampling interval fourier constrained phase
fs_PWC = 1/dt_PWC_PWC; %sampling frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement Butter Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the input PB_FC is angular frequency
PB_butter_FC = (1:1:PB_FC)*(w0/(2*pi))/(fs/2); % normalized Passband
d_FC_butter = zeros(length(PB_butter_FC),1);
psif = Weak_Dressed_H(PWC_FC,initial_state,target_state,spin,t_final,1,1);
Fid_FC = d_FC_butter;
Com_Fid_FC = Fid_FC;

PB_butter_PWC = (1:1:PB_PWC)*(w0/(2*pi))/(fs_PWC/2); % angular freq. normalized Passband (the wo/(2pi) is used to convert fs to ws, the angular sample frequency)
d_PWC_butter = zeros(length(PB_butter_PWC),1);
psif_PWC = Weak_Dressed_H(PWC_PWC,initial_state,target_state,spin,t_final,1,1);
Fid_PWC = d_PWC_butter;
Com_Fid_PWC = Fid_PWC;


for kk = 1:length(PB_butter_FC)
    [B,A] = butter(5,PB_butter_FC(kk)); % by default butter produces a lowpass filter
    FC_butter = filtfilt(B,A,PWC_FC);
    [psif_filt,Fid_FC(kk)] = Weak_Dressed_H(FC_butter,initial_state,target_state,spin,t_final,1,1);
    d_FC_butter(kk) = abs(psif_filt'*psif)^2;
    Com_Fid_FC(kk) = abs(psif_filt'*target_state)^2;
end
clear kk
%%%%%%%%%%%%%%%%%%%%%%PWC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk = 1:length(PB_butter_PWC)
    [B,A] = butter(5,PB_butter_PWC(kk)); % by default butter produces a lowpass filter
    PWC_butter = filtfilt(B,A,PWC_PWC);
    [psif_filt,Fid_PWC(kk)] = Weak_Dressed_H(PWC_butter,initial_state,target_state,spin,t_final,1,1);
    d_PWC_butter(kk) = abs(psif_filt'*psif_PWC)^2;
    Com_Fid_PWC(kk) = abs(psif_filt'*target_state)^2;
end
PB_FC_out = PB_butter_FC*(2*pi*fs_PWC/2);
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