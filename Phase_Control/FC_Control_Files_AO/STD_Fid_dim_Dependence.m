% The purpose for this function is to vary the dimension using optimization
% control function and find the dependence of the std on this parameter.

clear

max_spin = 1;
spin = 2;%0.5:0.5:max_spin;
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
          freq_FC,FC_power_spectrum,freq_PWC,PWC_power_spectrum,...
          FC_dt,PWC_dt,FC_phase,PWC_phase,aopt,...
          PWC_aopt,target_state,t_final,steps_FC,...
          initial_state,psi_f] =...
           Control_Optimization(spin(jj),1,1,1,20,1);     
end