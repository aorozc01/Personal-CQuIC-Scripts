% This function uses the Filtering_Function.m to calculate the passband
% frequency at which filtering still produces a fidelity of 0.99 or
% greater. Comparison of this value between the constrained and the
% unconstrained phases may give a universal relationship that is
% independent of the dimension.

clear
spin = (1:6)*(1/2);
ratio = zeros(1,length(spin));
for mm = 1:length(spin)
dim = 2*spin + 1;
PB_FC = 1.2*dim+0.8;
PB_PWC = 20*dim -30;
PB_FC_Interval = PB_FC-2:PB_FC+2;
PB_PWC_Interval = max(PB_FC,PB_PWC-40):PB_PWC+40;
for kk = 1:length(PB_FC_Interval)
    [PB_FC_out,d_FC_butter] = Filtering_Function_FC_Only(spin(mm),PB_FC_Interval(kk));
    if d_FC_butter(end) >= 0.99
        wmax_FC = PB_FC_Interval(kk);
        break
    end
end
for ll = 1:length(PB_PWC_Interval)
    [PB_PWC_out,d_PWC_butter] = Filtering_Function_PWC_In_Optim_script(spin(mm),PB_PWC_Interval(ll));
    if d_PWC_butter(end) >= 0.99
        wmax_PWC = PB_PWC_Interval(ll);
        break
    end
end
ratio(mm) = wmax_PWC/wmax_FC;
end