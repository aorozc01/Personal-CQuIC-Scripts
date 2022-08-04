% This function uses the Filtering_Function.m to calculate the passband
% frequency at which filtering still produces a fidelity of 0.99 or
% greater. Comparison of this value between the constrained and the
% unconstrained phases (control waveform) may give a universal relationship that is
% independent of the dimension.

clear
spin = 3;
PB_PWC = 72;%[linspace(69,72,4)];%linspace(50,60,10);linspace(70,80,10)];%...
%     linspace(130,140,10);linspace(130,160,10);linspace(140,170,10);linspace(120,170,10)];
for mm = 1:length(spin)
dim = 2*spin + 1;
max = 2*dim;
PB_FC = dim:max;
for kk = 1:length(PB_FC)
    [PB_FC_out,d_FC_butter,PB_PWC_out,d_PWC_butter] = Filtering_Function(spin(mm),PB_FC(kk),PB_PWC(kk));
    if d_FC_butter(end) >= 0.99
        wmax_FC = PB_FC(kk);
    elseif d_PB_butter(end)
    end
end
end
ratio(mm) = wmax_PWC/wmax_FC;
end