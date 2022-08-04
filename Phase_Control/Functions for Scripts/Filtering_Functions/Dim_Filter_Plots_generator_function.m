% This script used the filtering function to generate a contrasting plot
% between the filtering effects of the phase for different dimensional
% hilbert spaces.


spin = 5:7;%(1:6)/2;
PB_FC_final = 20;
PB_PWC_final = 180;
% figure; grid on;

for mm = 1:length(spin)
    [PB_FC,d_FC,PB_PWC,d_PWC] = Filtering_Function(spin(mm),PB_FC_final,PB_PWC_final);
%     hold on;plot(PB_FC,d_FC,'o',PB_PWC,d_PWC,'o')
end
for kk = 1:length(d_FC)
    if d_FC(kk)>=0.99
        FC_wmax = kk;
        break
    end
end
for ll = 1:length(d_PWC)
    if d_PWC(ll)>=0.99
        PWC_wmax = ll;
        break
    end
end
% xlabel('\omega_0')
% ylabel('Power spectrum')
