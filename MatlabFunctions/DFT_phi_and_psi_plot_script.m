function DFT_phi_and_psi_plot_script(aoptF,phiF,phiG,dt,dtG) 
%This script uses the analytical solution to the fourier transform for a
% piecewise constant function and a fourier expanded function. To be able
% to use the scritp first run a program like
% PI_search_w_grad_single_iteration or similar script

% VARIABLES NEEDED: {aopt,psi,dt}
[psi,T] = piecewise_repmat_leftmost_point(phiF,dt);
[psiG,TG] = piecewise_repmat_leftmost_point(phiG',dtG);
figure;
[wpsi,Fpsi] = FTSetOfPoints(psi,abs(T(2)-T(1)));
[wGpsi,Gpsi] = FTSetOfPoints(psiG,abs(TG(2)-TG(1)));
PF = length(phiF);
PG = length(phiG);
plot((0:PF-1)*dt,phiF,'*r',T,psi,'b',(0:PG-1)*dtG,phiG,'*k',TG,psiG,'m')
title('Phase as a function of time for Spin 3')
legend('Sample Points','Fourier Constrained Phase','Points','Peicewise Constant')
xlabel('time [1/Rabi]');ylabel('Phase Amplitude')
a = DFT_FourierExxpandedFunction(aoptF,dt);
k = 1:length(a);
figure;
stem(k,a/max(a),':r','Linewidth',1)
hold on;stem(-k,a/max(a),':r','linewidth',1)
plot(wpsi,Fpsi,'-b',wGpsi,Gpsi,'-m')
title('Fourier Spectrum')
legend('a','b','c','d')
ylabel('Amplitude');xlabel('frequency [Rabi Frequency]')
set(gca,'XTick',-40:5:40)
axis([-50,50,0,1])
figure;
plot(wpsi,Fpsi,'-b',wGpsi,Gpsi,'-m')
title('Fourier Spectrum')
legend('Fourier Constrained','Piecewise Constant')
ylabel('Amplitude');xlabel('frequency [Rabi Frequency]')
set(gca,'XTick',-40:5:40)
axis([-20,20,0,1])
end