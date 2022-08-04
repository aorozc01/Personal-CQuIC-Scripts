function [wpsi,Fpsi,wGpsi,Gpsi] = Spectrum_Plots(aoptF,phiF,phiG,dt,dtG) 
% This script uses the analytical solution to the fourier transform for a
% piecewise constant function and a fourier expanded function. To be able
% to use the script first run a program like
% PI_search_w_grad_single_iteration or similar script

% VARIABLES NEEDED: {aopt,psi,dt}

%---calculate data for piecewise constant plot
[psi,T] = piecewise_repmat_leftmost_point(phiF,dt);
[psiG,TG] = piecewise_repmat_leftmost_point(phiG,dtG);

%---calculate piecewise constant DFT
[wpsi,Fpsi] = FTSetOfPoints(psi,abs(T(2)-T(1)));
[wGpsi,Gpsi] = FTSetOfPoints(psiG,abs(TG(2)-TG(1)));
a = DFT_FourierExxpandedFunction(aoptF,dt);

%---plot piecewise constant functions
k = 1:length(a);
PF = length(phiF);
PG = length(phiG);
figure;plot((0:PF-1)*dt,phiF,'*r',T,psi,'b',(0:PG-1)*dtG,phiG,'*k',TG,psiG,'m')
title('Phase as a function of time for Spin 3')
legend('Sample Points','Fourier Constrained Phase','Points','Peicewise Constant')
xlabel('time [1/Rabi]');ylabel('Phase Amplitude')

%---plot fourier coefficients and Fourier constrained piecwise constant spectrums
figure;subplot(1,2,1)
stem(k,a/max(a),':r','Linewidth',1)
hold on;stem(-k,a/max(a),':r','linewidth',1)
plot(wpsi,Fpsi,'-b')
title('Fourier Spectrum')
legend('Fourier Coefficients','FC Phase')
ylabel('Amplitude');xlabel('frequency [Rabi Frequency]')
set(gca,'XTick',-40:5:40)
axis([-20,20,0,1])

%---plot spectrum of piecewise constants only
subplot(1,2,2)
plot(wpsi,Fpsi,'-b',wGpsi,Gpsi,'-m')
title('Fourier Spectrum')
legend('Fourier Constrained','Piecewise Constant')
ylabel('Amplitude');xlabel('frequency [Rabi Frequency]')
set(gca,'XTick',-40:5:40)
axis([-20,20,0,1])
end