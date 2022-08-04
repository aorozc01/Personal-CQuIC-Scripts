function [ang_freq_Psi,Psi,FT_phiF,PWC_Func,Time,...
    PWC_phiF_fft,df_FC,aoptF] = ...
    piecewise_function_generator(aoptF,phiF,dt,m) 
% This script is a general function that generates a piecewise constant
% function.

%--create piecewise constant phase functions
[PWC_Func,Time] = piecewise_repmat_leftmost_point(phiF,dt,dt/m);

%--calculate Fourier Transform of the piecewise constant phase functions
%with analytic FT calculationFC
PWC_FC_dt = median(diff(Time));%(2)-T_PWC_phiF(1);

% frequency output of FTSetOfPoints is w not f
[ang_freq_Psi,Psi,FT_phiF] = FTSetOfPoints(PWC_Func,PWC_FC_dt);
%--calculate Fourier Transform of the piecewise constant phase functions
%using fft for comparison with function above

PWC_phiF_fft = fft(PWC_Func); % fft

% frequency interval
df_FC = 1/(length(PWC_Func)*PWC_FC_dt);
end