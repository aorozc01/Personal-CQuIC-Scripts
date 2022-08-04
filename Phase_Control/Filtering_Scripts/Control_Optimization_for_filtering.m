function [Avg_FC_Fidelity,Avg_PWC_Fidelity,freq_PWC,PWC_power_spectrum,...
    PWC_FC,T_PWC_FC,PWC_PWC,T_PWC_PWC,FC_dt,PWC_dt,FC_FT,PWC_FT,FC_phase,...
    PWC_phase,optimal_coeff,PWC_phiF_fft,PWC_phiG_fft,df,optimal_PWC_coeff,...
    target_state,t_final,steps_FC,initial_state,psi_f_FC,psi_f_PWC] = ...
    Control_Optimization_for_filtering(spin,beta,Omega,w0,NumAvgs,plots,n)
% Opimization of a piecewise constant and fourier contrained phase
% function. This function also plots the functions and their power spectra.

rng default % sets random generator to default to reproduce results

%---set parameters
%spin = 1;%input('Spin?\n');
dim = 2*spin+1;
t_final = n*2*pi; %one period of the fundamental frequency [T=(2pi)/w0], w0 = Omega = 1
initial_state = zeros(1,dim)'; initial_state(1) = 1; % initial state is pointing in postive z-direction
initial_coeff = rand(1,2*dim); % initial random Fourier Amplitude guess
steps_FC  = n*100; %used in FC functions to determing dt size
%  A sampling rate is not needed here because we are constructing optimal
%  phase function to create the actual functions we are going to sample
%  from.
PWC_steps = 2*dim;

initial_PWC_coeff = randn(1,PWC_steps);

%---Parameters used to average results
FC_Fidelity = zeros(1,NumAvgs);
PWC_Fidelity = zeros(1,NumAvgs);

options = optimset(...  % these are the settings from unitary control search
    'TolX',            1e-16,...    % related to minimum tolerance for change in x
    'TolFun',          1e-8,...    % minumum tolerance for change in the function
    'MaxIter',         2000,...   % maximum number of iterations
    'DerivativeCheck', 'off',...    % compare analytic gradient to numerical estimation (off)
    'GradObj',         'on',...         % tells matlab whether the gradient is supplied to fminunc
    'LargeScale',       'off', ... % when turned off will increase calculation speed
    'Display',         'off',...          % output type in matlab command window, USE 'iter' to turn on
    'MaxFunEvals',     10^6,...         % maximum number of function calls'MaxFunEvals', 500);
    'ObjectiveLimit',  -0.99999);    % once objective function reaches this value, fminunc ...
                                        %stops (ONLY works for quasi-Newton method)
for jj = 1:NumAvgs
    target_state= randn_target_state(dim); % creates a normalized vector with random valued component
    
    %---Fourier Contrained Phase Optimization
    F = @(coeff) FourierConstrainedPhase(coeff,initial_state,target_state,...
        spin,t_final,steps_FC,beta,Omega,w0);
    optimal_coeff = fminunc(F,initial_coeff,options); % aopt changes every iteration  
    [FC_Fidelity(jj),~,FC_phase,FC_dt,psi_f_FC] = FourierConstrainedPhase...
        (optimal_coeff,initial_state,...
        target_state,spin,t_final,steps_FC,1,1,1);   

    %---Piecewise Constant Phase Optimization
    G = @(PWC_steps) PWC_Phase(PWC_steps,initial_state,target_state,spin,...
        t_final,beta,Omega);
    optimal_PWC_coeff = fminunc(G,initial_PWC_coeff,options); % aopt changes every iteration
    [PWC_Fidelity(jj),~,PWC_phase,PWC_dt,psi_f_PWC] = PWC_Phase...
        (optimal_PWC_coeff,initial_state,...
        target_state,spin,t_final,1,1);   
end
Avg_FC_Fidelity = sum(FC_Fidelity)/jj;
Avg_PWC_Fidelity = sum(PWC_Fidelity)/jj;clear jj

%---Plot the spectrum of the piecewise constant phase
[~,~,freq_PWC,PWC_power_spectrum,FC_FT,PWC_FT,PWC_FC,T_PWC_FC,...
    PWC_PWC,T_PWC_PWC,PWC_phiF_fft,PWC_phiG_fft,df,aopt] = ...
    Spectrum_Plots(optimal_coeff,FC_phase,...
    PWC_phase,t_final/steps_FC,t_final/PWC_steps,spin,plots,dim);
end