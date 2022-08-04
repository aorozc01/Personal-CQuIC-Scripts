function [Avg_PWC_Fidelity,PWC_Fidelity,PWC_PWC,T_PWC_PWC,PWC_dt,PWC_phase,optimal_PWC_coeff,...
    target_state,t_final,initial_state,psi_f_PWC,Avg_wmax,wmax_PWC,d_PWC_butter] = ...
    Control_Optimization_PWC_Only_filtering(spin,beta,Omega,NumAvgs,n)
% Opimization of a piecewise constant and fourier contrained phase
% function. This function also plots the functions and their power spectra.

rng default % sets random generator to default to reproduce results

%---set parameters
%spin = 1;%input('Spin?\n');
dim = 2*spin+1;
t_final = n*2*pi; %one period of the fundamental frequency [T=(2pi)/w0], w0 = Omega = 1
Avg_wmax = zeros(1,length(spin));
Avg_PWC_Fidelity = zeros(1,length(spin));
wmax_PWC = zeros(NumAvgs,length(spin));
PWC_Fidelity = zeros(NumAvgs,length(spin));
m = 33; 
for ii = 1:length(spin)
initial_state = zeros(1,dim(ii))'; initial_state(1) = 1; % initial state is pointing in postive z-direction
PWC_steps = 2*dim(ii);
initial_PWC_coeff = rand(1,PWC_steps);
%---Parameters used to average results
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
PB_PWC = 20*dim(ii) -30;
PB_PWC_0 = 1;%max(dim(ii),PB_PWC-80);
PB_PWC_f = PB_PWC+50;
for jj = 1:NumAvgs
    target_state= randn_target_state(dim(ii)); % creates a normalized vector with random valued component

    %---Piecewise Constant Phase Optimization
    G = @(PWC_steps) PWC_Phase(PWC_steps,initial_state,target_state,spin(ii),...
        t_final,beta,Omega);
    optimal_PWC_coeff = fminunc(G,initial_PWC_coeff,options); % aopt changes every iteration
    [PWC_Fidelity(jj,ii),~,PWC_phase,PWC_dt,psi_f_PWC] = PWC_Phase...
        (optimal_PWC_coeff,initial_state,...
        target_state,spin(ii),t_final,1,1);  
    
    [PWC_PWC,T_PWC_PWC] = piecewise_repmat_leftmost_point(PWC_phase,PWC_dt,PWC_dt/m); %create piecewise constant function
   
    [wmax_PWC(jj,ii),d_PWC_butter]= Filtering_Function_PWC_In_Optim_Func(initial_state,...
            target_state,PWC_PWC,T_PWC_PWC,t_final,spin(ii),PB_PWC_0,PB_PWC_f);
end
Avg_PWC_Fidelity(ii) = sum(PWC_Fidelity(:,ii))/NumAvgs;
Avg_wmax(ii) =  (1/NumAvgs)*sum(wmax_PWC(:,ii));
end
end