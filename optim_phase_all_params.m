function [Fid_Avg,GFid_Avg] = optim_phase_all_params(spin,tf,PF,PC,NumOfAvgs,B,O,w0,plotfunc,opts)
% tf final time, PF number of steps for Fourier expanded phase, PG number of
% steps for piecewise constant phase, NumOfAvgs is the number of random final 
% target states we average over
% This script is used to find the fidelity of using the Fourier constrained
% and the piecewise constant phase
tic
dim = 2*spin+1; 
L_PF = length(PF);
L_PC = length(PC);
Fid_Avg = zeros(L_PF,1);
GFid_Avg = zeros(L_PC,1);
Fid_Vector = zeros(NumOfAvgs,1);
GFid_Vector = zeros(NumOfAvgs,1);


for ll = 1:L_PF % iterates over different number of steps the user specifies
dt = tf/PF(ll);
options = optimset(...  % these are the settings from unitary control search
    'TolX',            1e-16,...    % related to minimum tolerance for change in x
    'TolFun',          1e-18,...    % minumum tolerance for change in the function
    'MaxIter',         2000,...   % maximum number of iterations
    'DerivativeCheck', 'off',...    % compare analytic gradient to numerical estimation (off)
    'GradObj',         'on',...         % tells matlab whether the gradient is supplied to fminunc
    'LargeScale',       'off',...
    'OutputFcn',    @searchStopFn,...   % stopping function (stops when f(x) < 10^-8)
    'Display',         'iter',...          % output type in matlab command window, USE 'iter' to turn on
    'MaxFunEvals',     10^6);%,...         % maximum number of function calls'MaxFunEvals', 500);
    %'ObjectiveLimit',  -0.99999);    % once objective function reaches this value, fminunc stops (ONLY works for quasi-Newton method)

psi0 = zeros(1,dim)'; % initial state is pointing in postive z-direction
psi0(1) = 1;
for jj = 1:NumOfAvgs;
    a0 = rand(1,2*dim); % initial random Fourier Amplitude guess
    psitMatrix= NormRandNVector(dim); % creates a normalized vector with random valued component
    
    F = @(a) PI_abs2_GRAPE_Sym_Ryd_Contin_Phase_2016_09_18_w_dt_B_O_w0_input(a,psi0,psitMatrix,spin,tf,dt,B,O,w0);
    [aoptF,FFval] = fminunc(F,a0,options); % aopt changes every iteration
    [Fid,grad,phiF] = PI_abs2_GRAPE_Sym_Ryd_Contin_Phase_2016_09_18_w_dt_B_O_w0_input(aoptF,psi0,psitMatrix,spin,tf,dt,B,O,w0);

    Fid_Vector(jj) = Fid;
end
Fid_Avg(ll) = sum(Fid_Vector)/jj;
end
for ll = 1:L_PC % iterates over different number of steps the user specifies
options = optimset(...  % these are the settings from unitary control search
    'TolX',            1e-16,...    % related to minimum tolerance for change in x
    'TolFun',          1e-18,...    % minumum tolerance for change in the function
    'MaxIter',         2000,...   % maximum number of iterations
    'DerivativeCheck', 'off',...    % compare analytic gradient to numerical estimation (off)
    'GradObj',         'on',...         % tells matlab whether the gradient is supplied to fminunc
    'LargeScale',       'off',...
    'OutputFcn',    @searchStopFn,...   % stopping function (stops when f(x) < 10^-8)
    'Display',         'iter',...          % output type in matlab command window, USE 'iter' to turn on
    'MaxFunEvals',     10^6);%,...         % maximum number of function calls'MaxFunEvals', 500);
    %'ObjectiveLimit',  -0.99999);    % once objective function reaches this value, fminunc stops (ONLY works for quasi-Newton method)

psi0 = zeros(1,dim)'; % initial state is pointing in postive z-direction
psi0(1) = 1;
for jj = 1:NumOfAvgs;
    wf = randn(PC,1); % initial guess for contants
    psitMatrix= NormRandNVector(dim); % creates a normalized vector with random valued component
    
    G = @(w) PI_PWC_Phase_all_params(w,psi0,psitMatrix,spin,tf,B,O);
    [boptG,GFval] = fminunc(G,wf,options); % aopt changes every iteration
    [GFid,Ggrad,phiG,dtG] = PI_PWC_Phase_all_params(boptG,psi0,psitMatrix,spin,tf,B,O);
    GFid_Vector(jj) = GFid;
end
GFid_Avg(ll) = sum(GFid_Vector)/jj;
end
if plotfunc ==1;
    DFT_phi_and_psi_plot_script(aoptF,phiF,phiG,dt,dtG) 
end
toc 
end