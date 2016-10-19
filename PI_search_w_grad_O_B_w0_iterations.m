% Evaluation of PI_abs2_GRAPE_Sym_Ryd_Contin_Phase with Gradient 
% This script is mostly used to figure out the relative dependence of the
% parameter w0, beta, and Omega
clearvars -except a0 psitMatrix
tic
spin = 3;%input('Spin?\n');
dim = 2*spin+1;
tf = 2*pi; %with O~2*pi*MHz then tf ~6us
steps = 28;
dt =tf/steps;
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
a0 = rand(1,2*dim); % initial random Fourier Amplitude guess
Omega = 1;
beta = linspace(0,5,10)*Omega;
w0 = linspace(0,5,10)*Omega;
NumOfAvgs= 5; % Number of random target states to average over
FFid = zeros(length(w0),length(beta));
GFid = zeros(1,NumOfAvgs);

for mm = 1:length(w0);
for ii = 1:length(beta);
    for kk = 1:NumOfAvgs
        psitMatrix = NormRandNVector(dim); % creates a normalized vector with random valued components
        F = @(a) PI_abs2_GRAPE_Sym_Ryd_Contin_Phase_2016_09_18_w_dt_B_O_w0_input(a,psi0,psitMatrix,spin,tf,dt,beta(ii),Omega,w0(mm));
        [aoptF,FFval] = fminunc(F,a0,options); % aopt changes every iteration
        [Fid,grad,phi] = PI_abs2_GRAPE_Sym_Ryd_Contin_Phase_2016_09_18_w_dt_B_O_w0_input(aoptF,psi0,psitMatrix,spin,tf,dt,beta(ii),Omega,w0(mm));
        GFid(kk) = Fid;
    end
    FFid(mm,ii) = sum(GFid)/NumOfAvgs;
end
end
figure;plot((0:(steps-1))*dt,phi,'*b')
figure;plot(beta,-FFid,'*r')
title(sprintf('Spin %0.2f, O =  %0.2f, w0 = %0.2f\n',spin,Omega,w0))
ylabel('Fidelity')
xlabel('beta [Omega]')

[W0,Beta] = meshgrid(w0,beta);
surf(W0,Beta,-FFid)

toc 