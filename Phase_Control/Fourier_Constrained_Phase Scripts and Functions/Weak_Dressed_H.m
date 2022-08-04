function [psi_f,Fidelity] = Weak_Dressed_H(phi,initial_state,target_state,spin,t_final,beta,Omega)
% all vector inputs are assumed to be column vectors! 
% Hamiltonian being considered 
%   H = B*J3^2 + O*(J1*cos(phi(t))+J2*sin(phi(t)))

% constants
[J1,J2,J3] = make_fs(spin); % create angular momentum matrices in symmetric basis
dt = t_final/(length(phi)); % time step width

% initial condition
psi_f = initial_state; % first entry of the set above is the initial state

% normalize the input wave functions
norm_target = target_state'*target_state;
norm_initial = initial_state'*initial_state;
if abs(norm_target-1) >= 1e-08
    target_state = target_state/sqrt(norm_target);
    disp('Target state not normalized')
end
if norm_initial ~= 1
    initial_state = initial_state/sqrt(norm_initial);
    disp('Initial state not normalized')
end



%--------------This portion calculates the fidelity--------------------

% create U_tot by creating U_j's where j is the time step
for jj = 1:length(phi) % This repeats for all time steps
        H_tot = beta*(J3^2)+Omega*(cos(phi(jj))*J1 + sin(phi(jj))*J2); % Hamiltonian evaluated for each time step
        U_jth_step = expm(-1i*H_tot*dt);
        psi_f = U_jth_step*psi_f; % generate time final states for each step
end % end of main for loop

clear ll kk mm pp  % clear for loop variables
    trace_term = conj(target_state'*psi_f);
    Fidelity = -abs(trace_term)^2; % the negative sign is used to find the maximum fidelity using fminunc which finds the minumum of a function
   
end