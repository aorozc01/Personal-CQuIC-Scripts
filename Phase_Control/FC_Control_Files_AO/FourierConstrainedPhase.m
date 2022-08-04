
function [Fidelity,Gradient,phase,dt,psi_f] = FourierConstrainedPhase(coeff,initial_state,target_state,spin,t_final,steps,beta,Omega,w0)
% all vector inputs are assumed to be column vectors! 
% This function is called by FourierPhase_ArbSpin_ParforIMP.m to solve for the fidelity
% function and also the gradient to be used with fminunc with 1st order RWA
%
% Hamiltonian being considered 
%   H = B*J3^2 + O*(J1*cos(phi(t))+J2*sin(phi(t)))

% constants
[J1,J2,J3] = make_fs(spin); % create angular momentum matrices in symmetric basis
dim = 2*spin+1; % Hilbert space dimension
w_vector = w0*(1:1:dim); % vector of harmonics
length_coeff = length(coeff);

dt = t_final/steps; % time step width
phase = zeros(1,steps); % row vector of zeros for phi

% initializes unitaries to be used in GRAPE algorithm
U_tot = eye(dim);
U_jth_step = zeros(dim,dim,steps+1); % create dim x dim matrices T_length times
U_jth_step(:,:,steps+1) = eye(dim);
dU_dC = zeros(dim,dim, steps, 2*dim);

Gradient = zeros(length_coeff,1);

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

psi_f = initial_state; % first entry of the set above is the initial state

%--------------This portion calculates the fidelity--------------------

% create U_tot by creating U_j's where j is the time step
Amp_constraint = 1/(2*pi);
for jj = 1:steps % This repeats for all time steps
        phase(jj) = Amp_constraint*(coeff(1,1:length_coeff/2)*(sin(w_vector*jj*dt))' + coeff(1,length_coeff/2+1:length_coeff)*(cos(w_vector*jj*dt))'); % Fourier Sum for the phase for each time step
        H_tot = beta*(J3^2)+Omega*(cos(phase(jj))*J1 + sin(phase(jj))*J2); % Hamiltonian evaluated for each time step
       
        % eigendecomposition of H_tot
        [VH,DH] = eig(H_tot); % VH is a unitary that takes the eigenbasis to the symmetric basis, Dh is a matrix is with diagonal elements=eigenvalues of H_tot
        D_vec = diag(DH); % creates a vector with components = to the eigenvalues of H_tot
        VH_ct = ctranspose(VH);
   
        exp_eig_factor = -1i*dt*diag(exp(-1i*dt*D_vec)); % matrix with diagonal values equal to eigenvalues of H_tot
        
        % solve for derivative of U w.r.t control fields (in H_tot eigenbasis)
        
        for pp = 1:dim % Only creates off diagonal elements
            for mm = (pp + 1):dim
                exp_eig_factor(pp,mm) = ( exp(-1i*dt*D_vec(pp)) - exp(-1i*dt*D_vec(mm)) ) / ( D_vec(pp)-D_vec(mm) );
                exp_eig_factor(mm,pp) = exp_eig_factor(pp,mm); % this matrix is symmetric
            end
        end
        clear pp
        % derivative of H with respect to each control field
        
        for pp = 1:length_coeff/2; % pp index corresponds to parameter subscript in Fourier expansion 
            dH_A = Omega*(-J1*sin(phase(jj))+J2*cos(phase(jj)))*sin(w_vector(pp)*jj*dt)*Amp_constraint; % these matrices are written in the symmetric basis of H_tot
            dH_B = Omega*(-J1*sin(phase(jj))+J2*cos(phase(jj)))*cos(w_vector(pp)*jj*dt)*Amp_constraint; % these matrices are written in the symmetric basis
        
        
        % derivative of the unitaries and transform them to H_tot eigenbasis (O is the rabi frequency)
            dU_dC(:,:,jj,pp) = VH * ((VH_ct*dH_A*VH).*exp_eig_factor) * VH_ct; % Derivative for two different parameters: the fourier coefficients A and B
            dU_dC(:,:,jj,pp+length_coeff/2) = VH * ((VH_ct*dH_B*VH).*exp_eig_factor) * VH_ct; % There are 2d total parameter parameters A_k and B_k for k=1:d
        end
        U_jth_step(:,:,jj) = VH*diag(exp(-1i*dt*D_vec))*VH_ct; % calculate U_tot in eigenbasie and then transform back from H_tot eigenbasis to symmetric basis
        psi_f(:,:,jj + 1) = U_jth_step(:,:,jj)*psi_f(:,:,jj); % generate time final states for each step
        
end % end of main for loop

clear ll kk mm pp  % clear for loop variables
    trace_term = conj(target_state'*psi_f(:,:,steps+1));  % step+1 is used because the initial state is stored as first entry
    Fidelity = -abs(trace_term)^2; % the negative sign is used to find the maximum fidelity using fminunc which finds the minumum of a function

    %Construct total Unitary
for jj = 1:steps;
        jp = steps - jj + 1;
        U_tot = U_tot*U_jth_step(:,:,jp); % construct total unitary
end    
    
%-----------------This portion calculates the gradient---------------------

if nargout > 1
        % solve for psi_p's
        for pp = 1:length_coeff; % number of parameters
            %W = eye(dim);
            psiT = target_state';
            %sumG = zeros(dim,dim);
            sumG = 0;
            for jj = 1:steps % for loop for each time step used to calculate each U_tot derivative term
                jp = steps - jj +1; % trick used to switch starting point from 1 to T_length to T_length to 1!
                psiT = psiT*U_jth_step(:,:,jp+1);
                G = psiT*dU_dC(:,:,jp,pp)*psi_f(:,:,jp);
                sumG = G + sumG;
            end 
            %gradF(pp) = -2*real(trace_term*psi_tf'*sumG(:,:)*psi0);
            Gradient(pp) = -2*real(trace_term*sumG);
        end
         
end