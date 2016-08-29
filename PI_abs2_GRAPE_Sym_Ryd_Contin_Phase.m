function [fid,grad_tot] = PI_abs2_GRAPE_Sym_Ryd_Contin_Phase(a,psi0,psi_tf,J,tf)
% all vector inputs are assumed to be column vectors! 
% This function is called by FourierPhase_ArbSpin_ParforIMP.m to solve for the fidelity
% function and also the gradient to be used with fminunc with 1st order RWA
%
% Hamiltonian being considered 
%   H = B*J3^2 + O*(J1*cos(phi(t))+J2*sin(phi(t)))

% constants
[J1,J2,J3] = make_fs(J); % create angular momentum matrices in symmetric basis
dim = 2*J+1; % Hilbert space dimension
O = 1; % Rabi Frequency
B = O; % Entangling coupling stregth
w0 = O; % fundamental frequency for Fourier Expansion
w_vector = w0*(1:1:dim); % vector of harmonics

dt = 0.001; % time step width
T = 0:dt:tf; 
T_length = length(T); % total number of time steps
phi = zeros(1,T_length); % row vector of zeros for phi

% initializes unitaries to be used in GRAPE algorithm
U_tot = zeros(dim,dim,T_length); % create dim x dim matrices T_length times
U_tot(:,:,1) = eye(dim);% initial unitary is the identity
dU_dA = zeros(dim,dim, T_length);% create set of dim x dim matrices T_length times
dU_dB = zeros(dim,dim, T_length);% "    "    "

dH_J1_A = zeros(dim,dim,T_length);
dH_J1_B = zeros(dim,dim,T_length);
dH_J2_A = zeros(dim,dim,T_length);
dH_J2_B = zeros(dim,dim,T_length);

grad = zeros(T_length,2,dim);
grad_tot = zeros(2*dim,1); % initial value

% normalize the input wave functions
norm_psi_tf = psi_tf'*psi_tf;
norm_psi0 = psi0'*psi0;
if norm_psi_tf ~=1
    psi_tf = psi_tf/sqrt(norm_psi_tf);
end
if norm_psi0 ~= 1
    psi0 = psi0/sqrt(norm_psi0);
end

psi_f = zeros(dim,1,(T_length + 1)); % T_length + 1 column vectors to store consecutively evolved final states
psi_f(:,:,1) = psi0; % first entry of the set above is the initial state

%--------------This portion calculates the fidelity--------------------

% create U_tot by creating U_j's where j is the time step
for jj = 1:T_length % This repeats for all time steps
        phi(jj) = a(1,:)*(sin(w_vector*T(jj)))' + a(2,:)*(cos(w_vector*T(jj)))'; % Fourier Sum for the phase for each time step
        H_tot = B*((J3)^2)+O*(cos(phi(jj))*J1 + sin(phi(jj))*J2); % Hamiltonian evaluated for each time step
       
        % eigendecomposition of H_tot
        [VH,DH] = eig(H_tot); % VH is a unitary that takes the eigenbasis to the symmetric basis, Dh is a matrix is with diagonal elements=eigenvalues of H_tot
        D_vec = diag(DH); % creates a vector with components = to the eigenvalues of H_tot
        VH_ct = ctranspose(VH);
   
        exp_eig_factor = -1i*dt*diag(exp(-1i*dt*D_vec)); % matrix with diagonal values equal to eigenvalues of H_tot
        
        % solve for derivative of U w.r.t control fields (in H_tot eigenbasis)
        
        for kk = 1:dim % Only creates off diagonal elements
            for mm = (kk + 1):dim
                exp_eig_factor(kk,mm) = ( exp(-1i*dt*D_vec(kk)) - exp(-1i*dt*D_vec(mm)) ) / ( D_vec(kk)-D_vec(mm) );
                exp_eig_factor(mm,kk) = exp_eig_factor(kk,mm); % this matrix is symmetric
            end
        end
        
        % derivative of H with respect to each control field
        
        for pp = 1:length(w_vector); % pp index corresponds to parameter subscript in Fourier expanision 
            dH_J1_A(:,:,pp) = -J1*sin(phi(jj))*sin(w_vector(pp)*T(jj)); % these matrices are written in the symmetric basis of H_tot
            dH_J1_B(:,:,pp) = -J1*sin(phi(jj))*cos(w_vector(pp)*T(jj)); % these matrices are written in the symmetric basis
            dH_J2_A(:,:,pp) = J2*cos(phi(jj))*sin(w_vector(pp)*T(jj));
            dH_J2_B(:,:,pp) = J2*cos(phi(jj))*cos(w_vector(pp)*T(jj));
        end
        
        % derivative of the unitaries and transform them to H_tot eigenbasis (O is the rabi frequency)
        dU_dA(:,:,jj) = VH * ((VH_ct*(O*(dH_J1_A+dH_J2_A))*VH).*exp_eig_factor) * VH_ct; % Derivative for two different parameters: the fourier coefficients A and B
        dU_dB(:,:,jj) = VH * ((VH_ct*(O*(dH_J1_B+dH_J2_B))*VH).*exp_eig_factor) * VH_ct; % There are 2d total parameter parameters A_k and B_k for k=1:d
        
        U_tot(:,:,jj) = (VH*diag(exp(-1i*dt*D_vec))*VH_ct); % calculate U_tot and then transform back from H_tot eigenbasis to symmetric basis
        psi_f(:,:,jj + 1) = U_tot(:,:,jj)*psi_f(:,:,jj); % generate time final states for each step

end % end of main for loop

clear ll kk mm   % clear for loop variables
    trace_term_U_tot = conj(psi_tf'*psi_f(:,:,T_length+1));
    fid = abs(trace_term_U_tot)^2; % the negative sign is used to find the maximum fidelity using fminunc which finds the minumum of a function
    
%-----------------This portion calculates the gradient---------------------

if nargout > 1
        psi_p = zeros(1, dim, T_length + 1); % bra (row vectors) states to be used in calculation
        psi_p(:,:,T_length + 1) = psi_tf'; % the last vector is the function input target state psi_tf
        
        % solve for psi_p's
        for kk = 1:dim;
            for jj = 1:T_length % for loop for each time step used to calculate each U_tot derivative term
                % jp is the number of control time steps
                jp = T_length - jj + 1; % trick used to switch starting point from 1 to T_length to T_length to 1!
                psi_p(:,:,jp) = psi_p(:,:,jp + 1)*U_tot(:,:,jp); % creates consecutive bras       

                % create control wave form update method
                grad(jp, 1,k) = -2*(real(abs(psi_p(:,:,jp+1)*dU_dA(:,:,jp)*psi_f(:,:,jp))*trace_term))+grad(jp-1,1); % vector of fidelities for A coefficients
                grad(jp, 2,k) = -2*real(abs(psi_p(:,:,jp+1)*dU_dB(:,:,jp)*psi_f(:,:,jp))*trace_term); % vector of fidelities for B coefficients
            end  
        end
        
        g=size(grad);v = zeros(dim,dim);
        for oo = 1:g(end);
            grad_tot(oo,1) = grad(:,:,oo)+ v;
            grad_tot(2*oo,1) = grad(:,:,2*oo)+ v;
        end
end