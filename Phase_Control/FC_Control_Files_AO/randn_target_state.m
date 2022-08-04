function psit = randn_target_state(dim)
% randn_target_state(dim) creates a vector with number of rows equal to dim. 
% The components are chosen from a normalized gaussian distribution.
vector = randn(dim,1);
psit = vector/sqrt(abs(vector'*vector));
end