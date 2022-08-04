function out = variance(x,f,mu)
% variance(x,f,mu)finds the variance of the data x with average mu weighted 
% by the values f. Both x and f must be column vectors and mu a scalar.
out = ((x-mu).^2)'*f;
end
