function [out] = mean_data(x,f)
% mean_data(x,f) finds the mean of the data x weighted by the values f.
% Both inputs must be row vectors.
out = x*f';
end