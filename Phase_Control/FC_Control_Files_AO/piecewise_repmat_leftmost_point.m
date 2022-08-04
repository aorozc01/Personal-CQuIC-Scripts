function [stack,T] = piecewise_repmat_leftmost_point(phi,dt)
%this function plots the vector of points f as a step function with m (integer) points
%used in between steps that are dt wide. The time interval in which f is plotted is t0 to tf.
%The total number of points used for the plot is m*length(f).

m = 10; % number of points used to create piecewise constant plot of phi
N = length(phi);
T = (0:(1/m):N-(1/m))*dt;
copy1 = repmat(phi(1),m,1);
copy2 = repmat(phi(2:end),m,1);
%copy3 = repmat(phi(end),m/2,1);
copy2 = copy2(:);
stack = [copy1;copy2];
%stack = [copy1;copy2;copy3];
return