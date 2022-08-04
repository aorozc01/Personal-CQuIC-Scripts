function [t_samples,samples] = samplegen(func,t_func,dt_samp,m,dt_step)
M = length(func);
T = M*dt_step;
dt_func = dt_step/m;
p = round(dt_samp/dt_func);


samples = zeros(length(M/p),1);
t_samples = zeros(length(M/p),1);
for jj = 1:(M/p);
    samples(jj) = func(jj*p);
    t_samples(jj) = t_func(jj*p);
end
end