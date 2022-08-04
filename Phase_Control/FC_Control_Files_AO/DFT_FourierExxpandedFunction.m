function a = DFT_FourierExxpandedFunction(aopt,dt)
% This function calculates the discrete Fourier transform of the given set
% of points f with sample point spacing dt. 
M = length(aopt);
c1 = aopt(1:M/2);
c2 = aopt(M/2+1:end);

%{
if length(c1)~=length(c2)
    disp('The input vector does not have an even number of elements')
    return
end
%}

a = zeros(length(c1),1); % vector for complex amplitudes in Fourier Series

for jj = 1:length(c1); % evaluate the modulus squared of the complex amplitudes using the input amplitudes
    a(jj) = abs((c1(jj) + 1i*c2(jj))/2)^2;
end

% The following plots the results from above
%{
k = 1:length(a);
figure;stem(k,a,'b')
hold on;stem(-k,a,'b')
%}