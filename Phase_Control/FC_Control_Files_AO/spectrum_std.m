function [std_spectrum1,std_spectrum2]=spectrum_std(freq1,power_spectrum1,freq2,power_spectrum2)
% This function calculates the standard deviation (std) of the
% power_spectrum as a function of frequency using the internal matlab
% function var. var(x,f(x)) calculates the variance of the sequence x 
% weighted by f(x). See var for more information on the function.
std_spectrum1 = sqrt(var(freq1,power_spectrum1));
std_spectrum2 = sqrt(var(freq2,power_spectrum2));
end
