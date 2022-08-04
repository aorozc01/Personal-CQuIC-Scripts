function [freq_pos1,SS1,freq_pos2,SS2] = TSS_to_SSS(freq1,spectrum1,freq2,spectrum2)
% TSS_to_SSS(freq1,spectrum1,freq2,spectrum2) converts the freq and
% spectrum inputs froma two-sided spectrum into a single-sided spectrum;
m1 = (length(freq1)-1)/2;
n1 = (length(spectrum1)-1)/2;
freq_pos1 = freq1(m1+1:end); % includes postive freq and DC component
spectrum_DC1 = spectrum1(m1+1);
spectrum_pos1 = 2*spectrum1(m1+2:end); 
SS1 = [spectrum_DC1;spectrum_pos1]/max(spectrum_pos1);

if m1~=n1
    disp('check vector lengths are equal')
    return
end

if nargin > 2
m2 = (length(freq2)-1)/2;
n2 = (length(spectrum2)-1)/2;
freq_pos2 = freq2(m2+1:end);
spectrum_DC2 = spectrum1(m2+1);
spectrum_pos2 = 2*spectrum2(m2+2:end);
SS2 = [spectrum_DC2;spectrum_pos2]/max(spectrum_pos2);
end

end