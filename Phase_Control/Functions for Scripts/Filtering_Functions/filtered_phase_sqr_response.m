function f_filt = filtered_phase_sqr_response(phase,wp)
Lf = length(phase);
f_sum = zeros(length(phase),1);
f_filt = f_sum;
for kk = 1:Lf
    f_sum = zeros(length(phase),1);
for jj = 1:Lf
    f_sum(jj) = f_sum(jj) + (2/sqrt(2*pi))*phase(jj)*wp*sinc(wp*(kk-jj));
end
f_filt(kk) = sum(f_sum);
end

k = 1:length(phase);
j = 0:(length(phase)-1);
[K,J] = m
n = 0:(length(phase)-1)*0.06;
O = (2/sqrt(2*pi))*wp*sinc(wp*n);
A = conv(phase,O,'same');