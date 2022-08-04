function [out,A,D] = phase_filter(data,PB,SB)
% phase_filter(F,fc,dt) filters F data cuttoff frequency fc and time steps dt
% unsure if fc should be multiplied by 2*pi

% data creation
% The phase is the function recieved from the control optimization. It
% needs to be created as a piecewise constant function. The interval
% dt is the spacing between piecewise constant steps and is not related to
% the sampling interval or sampling frequency needed for filtering

% filter
A = designfilt('lowpassfir', 'PassbandFrequency', PB, 'StopbandFrequency', ...
    SB,'PassbandRipple', 0.1, 'StopbandAttenuation', 80);%, 'SampleRate', fs);
% compensate for constant digital filter delay
D = floor(median(grpdelay(A)));
out = filter(A,[data;zeros(D,1)]);
end