function [out,A] = phase_filter(data,PB,SB)
% phase_filter(F,fc,dt) filters F data cuttoff frequency fc and time steps dt
% unsure if fc should be multiplied by 2*pi

% data creation
% The phase is the function recieved from the control optimization. It
% needs to be created as a piecewise constant function. The interval
% dt is the spacing between piecewise constant steps and is not related to
% the sampling interval or sampling frequency needed for filtering

% filter
A = designfilt('lowpassiir', 'PassbandFrequency', PB, 'StopbandFrequency', ...
    SB,'PassbandRipple', 0.1, 'StopbandAttenuation', 80);%, 'SampleRate', fs);
% compensate for frequency dependent digital filter delay using filtfilt
% function
out = filtfilt(A,data);
end