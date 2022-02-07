function [y] = moved_idealLP(x, omegaLow, omegaHigh, cutoutN)
% filter the signal x with a shifted ideal LP
% x - input signal
% omegaLow - the lower cutoff Freq
% omegaHigh - the higher cutoff Freq
% cutoutN - the cut out N in the time domain for the ideal pass filter

Omega = (omegaHigh - omegaLow)/2;
shiftOmega = (omegaHigh + omegaLow)/2;
smp_n = -cutoutN : cutoutN;
% filter_response = sin(Omega * smp_n)./(pi * smp_n + eps) .* cos(shiftOmega * smp_n);
filter_response = 2 * Omega/pi * sinc(Omega/pi * smp_n) .* cos(shiftOmega * smp_n);
x = columnVector(x)';
% image_x = 
y = imfilter(x,filter_response);
