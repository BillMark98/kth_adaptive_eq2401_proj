function [y] = raisedCosShiftFilter(x, omegaLow,omegaHigh,beta, cutoutN)
% shifted raised cosine

sps = floor(2 * pi / (omegaHigh - omegaLow));
span = ceil(2 * cutoutN / sps);
if (mod(span * sps, 2) == 1)
    span = span + 1
end
rcos = rcosdesign(beta,span,sps);
shiftOmega = (omegaHigh + omegaLow)/2;
% empirical to make pass band 1
rcos = rcos * span;
oneSideLen = (length(rcos)-1)/2;
smp_n = -oneSideLen:oneSideLen;
filter_response = rcos .* cos(shiftOmega * smp_n);
x = columnVector(x)';
y = imfilter(x,filter_response);