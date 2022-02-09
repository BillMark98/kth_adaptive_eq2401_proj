function [spec] = conjugateFlip(onesidedSpec)
% given a onesided spectrum (0 to pi), extend to 2pi so that 
% corresponding to a real signal
onesidedSpec = columnVector(onesidedSpec)';
if (mod(length(onesidedSpec),2))
    % simply flip conjugate
    spec = [onesidedSpec(1), onesidedSpec(2:end), conj(flip(onesidedSpec(2:end)))];
else
    spec = [onesidedSpec(1), onesidedSpec(2:end-1), real(onesidedSpec(end)), conj(flip(onesidedSpec(2:end-1)))];
end