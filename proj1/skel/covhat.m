function SigmaYYhat = covhat(y, N)

%
% SigmaYYhat = covhat(y,N)
%
%	y			- Data sequence
%	N			- Size of covariance matrix
%
%  covhat: Estimates SigmaYY=E[Y(n)Y^{T}(n)]
%
%		where 
%
%	   	Y(n)=[y(n) y(n-1) ... y(n-N+1)]^{T}
%
%     
%     Author: Panwei Hu
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the length of signal
M = length(y);
r_yy = zeros(N,1);
for k = 1 : N
    r_yy(k) = sum(y(1:M-k+1) .* y(k : end))/M;
end
SigmaYYhat = toeplitz(r_yy);