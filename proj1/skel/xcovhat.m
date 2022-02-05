function SigmaYxhat = xcovhat(x, y, N)

%
% SigmaYxhat = xcovhat(x, y, N)
%
%	y, x			- Data sequences
%	N			- Size of covariance matrix
%
%  xcovhat: Estimates SigmaYx=E[Y(n)x(n)]
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
M_1 = length(x);
if (M ~= M_1)
    disp("length unequal, choose the minimu");
end
M = min(M,M_1);
r_xy = zeros(N,1);
for k = 1 : N
    r_xy(k) = sum(y(1:M-k+1) .* x(k : M))/M;
end
% size is N * 1
% Note that E Yx  has entry y(n), y(n-1) ... y(n-N+1) x(n)
% so corresponds to E(y(n-k)x(n) = E(x(n)y(n-k) = rxy(k)
SigmaYxhat = r_xy;
