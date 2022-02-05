function [xhat,theta] = firw(y, SigmaYx, SigmaYY)

%
% [xhat,theta] = firw(y, SigmaYx, SigmaYY)
%	
%	y	- y(n)=x(n)+v(n)
% 	SigmaYY	- E[Y(n) (Y(n))']
%	SigmaYx	- E[Y(n) x(n)]
%	
% 	xhat	- FIR Wiener estimate of x(n) from y(n)
% 	theta	- FIR Wiener filter.
%	
%
%  firw: FIR Wiener estimate of x(n) from y(n)
%     
%     
%     Author: Panwei Hu
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% theta = inv(SigmaYY) SigmaYY
theta = SigmaYY\SigmaYx;
[N,N] = size(SigmaYY);
% topelitz Y matrix, upper triangle 0  assume y(n) = 0 for n leq 0
Gamma_y = toeplitz(y, [y(1), zeros(1, N-1)]);

% xhat = columnVector(y)' * theta;
xhat = Gamma_y * theta;