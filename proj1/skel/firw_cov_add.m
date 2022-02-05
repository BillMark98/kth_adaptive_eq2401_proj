function [SigmaYY,SigmaYx] = ...
           firw_cov_add(A, sigma2, Anoise, sigma2noise, N)

%
% [SigmaYY,SigmaYx] = firw_cov_add(A, sigma2, Anoise, sigma2noise, N)
%	
%	A		- AR model for the signal x(n), A(q)x(n)=w(n)
%	sigma2		- E[w(n)*w(n)]
%	Anoise		- AR model for the noise v(n), Anoise(q)v(n)=e(n)
%	sigma2noise	- E[e(n)*e(n)]
%	N  	- Length of Y(n)
%	
% 	SigmaYY		- E[Y(n) (Y(n))']
%	SigmaYx		- E[Y(n) x(n)]
%
%  firw_cov_add: Calculate covariance and cross-covariance for
%     Y(n)=[y(n), y(n-1),...,y(n-N+1)]' where y(n)=x(n)+v(n)
%     
%     Author: Panwei Hu
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y = x + v
% SigmaYx = [rx(0), rx(1),..., rx(N-1)]'
rx_N = ar2cov(A,sigma2,N-1);
SigmaYx = rx_N;
% SigmaYY ,row1 [rx(0) + rv(0), rx(1) + rv(1)...]
rv_N = ar2cov(Anoise,sigma2noise,N-1);
SigmaYY = toeplitz(rx_N + rv_N);