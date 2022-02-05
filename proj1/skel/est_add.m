function [xhatnc,xhatc,xhatfir,numnc,dennc,numc,denc,thetahatfir] =...
    est_add(x, v, N, Ahat, sigma2hat, Anoisehat,...
    sigma2noisehat,SigmaYxhat, SigmaYYhat)

%
% [xhatnc,xhatc,xhatfir,numnc,dennc,numc,denc,thetahatfir] = 
%     est_add(x, v, N, Ahat, sigma2hat, Anoisehat,
%       sigma2noisehat,SigmaYxhat, SigmaYYhat)
%	
%	x			 - AR Signal
%	v			 - AR Noise, y(n)=x(n)+v(n)
%	N			 - Length of the FIR Wiener filter
%	Ahat,sigma2hat 		 - Estimated or true parameters of x
%	Anoisehat,sigma2noisehat - Estimated or true parameters of v
%	SigmaYxhat		 - E[Y(n) x(n)]
% 	SigmaYYhat		 - E[Y(n) (Y(n))']
%	
% 	xhatnc		- Non-causal Wiener estimate of x
% 	xhatc		- Causal Wiener estimate of x
% 	xhatfir		- FIR Wiener estimate of x
% 	numnc,dennc	- Non-causal Wiener filter
% 	numc,denc	- Causal Wiener filter
% 	thetahatfir	- FIR Wiener filter
%	
%
%  est_add: Estimate using the three different Wiener filters.
%     Plot the 30 first samples of each estimate, x and y.
%     Calculate the MSE of the estimates normalized with the
%     inverse of the noise variance.
%     
%     Author: Panwei Hu
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% generate y
y = x + v;
[xhatfir, thetahatfir] = firw(y, SigmaYxhat, SigmaYYhat);


% % alternative
% [MySigmaYY, MySigmaYx] = firw_cov_add(A, sigma2, Anoise, sigma2noisehat,N);

% NCW
[PhixyNum,PhixyDen,PhiyyNum,PhiyyDen] = ...
               spec_add(Ahat, sigma2hat, Anoisehat, sigma2noisehat);
[xhatnc, numnc, dennc] = ncw(y, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen);

% default m = 0
m = 0;
[xhatc, numc, denc] = cw(y, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen, m);
