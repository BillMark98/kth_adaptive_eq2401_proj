function [xhat,num,den] = ncw(y, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen)

%
% [xhat,num,den] = ncw(y, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen)
%	
%	y			- y(n)=x(n)+v(n)
% 	PhixyNum,PhixyDen	- Cross-spectrum between x(n) and y(n)
% 	PhiyyNum,PhiyyDen	- Spectrum of y(n)
%	
% 	xhat		- Non-causal Wiener estimate of x(n) from y(n)
% 	num,den		- The Non-causal Wiener filter
%
%  ncw: Non-causal Wiener filtering.
%     
%     
%     Author: Panwei Hu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num = conv(PhixyNum, PhiyyDen);
den = conv(PhixyDen, PhiyyNum);
xhat = ncfilt(num, den, y);