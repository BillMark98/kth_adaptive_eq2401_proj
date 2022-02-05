function filterCompare(N, A, sigma2, Anoise,...
    sigma2noise, signalLen, N_sample)
% compare filters
%   N   - FIR filter length
%   A   - AR coeff for x
%   Anoise  AR coeff for v
%   signalLen  the generated signal length
%   N_sample the plotted sample length
K = signalLen;
% random noise
e = sqrt(sigma2) * randn(K,1);

x = filter(1,A,e);
w = sqrt(sigma2) * randn(K,1);
v = filter(1,Anoise,w);

y = x + v;

[SigmaYY, SigmaYx] = firw_cov_add(A,sigma2,Anoise, sigma2noise,N);

[xhatnc,xhatc,xhatfir,numnc,dennc,numc,denc,thetahatfir] =...
    est_add(x, v, N, A, sigma2, Anoise,...
    sigma2noise,SigmaYx, SigmaYY);

% calculate mse
% mse_fir = myMSE(xhatfir, x, 1/sigma2noise);
% mse_nc = myMSE(xhatnc, x, 1/sigma2noise);
% mse_c = myMSE(xhatc, x, 1/sigma2noise);

mse_fir = myMSE(xhatfir, x, sigma2noise);
mse_nc = myMSE(xhatnc, x, sigma2noise);
mse_c = myMSE(xhatc, x, sigma2noise);
% legends
fir_legend = sprintf("FIR(%.4f)", mse_fir);
nc_legend = sprintf("Non-causal(%.4f)", mse_nc);
c_legend = sprintf("Causal(%.4f)", mse_c);
% plot the first 30 samples


indices = 1:1:N_sample;
figure;
plot(indices, y(indices), 'b');
hold on;
plot(indices, x(indices), ':k');
plot(indices, xhatnc(indices),'--.r');
plot(indices, xhatc(indices), '--g');
plot(indices, xhatfir(indices),'-m');

legend("Observations y", "Signal", nc_legend, c_legend, fir_legend);

% compare spectra
figure
spec_comp(A, sigma2, Anoise, sigma2noise, numnc, dennc, numc, denc, thetahatfir);