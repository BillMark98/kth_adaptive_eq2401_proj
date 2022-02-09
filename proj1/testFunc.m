
%%  3.3

load wrefdata;
N = 5;
[MySigmaYY, MySigmaYx] = firw_cov_add(A, sigma2, Anoise, sigma2noisehat,N);
[xhat, theta] = firw(y, MySigmaYx, MySigmaYY);

plot(1:1:1000, xhat, '-b')
hold on
plot(1:1:1000, xhatfir, '-r')
legend('xhat','xhatfir')

%% 3.4

load wrefdata;
[MyPhixyNum,MyPhixyDen,MyPhiyyNum,MyPhiyyDen] = ...
               spec_add(A, sigma2, Anoise, sigma2noise)
[xhat, myNum, myDen] = ncw(y, MyPhixyNum, MyPhixyDen, MyPhiyyNum, MyPhiyyDen)


%% 3.6
clear;
load wrefdata;

[xhatnc,xhatc,xhatfir,numnc,dennc,numc,denc,thetahatfir] =...
    est_add(x, v, N, Ahat, sigma2hat, Anoisehat,...
    sigma2noisehat,SigmaYxhat, SigmaYYhat);

% calculate mse
% mse_fir = myMSE(xhatfir, x, 1/sigma2noisehat);
% mse_nc = myMSE(xhatnc, x, 1/sigma2noisehat);
% mse_c = myMSE(xhatc, x, 1/sigma2noisehat);

mse_fir = myMSE(xhatfir, x, 1/sqrt(sigma2noisehat));
mse_nc = myMSE(xhatnc, x, 1/sqrt(sigma2noisehat));
mse_c = myMSE(xhatc, x, 1/sqrt(sigma2noisehat));
% legends
fir_legend = sprintf("FIR(%.4f)", mse_fir);
nc_legend = sprintf("Non-causal(%.4f)", mse_nc);
c_legend = sprintf("Causal(%.4f)", mse_c);
% plot the first 30 samples

N_sample = 30;

indices = 1:1:N_sample;
figure;
plot(indices, y(indices), 'b');
hold on;
plot(indices, x(indices), ':k');
plot(indices, xhatnc(indices),'--.r');
plot(indices, xhatc(indices), '--g');
plot(indices, xhatfir(indices),'-m');

legend("Observations y", "Signal", nc_legend, c_legend, fir_legend);

figure
spec_comp(A, sigma2, Anoise, sigma2noise, numnc, dennc, numc, denc, thetahatfir);

%% 4.1
F = [1,2;1.5,0];
G = [1,1;0,1];
R1 = [2,0;0,3];
H = [0,1];
R2 = 1;
x0 = [1;0];
K = 4;
[x,y] = gen_ss(F,G,H,R1,R2,x0,K)

%% 4.2

clear;
load krefdata
[yhat,xhatfilt,xhatpred,P,Q] = kalman(y, F, G, H, R1, R2, x0, Q0)
diff_max_xhatfilt = max(xhatfiltref - xhatfilt);
diff_max_yhat = max(yhatref - yhat);
diff_max_xhatpred = max(xhatpredref - xhatpred);

%% 4.3
clear;
[x,y] = pol2cart(pi/4, 0.98);
A = poly([x + 1i * y, x - 1i * y]);
sigma2noise = 1;
L = 1000;
e = randn(L,1);
y = sqrt(sigma2noise) * filter(1,A,e);

x0 = [0,0]';
Q0 = zeros(2);
F = [-A(2), -A(3); 1, 0];
G = [1,0]';
H = [1, 0];
R1 = sigma2noise;
R2 = 0;
[yhat,xhatfilt,xhatpred,P,Q] = kalman(y, F, G, H, R1, R2, x0, Q0);
ykalman_error = y - yhat;

% causal
[PhixyNum, PhixyDen] = filtspec(1,A,sigma2noise);
% should shift by one, because z/A(z)A(z^-1)
PhixyDen = [0, PhixyDen];
[PhiyyNum, PhiyyDen] = filtspec(1,A,sigma2noise);
[yhatc,num,den] = cw(y,PhixyNum,PhixyDen,PhiyyNum,PhiyyDen,1);
yc_error = y - yhatc;

[yhatnc,num,den] = ncw(y,PhixyNum, PhixyDen, PhiyyNum, PhiyyDen);
ync_error = y - yhatnc;
figure
plot(ykalman_error);
hold on
plot(yc_error); 
plot(ync_error);
legend('kalman', 'causal', 'non-causal')
%% 5.1
test_r = randn(1000,1);
A = [1, -0.5];

e = filter(1,A,test_r);
N = 5;
sigmayy_hat = covhat(e,N);
sigmayy_hat(1,:)
r_should = ar2cov(A,1,N)


%% 5.3
clear;
test_r = randn(10000,1);
A = poly([0.5,0.2]);
N = length(A) - 1;
y = filter(1,A,test_r);
% calculate the should r_yy
r_yy = ar2cov(A,1,N);

r_yx = r_yy(2:end);
r_yy = r_yy(1:end-1);
sigmaYY_should = toeplitz(r_yy);
a_should = -sigmaYY_should\r_yx

% sigmayy_hat = covhat(y,N);
% sigmayy_hat(1,:)
% r_should = ar2cov(A,1,N)

[ahat, sigma2hat] = ar_id(y,N)
