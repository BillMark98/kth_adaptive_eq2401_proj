clear; clc;
% add path
addpath(genpath(pwd))
% load audio
[z,fs] = audioread('EQ2401project1data2022.wav');

%% time domain
plot(z)
%% raw sound
soundsc(z,fs);

% plot dft of sound
figure
plotDFT_rad(z,'raw signal');
%% get background noise and sound 

% guess
noise_sample = z(end-4000:end);
sound_sample = z(100:end-4001);

%% dft noise

figure
plotDFT_rad(noise_sample,"Noise");
soundsc(noise_sample, fs);
%% dft sound

figure
plotDFT_rad(sound_sample, "Sound+noise");
soundsc(sound_sample,fs);
% estimate 1st, nearly zero
m_sound = mean(sound_sample)
%% FIR
fir_length = 1000;

y = sound_sample;
SigmaYyhat = xcovhat(y,y, fir_length);
SigmaVvhat = xcovhat(noise_sample, noise_sample,fir_length);
SigmaYxhat = SigmaYyhat - SigmaVvhat;
SigmaYYhat = covhat(y,fir_length);
[xhatfir, thetahatfir] = firw(y, SigmaYxhat, SigmaYYhat);
plotDFT_rad(xhatfir, ['FIR(' int2str(fir_length) ') Filter '])

%% sound fir
soundsc(xhatfir, fs);

% noncausal and causal
%% estimate noise
N_noise = 6;
[Anoisehat, sigma2noisehat] = ar_id(noise_sample,N_noise);
w=linspace(0,pi);
figure
[mags,phases,ws]=dbode(1,Anoisehat,1,w);
semilogy(ws,mags.^2*sigma2noisehat);
legend('Noise');
title('Spectra')
xlabel('Frequency (rad/s)')
ylabel('Magnitude')

% generate the noise
e_n = rand(length(noise_sample),1) * sqrt(sigma2noisehat);
noise_gen = filter(1, Anoisehat, e_n);
% plot noise difference
noise_diff = noise_sample - noise_gen;
figure
plot(noise_diff);

mean(noise_diff)
var(noise_diff)
%% test noise generated
plotDFT_rad(noise_gen, 'noise gen')
soundsc(noise_gen, fs);

%% get signal x(n) = y(n) - v_hat(n)
e_n = rand(length(sound_sample),1) * sqrt(sigma2noisehat);
vhat = filter(1,Anoisehat,e_n);
xhat = sound_sample - vhat;
soundsc(xhat,fs);
%% estimate sound
N_sound = 8;
[Ahat, sigma2hat] = ar_id(xhat,N_sound);
plotSpec(1,Ahat,1, sigma2hat, 'estimate Sound');

%% NCW
[PhixyNum,PhixyDen,PhiyyNum,PhiyyDen] = ...
               spec_add(Ahat, sigma2hat, Anoisehat, sigma2noisehat);
[xhatnc, numnc, dennc] = ncw(z, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen);

%% plot NC
figure
plotDFT_rad(xhatnc,'Non-causal')

%% sound nc
soundsc(xhatnc,fs);

%% CW
m = 0;
[xhatc, numc, denc] = cw(z, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen, m);
figure
plotDFT_rad(xhatc,'Causal')
%% sound cw
soundsc(xhatc,fs);

%% compare
spec_comp(Ahat,sigma2hat,Anoisehat,sigma2noisehat,numnc,dennc,numc,denc,thetahatfir);

%% alternatives


%% first filter out to get sound
omegaLow = 0.1;
omegaHigh = 0.8;
cutoutN = 200;


filter_sound = moved_idealLP(sound_sample, omegaLow, omegaHigh, cutoutN);
filter_sound = columnVector(filter_sound);
figure
plotDFT_rad(filter_sound, 'filtered signal');

% estimate 1st, nearly zero
m_sound = mean(sound_sample);
m_filter_sound = mean(filter_sound);

soundsc(filter_sound,fs);

%% raised cosine
omegaLow = 0.01;
omegaHigh = 1.002;
cutoutN = 400;
beta = 0.25;
filter_sound_rcos = raisedCosShiftFilter(sound_sample,omegaLow, omegaHigh, beta, cutoutN);
filter_sound_rcos = columnVector(filter_sound_rcos);
figure
plotDFT_rad(filter_sound_rcos, 'filtered signal');

soundsc(filter_sound_rcos,fs);

%% estimate sound, use raised cosine for example
N_sound = 8;
[Ahat, sigma2hat] = ar_id(filter_sound_rcos,N_sound);
plotSpec(1,Ahat,1, sigma2hat, 'estimate Sound');

%% NCW
[PhixyNum,PhixyDen,PhiyyNum,PhiyyDen] = ...
               spec_add(Ahat, sigma2hat, Anoisehat, sigma2noisehat);
[xhatnc, numnc, dennc] = ncw(z, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen);

% plot NC
figure
plotDFT_rad(xhatnc,'Non-causal')

% CW
m = 0;
[xhatc, numc, denc] = cw(z, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen, m);
figure
plotDFT_rad(xhatc,'Causal')
%% sound cw
soundsc(xhatc,fs);

%% sound nc
soundsc(xhatnc,fs);

%% compare
spec_comp(Ahat,sigma2hat,Anoisehat,sigma2noisehat,numnc,dennc,numc,denc,thetahatfir);


%% Use estimate of psd to calculate nc
%% alt  phixy / phiyy

[phixyEst, w_xy] = cpsd(filter_sound, sound_sample, [],[], length(sound_sample));

[phiyyEst, w_yy] = pwelch(sound_sample, [],[], length(sound_sample));

H_est = phixyEst ./ phiyyEst;   
N_len = 2 * length(H_est) - 1;

z_trans = fft(z, N_len);
z_trans_onesided = z_trans(1:length(H_est));
x_trans = conjugateFlip(H_est .* z_trans_onesided);
x = ifft(x_trans);
% weird sound
plotDFT_rad(x,'xconv');
soundsc(x, fs)


%% conv in time domain
H_est_freqShift = H_est .* (-1).^(0:length(H_est)-1)';
H_est_spec = conjugateFlip(H_est_freqShift) * length(H_est_freqShift);
h_est = columnVector(ifft(H_est_spec));
x_nc = imfilter(sound_sample, h_est);
plotDFT_rad(x_nc,'non-causal')
soundsc(x_nc,fs);
%% phixy = phiyy - phivv
[phivvEst, wvv] = pwelch(noise_sample,[],[], length(sound_sample));

H_est2 = 1 - phivvEst ./ phiyyEst;
% shift it by pi
H_est2_freqShift = H_est2 .* (-1).^(0:length(H_est2)-1)';
% conjugate flip to get 0 to 2pi spectrum
H_est2_spec = conjugateFlip(H_est2_freqShift);
% rescale to correct value for ifft
H_est2_spec = H_est2_spec * length(H_est2_spec);
h_est2 = columnVector(ifft(H_est2_spec));
x_nc = imfilter(sound_sample, h_est2);
plotDFT_rad(x_nc,'non-causal')
soundsc(x_nc,fs);

%% kalman forward

%% 1 step estimate sound
N_sound = 12;
[Ahat_rcos, sigma2hat_rcos] = ar_id(filter_sound_rcos,N_sound);

plotSpec(1,Ahat_rcos,1, sigma2hat_rcos, 'estimate Sound');

N_noise = 6;
[Anoisehat, sigma2noisehat] = ar_id(noise_sample,N_noise);

%% 2 step set state space model
F = diag(ones(N_sound-1,1),-1);
F(1,:) = -Ahat_rcos(2:end);
G = zeros(N_sound,1);
G(1) = 1;
H = zeros(1,N_sound);
H(1) = 1;
R1 = sigma2hat_rcos;
% first alternative, fixed R2
% R2 = E(signal^2) = E(total^2) - E(noise^2)
R2 = sigma2noisehat;

x0 = zeros(N_sound,1);
Q0 = zeros(N_sound);
y = sound_sample;
[yhat_k,xhatfilt_k,xhatpred_k,P_k,Q_k] = kalman(y, F, G, H, R1, R2, x0, Q0);

updateStep = 8;
updateFactor = 0.4;
[yhat_kad,xhatfilt_kad,xhatpred_kad,P,Q] = kalman_adapt(y, F, G, H, R1, R2, x0, Q0,updateFactor, updateStep);

sound_pred_k = xhatfilt_k(:,1);

sound_pred_kad = xhatfilt_kad(:,1);

%% kalman
plotDFT_rad(sound_pred_k, 'kalman')
soundsc(sound_pred_k, fs);

%% kalman adaptive
plotDFT_rad(sound_pred_kad, 'kalman')
soundsc(sound_pred_kad, fs);
