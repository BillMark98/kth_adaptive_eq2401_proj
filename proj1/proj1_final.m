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
m_sound = mean(sound_sample);
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
figure
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
