clear; clc;
% load audio
[z,fs] = audioread('EQ2401project1data2022.wav');

%% time domain
plot(z)
%% raw sound
soundsc(z,fs);

%% DFT
y = fft(z);                               % Compute DFT of x
m = abs(y);                               % Magnitude
y(m<1e-6) = 0;
p = unwrap(angle(y));  

f = (0:length(y)-1)*100/length(y);        % Frequency vector

subplot(2,1,1)
plot(f,m)
title('Magnitude')
ax = gca;
ax.XTick = [15 40 60 85];

subplot(2,1,2)
plot(f,p*180/pi)
title('Phase')
ax = gca;
ax.XTick = [15 40 60 85];   

%% get background noise and sound 

% guess
noise_sample = z(end-3000:end);
sound_sample = z(100:10000);
%% filters

% estimate sigmayx, yy

% FIR filter

% estimate noise
N_noise = 5;
[Anoisehat, sigma2noisehat] = ar_id(noise_sample,N_noise);

N_sound = 8;
[Ahat, sigma2hat] = ar_id(sound_sample,N_sound);

[SigmaYYhat, SigmaYxhat] = firw_cov_add(Ahat,sigma2hat,Anoisehat, sigma2noisehat,N_sound);

[xhatfir, thetahatfir] = firw(z, SigmaYxhat, SigmaYYhat);


% % alternative
% [MySigmaYY, MySigmaYx] = firw_cov_add(A, sigma2, Anoise, sigma2noisehat,N);

% NCW
[PhixyNum,PhixyDen,PhiyyNum,PhiyyDen] = ...
               spec_add(Ahat, sigma2hat, Anoisehat, sigma2noisehat);
[xhatnc, numnc, dennc] = ncw(z, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen);

% default m = 0
m = 0;
[xhatc, numc, denc] = cw(z, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen, m);

%% DFTs
figure
plotDFT(xhatfir, 'FIR')

figure
plotDFT(xhatnc, 'Non-causal')

figure
plotDFT(xhatc, 'Causal')

% filtered sound

%% fir
soundsc(xhatfir,fs);
%% nc
soundsc(xhatnc, fs);
%% causal
soundsc(xhatc, fs);