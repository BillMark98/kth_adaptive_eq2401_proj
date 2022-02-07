clear; clc;
% load audio
[z,fs] = audioread('EQ2401project1data2022.wav');

%% time domain
plot(z)
%% raw sound
soundsc(z,fs);

%% DFT
% y = fft(z);                               % Compute DFT of x
% m = abs(y);                               % Magnitude
% y(m<1e-6) = 0;
% p = unwrap(angle(y));  
% 
% f = (0:length(y)-1)*100/length(y);        % Frequency vector
% 
% subplot(2,1,1)
% plot(f,m)
% title('Magnitude')
% ax = gca;
% ax.XTick = [15 40 60 85];
% 
% subplot(2,1,2)
% plot(f,p*180/pi)
% title('Phase')
% ax = gca;
% ax.XTick = [15 40 60 85];   
% 
% plotDFT(z,fs,'raw signal');


plotDFT_rad(z,'raw signal');
%% get background noise and sound 

% guess
noise_sample = z(end-4000:end);
sound_sample = z(100:10000);

%% check noise and sound

soundsc(noise_sample, fs);


%% dft noise

% plotDFT(noise_sample, fs, "Noise");
figure
plotDFT_rad(noise_sample,"Noise");
%% 
soundsc(sound_sample, fs);

%% dft sound
% 
% plotDFT(sound_sample,fs, "Sound+noise");
figure
plotDFT_rad(sound_sample, "Sound+noise");
%% filters

% estimate sigmayx, yy

% FIR filter

%% estimate noise
N_noise = 4;
[Anoisehat, sigma2noisehat] = ar_id(noise_sample,N_noise);
w=linspace(0,pi);
figure
[mags,phases,ws]=dbode(1,Anoisehat,1,w);
semilogy(ws,mags.^2*sigma2noisehat);
legend('Noise');
title('Spectra')
xlabel('Frequency (rad/s)')
ylabel('Magnitude')

% % tf function uses positive s power so flip the order
% tf_num = zeros(1,length(Anoisehat)-1);
% tf_num(1) = 1;
% bode(tf(tf_num, fliplr(Anoisehat)));  

%% first filter out to get sound
omegaLow = 0.06;
omegaHigh = 1.2;
cutoutN = 200;

% omegaLow = 0.2;
% omegaHigh = 0.8;
% cutoutN = 3;


filter_sound = moved_idealLP(sound_sample, omegaLow, omegaHigh, cutoutN);
filter_sound = columnVector(filter_sound);
figure
plotDFT_rad(filter_sound, 'filtered signal');

% estimate 1st, nearly zero
m_sound = mean(sound_sample);
m_filter_sound = mean(filter_sound);

%% filter sound 
soundsc(filter_sound,fs);

%% FIR
fir_length = 1000;
param_length = fir_length + 1; % add one for dc estimate
x = filter_sound;
y = sound_sample;
SigmaYxhat = xcovhat(x,y, fir_length);
SigmaYYhat = covhat(y,fir_length);
% SigmaYeYehat = [SigmaYYhat, m_sound * ones(fir_length,1); m_sound * ones(1,fir_length), 1];

% xhatfir = [SigmaYxhat', m_filter_sound] * (SigmaYeYehat \ y);
[xhatfir, theta] = firw(y, SigmaYxhat, SigmaYYhat);
plotDFT_rad(xhatfir, ['FIR(' int2str(fir_length) ') Filter '])

%% sound fir
soundsc(xhatfir, fs);

%% FIR alternative

%% FIR
fir_length = 1000;
param_length = fir_length + 1; % add one for dc estimate
x = filter_sound;
y = sound_sample;
SigmaYyhat = xcovhat(y,y, fir_length);
SigmaVvhat = xcovhat(noise_sample, noise_sample,fir_length);
SigmaYxhat = SigmaYyhat - SigmaVvhat;

SigmaYYhat = covhat(y,fir_length);
% SigmaYeYehat = [SigmaYYhat, m_sound * ones(fir_length,1); m_sound * ones(1,fir_length), 1];

% xhatfir = [SigmaYxhat', m_filter_sound] * (SigmaYeYehat \ y);
[xhatfir, theta] = firw(y, SigmaYxhat, SigmaYYhat);
plotDFT_rad(xhatfir, ['FIR(' int2str(fir_length) ') Filter '])


%% sound fir
soundsc(xhatfir, fs);

%% estimate sound
N_sound = 8;
[Ahat, sigma2hat] = ar_id(filter_sound,N_sound);
figure
plotSpec(1,Ahat,1, sigma2hat, 'estimate Sound');

% [SigmaYYhat, SigmaYxhat] = firw_cov_add(Ahat,sigma2hat,Anoisehat, sigma2noisehat,N_sound);

% [xhatfir, thetahatfir] = firw(z, SigmaYxhat, SigmaYYhat);


% % alternative
% [MySigmaYY, MySigmaYx] = firw_cov_add(A, sigma2, Anoise, sigma2noisehat,N);

%% NCW
[PhixyNum,PhixyDen,PhiyyNum,PhiyyDen] = ...
               spec_add(Ahat, sigma2hat, Anoisehat, sigma2noisehat);
[xhatnc, numnc, dennc] = ncw(z, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen);

%% plot NC
figure
plotDFT_rad(xhatnc,'Non-causal')

% psd method

%% phixy / phiyy
[phixyEst, w_xy] = cpsd(filter_sound, sound_sample);
[phiyyEst, w_yy] = pwelch(sound_sample);

H_est = phixyEst ./ phiyyEst;

%% sound
soundsc(xhatnc, fs);
%% default m = 0
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

%% other approaces
omegaLow = 0.06;
omegaHigh = 1.2;
cutoutN = 200;

% omegaLow = 0.2;
% omegaHigh = 0.8;
% cutoutN = 3;


filter_z = moved_idealLP(z, omegaLow, omegaHigh, cutoutN);
figure
plotDFT_rad(filter_z, 'filtered signal');

%% sound
soundsc(filter_z, fs);