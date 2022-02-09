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

figure
plotDFT_rad(z,'raw signal');
%% get background noise and sound 

% guess
noise_sample = z(end-4000:end);
sound_sample = z(100:end-4001);

plot(noise_sample)

%% check noise and sound

soundsc(noise_sample, fs);
%audiowrite("noise_sample.wav",noise_sample, fs);
%% xcorr2 
plot(xcorr2(noise_sample))

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
audiowrite("sound_noise.wav",sound_sample, fs);
%% filters

% estimate sigmayx, yy

% FIR filter

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
% % tf function uses positive s power so flip the order
% tf_num = zeros(1,length(Anoisehat)-1);
% tf_num(1) = 1;
% bode(tf(tf_num, fliplr(Anoisehat)));  
noise_diff = noise_sample - noise_gen;
figure
plot(noise_diff);

mean(noise_diff)
var(noise_diff)
%%
mean(noise_sample.^2)

%%
mean(noise_diff)
var(noise_diff)
%% test noise generated
plotDFT_rad(noise_gen, 'noise gen')
soundsc(noise_gen, fs);

%% first filter out to get sound
omegaLow = 0.1;
omegaHigh = 0.8;
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

soundsc(filter_sound,fs);
%audiowrite("filter_sound.wav",filter_sound, fs);

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
audiowrite("xhatfir.wav",xhatfir, fs);

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

%% sound nc
soundsc(xhatnc,fs);
% psd method
audiowrite("xhatnc.wav",xhatnc, fs);
%% alt 1 phixy / phiyy
[phixyEst, w_xy] = cpsd(filter_sound, sound_sample, [],[], floor(length(filter_sound)/2));

[phiyyEst, w_yy] = pwelch(sound_sample, [],[], floor(length(sound_sample)/2));

% extend phixy and phiyy
% phixy is complex conjugate flip
phixyEst = [phixyEst;conj(phixyEst(end:-1:1))];
phiyyEst = [phiyyEst;phiyyEst(end:-1:1)];
H_est = phixyEst ./ phiyyEst;
z_trans = fft(z, length(H_est));
x_trans = H_est .* z_trans;
x = ifft(x_trans, length(z));
x_conv = abs(x);
% 
% h_est = ifft(H_est);
% 
% x_conv = conv(abs(h_est), z); % does not work, h_est complex, meaning apprx in frequency domain not
%%
plotDFT_rad(x_conv,'xconv');
soundsc(x_conv, fs)

%%
z_trans = fftshift(fft(z)/length(z));

plot(abs(z_trans))
% 1/2*(length(z_trans)-1)
onesided=z_trans((1/2*(length(z_trans))-1):end)

plot(abs(onesided))
xfft=onesided(1:4951).*H_est
Hnew=[H_est(end:-1:2);H_est]

%% alt 2  phixy = phiyy - phivv
[phivvEst, wvv] = pwelch(noise_sample,[],[], length(filter_sound));

H_est2 = 1 - phivvEst ./ phiyyEst;
h_est2 = ifft(H_est2);
%% sound of filtered
soundsc(x_conv,fs);

%% default m = 0
m = 0;
[xhatc, numc, denc] = cw(z, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen, m);

%% sound
soundsc(xhatc, fs);
audiowrite("xhatc.wav",xhatc, fs);
figure
plotDFT_rad(xhatc, 'Causal')
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


%% kalman forward

%% 1 step estimate sound
N_sound = 8;
[Ahat_rcos, sigma2hat_rcos] = ar_id(filter_sound_rcos,N_sound);
figure
plotSpec(1,Ahat_rcos,1, sigma2hat_rcos, 'estimate Sound');

N_noise = 4;
[Anoisehat, sigma2noisehat] = ar_id(noise_sample,N_noise);
% w=linspace(0,pi);
% figure
% [mags,phases,ws]=dbode(1,Anoisehat,1,w);
% semilogy(ws,mags.^2*sigma2noisehat);
% legend('Noise');
% title('Spectra')
% xlabel('Frequency (rad/s)')
% ylabel('Magnitude')

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
[yhat,xhatfilt,xhatpred,P,Q] = kalman(y, F, G, H, R1, R2, x0, Q0);

sound_pred = xhatpred(:,1);

sound_diff = y - sound_pred;
%%
plotDFT_rad(sound_pred, 'kalman forward')
soundsc(sound_pred, fs);

%%
plotDFT_rad(sound_diff, 'kalman forward diff')
soundsc(sound_diff,fs);


%% kalman backward
%% step 1
% first time reversal
noise_sample_rev = flip(noise_sample);
sound_sample_rev = flip(sound_sample);
figure
plotDFT_rad(noise_sample_rev, 'noise reversed')
figure
plotDFT_rad(sound_sample_rev, 'sound reversed')

soundsc(noise_sample_rev, fs);
soundsc(sound_sample_rev, fs);

%% step 2
% estimate the AR model param ov v
N_order = 4;
[Anoisehat_rev, sigma2noisehat_rev] = ar_id(noise_sample_rev,N_order);
w=linspace(0,pi);
figure
[mags,phases,ws]=dbode(1,Anoisehat_rev,1,w);
semilogy(ws,mags.^2*sigma2noisehat_rev);
legend('Noise');
title('Spectra')
xlabel('Frequency (rad/s)')
ylabel('Magnitude')

%% step 3 constrcut kalman state space model 
F = diag(ones(N_order-1,1),-1);
F(1,:) = -Anoisehat_rev(2:end);
G = zeros(N_order,1);
G(1) = 1;
H = zeros(1,N_order);
H(1) = 1;
R1 = sigma2noisehat_rev;
% first alternative, fixed R2
% R2 = E(signal^2) = E(total^2) - E(noise^2)
R2 = mean(sum(sound_sample_rev.^2)) - sigma2noisehat_rev;

% initial vector 
% x0 = [sound_sample_rev(1); noise_sample_rev(end:-1:end-N_order+2)];
x0 = [noise_sample_rev(end:-1:end-N_order+1)];
y = sound_sample_rev;
Q0 = zeros(N_order);
Q0(1,1) = (sound_sample_rev(1) - F(1,:) * x0)^2 + sigma2noisehat_rev;

% [yhat,xhatfilt,xhatpred,P,Q] = kalman(y, F, G, H, R1, R2, x0, Q0);

[yhat,xhatfilt,xhatpred,P,Q] = kalman_varyR(y, F, G, H, R1, R2, x0, Q0);

yfinal = flip(y - xhatpred(:,1));
%%
soundsc(flip(yfinal),fs);
%%
soundsc(flip(xhatfilt(:,1)),fs)

%% 
soundsc(flip(y),fs);