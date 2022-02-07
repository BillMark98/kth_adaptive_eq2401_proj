function plotDFT_rad(z,plotTitle)
% given time domain sampled z, express in sampled omega
% plot the dft
% fs  the sampling frequency
y = fft(z);                               % Compute DFT of x

omega_samp = (0:length(y)-1)*2 * pi/length(y);        % Frequency vector
% only need to plot half
y = y(1:floor(length(y)/2));
omega_samp = omega_samp(1:length(y));
m = abs(y);                               % Magnitude
y(m<1e-6) = 0;
p = unwrap(angle(y));  

subplot(2,1,1)
plot(omega_samp,m)
title('Magnitude')
% ax = gca;
% ax.XTick = [15 40 60 85];
xlabel('Hz');
subplot(2,1,2)
plot(omega_samp,p*180/pi)
title('Phase')
% ax = gca;
% ax.XTick = [15 40 60 85];   
xlabel('Hz');

sgtitle(plotTitle)