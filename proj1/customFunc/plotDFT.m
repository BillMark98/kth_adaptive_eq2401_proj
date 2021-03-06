function plotDFT(z,fs, plotTitle)
% given time domain sampled z,
% plot the dft
% fs  the sampling frequency
y = fft(z);                               % Compute DFT of x

f = (0:length(y)-1)*fs/length(y);        % Frequency vector
% only need to plot half
y = y(1:floor(length(y)/2));
f = f(1:length(y));
m = abs(y);                               % Magnitude
y(m<1e-6) = 0;
p = unwrap(angle(y));  

subplot(2,1,1)
plot(f,m)
title('Magnitude')
% ax = gca;
% ax.XTick = [15 40 60 85];
xlabel('Hz');
subplot(2,1,2)
plot(f,p*180/pi)
title('Phase')
% ax = gca;
% ax.XTick = [15 40 60 85];   
xlabel('Hz');

sgtitle(plotTitle)