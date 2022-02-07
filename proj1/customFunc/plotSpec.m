function plotSpec(num, den,Ts, sigma2hat, legendName)
% plot the spectrum given the numerator (polynomial coeff in z^(-1)
% denom, sampling rate
w=linspace(0,pi);
figure
[mags,phases,ws]=dbode(num,den,Ts,w);
semilogy(ws,mags.^2*sigma2hat);
legend(legendName);
title('Spectra')
xlabel('Frequency (rad/s)')
ylabel('Magnitude')
