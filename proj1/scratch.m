Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector
X = sin(2*pi*120*t);
Y = fft(X);

Y_mag = abs(Y);
f_bins = (0:L-1)/L * Fs;
plot(f_bins, Y_mag);
