Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector
X = sin(2*pi*120*t) + 1;
Y = fft(X) / L;

Y_mag = abs(Y);
f_bins = (0:L-1)/L * Fs;
figure
plot(f_bins, Y_mag);

x_retrans = ifft(Y * L);
figure
plot(t, X);
hold on
plot(t, x_retrans);
legend('x_orig','x_retrans');
x_diff = x_retrans - X;
max(abs(x_diff))

%% fft shift


Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector
X = zeros(1,L);
X(L/2) = 1;
Y = fft(X) / L;

Y_mag = abs(Y);
f_bins = (0:L-1)/L * Fs;
figure
plot(f_bins, Y_mag);

x_retrans = ifft(Y * L);
figure
plot(t, X);
hold on
plot(t, x_retrans);
legend('x_orig','x_retrans');
x_diff = x_retrans - X;
max(abs(x_diff))


%% fft self imple, shift
clear;

Fs = 100;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector
f0 = 1;
X = sin(2*pi*f0 * t);
X(L/2) = 1;
Y = fft(X) / L;
Y = Y .* exp(-j * pi*Fs/(2*f0)).^((0:L-1)/L);

Z = sin(2 * pi * f0 /Fs * ((0:L-1) + Fs / (4 * f0)));

Y_mag = abs(Y);
f_bins = (0:L-1)/L * Fs;
figure
plot(f_bins, Y_mag);

x_retrans = ifft(Y * L);
figure
plot(t, X);
% hold on
figure
plot(t, x_retrans);
legend("x_retrans");

figure
plot(t, Z);
legend('Z')
% legend('x_orig','x_retrans');
% x_diff = x_retrans - X;
% max(abs(x_diff))

%% ifft shift

V = [1 2 3 4 5 6 7];
X = fft(V);
A = fftshift(X);

B = ifftshift(A);
ifft(B)
% Y = ifft(ifftshift(X))

%% cpsd

% openExample('signal/CPSDOfColoredNoiseSignalsExample')
r = randn(16384,1);
a = xcorr2(r);
plot((-16383):16383,a)
%% 
hx = fir1(30,0.2,rectwin(31));
x = filter(hx,1,r);

hy = ones(1,10)/sqrt(10);
y = filter(hy,1,r);

[a,b] = cpsd(x,y,triang(500),250,1024);

plot(b,abs(a))

%% ifft

ff_ = [1, 1+i, 2, 2, 1-i, 1];
tt = ifft(ff_)

%% rcosdesign

% first is rolloff factor
% second is the number of symbols
% third is the sampled freq, so 2 means that the span in frequency is -pi/2
% to pi/2
h = rcosdesign(0.2,10,4) * 10;
plotDFT_rad(h,'h')