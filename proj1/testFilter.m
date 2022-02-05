%% 3.7
clc; clear;
K = 10000;
sigma2 = 1;


[x,y] = pol2cart(pi/4, 0.98);
A = poly([x + 1i * y, x - 1i * y]);

sigma2noise = 1;
[x_n, y_n] = pol2cart(pi/3, 0.98);
Anoise = poly([x_n + 1i * y_n, x_n - 1i * y_n]);
N = 20;
filterCompare(N,A,sigma2,Anoise,sigma2noise,K,30);


%% 3.8
clc; clear;
K = 10000;
sigma2 = 1;
e = sqrt(sigma2) * randn(K,1);

[x,y] = pol2cart(pi/4, 0.98);
A = poly([x + 1i * y, x - 1i * y]);

sigma2noise = 1;
[x_n, y_n] = pol2cart(pi/3, 0.98);
Anoise = [1];

% order of FIR filter
N = 10;
filterCompare(N,A,sigma2,Anoise,sigma2noise,K,30);

