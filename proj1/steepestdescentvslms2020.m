% Illustration of steepest descent and LMS
%  
% 
%% Magnus Jansson
%%
close all,clear all, clc

N =  5000; 
a1 = [.9 -.3];
roots([1 -a1])

e = randn(N,1);
y =  filter(1,[1 -a1],e);

r=ar2cov([1 -a1],1,2); 
syy1 = toeplitz(r(1:2))
syx1 = r(2:3)

%N =  2000; 
a2 = [.9 -.5];
roots([1 -a2])
thtrue = [a1'*ones(1,N) a2'*ones(1,N)];

e = randn(N,1);
y = [y ;filter(1,[1 -a2],e)]; N=2*N;
figure(1), plot(y), title('Data series')
pause

figure(2)

r = ar2cov([1 -a2],1,2); 
syy2 = toeplitz(r(1:2))
syx2 = r(2:3)

lam1 = eig(syy1)
lam2 = eig(syy2)

switch 1
    case 1
        fs = .005
        fl = .005
    case 2
        fs = .05
        fl = .05
        
    case 3
        fs = 1      % marginally stable
        fl = .05
    case 4
        fs = 1.0001 % unstable
        fl = .2
        
    case 5
        fs = max(lam1)/sum(lam1), %fastest possible convergence for the first part of data
        fl = .05
end
 
mus1 = fs*2/max(lam1)
mul1 = fl*2/max(lam1)
mus2 = fs*2/max(lam2)
mul2 = fl*2/max(lam2)


ths = zeros(2,N);
thl = zeros(2,N);

plotinterval = 20;
for k = 3:N/2
    ths(:,k) = ths(:,k-1) - mus1*(syy1*ths(:,k-1) - syx1); 
    Y = y((k-1):-1:(k-2));
    thl(:,k) = thl(:,k-1) + mul1*Y*(y(k)-Y'*thl(:,k-1)); 

    if rem(k,plotinterval) == 0
      plot(ths','--b'),hold on, plot(thl','r'), plot(thtrue','k'),
      hold off, drawnow
    end
end
for k = N/2+1:N
    ths(:,k) = ths(:,k-1) - mus2*(syy2*ths(:,k-1) - syx2); 
    Y=y((k-1):-1:(k-2));
    thl(:,k) = thl(:,k-1) + mul2*Y*(y(k)-Y'*thl(:,k-1)); 

    if rem(k,plotinterval) == 0
     plot(ths','--b'),hold on, plot(thl','r'), plot(thtrue','k'),
     hold off, drawnow
    end
end

legend([ '\mu_{SD}=' num2str(fs) '*2/\lambda_1'] , '', ...
    [ '\mu_{LMS}=' num2str(fl) '*2/\lambda_1'], '' , ...
 [ 'true '] )

pause
%% 
plotlms = 1;  % == 1 plot LMS estimates

nop = 200;
TH1 = linspace(-.0 ,1.2,nop);
TH2 = linspace(-.8 ,0.4,nop);
MSE1 = zeros(nop,nop);
MSE2 = zeros(nop,nop);
for i=1:nop
for j=1:nop
 TH = [TH1(i); TH2(j)];
 MSE1(i,j) = TH'*syy1*TH -2*TH'*syx1 + syy1(1,1);
 MSE2(i,j) = TH'*syy2*TH -2*TH'*syx2 + syy2(1,1);
end
end

figure(3)

contour(TH1,TH2,MSE1',30,'b')
axis('square')
xlabel('\theta_1'),ylabel('\theta_2')
title('Contour plot of MSE cost function')

for k = 1:11:2000
    hold on
    plot(ths(1,k),ths(2,k),'x'),
    if plotlms==1, plot(thl(1,k),thl(2,k),'o'), end
    if k<100, text(ths(1,k),ths(2,k),int2str(k)), end
    drawnow
end
hold off, 

pause
hold on
contour(TH1,TH2,MSE2',30,'--g')
pause
for k = N/2+1 :11: N/2+1000
    hold on
    plot(ths(1,k),ths(2,k),'x'),
    if plotlms==1, plot(thl(1,k),thl(2,k),'o'), end 
    drawnow
   %text(ths(1,k),ths(2,k),int2str(k))
end
hold off

