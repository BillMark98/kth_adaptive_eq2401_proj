% Illustration of steepest descent, LMS and RLS
%  Prediction of AR(2) process where params change
%
%  added RLS 210217 , modifications 220201
% 
%% Magnus Jansson
%%
%close all,clear all, clc

plotinterval = 1;%odd number may be good in case of marginal stab.
plotlms  = 0; %=1 plot LMS 
plotrls = 0;  %=1 plot RLS 

% First model
N =  120; 
a1 = [.9 -.3];
%roots([1 -a1])

e = randn(N,1);
y =  filter(1,[1 -a1],e);

r=ar2cov([1 -a1],1,2); 
syy1 = toeplitz(r(1:2));
syx1 = r(2:3);

% Second model
a2 = [.9 -.5];
%roots([1 -a2])
thtrue = [a1'*ones(1,N) a2'*ones(1,N)];

e = randn(N,1);
y = [y ;filter(1,[1 -a2],e)]; 
N = 2*N;

r = ar2cov([1 -a2],1,2); 
syy2 = toeplitz(r(1:2));
syx2 = r(2:3);


figure(1), plot(y), title('Data series')
%pause



lam1 = eig(syy1);
lam2 = eig(syy2);

switch 5
    case 1
        fs = .001
        fl = .001
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

% RLS parameters
lambda = .999;
P = .05*eye(2);

ths = zeros(2,N);
thl = zeros(2,N);
thr = zeros(2,N);


for k = 3:N/2
    % Steepest descent
    ths(:,k) = ths(:,k-1) - mus1*(syy1*ths(:,k-1) - syx1); 
    
    % LMS
    Y = y((k-1):-1:(k-2));
    thl(:,k) = thl(:,k-1) + mul1*Y*(y(k)-Y'*thl(:,k-1));   

% RLS
  % Update K
  K=P*Y/(lambda+Y'*P*Y);  
  % Update P
  P=1/lambda*(P-K*Y'*P);
  thr(:,k)=thr(:,k-1)+K*(y(k)-Y'*thr(:,k-1));

  if rem(k,plotinterval) == 0
      figure(2), 
      plot(ths','--b','DisplayName',['\mu_{SD}=' num2str(fs) '*2/\lambda_1']),
      hold on,
      if plotlms, plot(thl','r','DisplayName',['\mu_{LMS}=' num2str(fl) '*2/\lambda_1'] ), end 
      if plotrls, plot(thr','g','DisplayName',[ '\lambda_{RLS}=' num2str(lambda) ] ), end
      plot(thtrue','k','DisplayName','true '),
      hold off,
      xlabel('n'),ylabel('$\hat \theta$','Interpreter','latex')
      drawnow
  end
end
% pause
for k = N/2+1:N
    %Steepest descent
    ths(:,k) = ths(:,k-1) - mus2*(syy2*ths(:,k-1) - syx2); 
    % LMS
    Y=y((k-1):-1:(k-2));
    thl(:,k) = thl(:,k-1) + mul2*Y*(y(k)-Y'*thl(:,k-1)); 
% RLS
  K=P*Y/(lambda+Y'*P*Y);  
  P=1/lambda*(P-K*Y'*P);
  thr(:,k)=thr(:,k-1)+K*(y(k)-Y'*thr(:,k-1));
%  
  if rem(k,plotinterval) == 0
      plot(ths','--b','DisplayName',['\mu_{SD}=' num2str(fs) '*2/\lambda_1']),
      hold on,
      if plotlms, plot(thl','r','DisplayName',['\mu_{LMS}=' num2str(fl) '*2/\lambda_1'] ), end 
      if plotrls, plot(thr','g','DisplayName',[ '\lambda_{RLS}=' num2str(lambda) ] ), end
      plot(thtrue','k','DisplayName','true '),
      hold off, 
      xlabel('n'),ylabel('$\hat \theta$','Interpreter','latex')
      drawnow
  end
end

legend
%legend([ '\mu_{SD}=' num2str(fs) '*2/\lambda_1'] , '', ...
%    [ '\mu_{LMS}=' num2str(fl) '*2/\lambda_1'], '' , ...
%     [ '\lambda_{RLS}=' num2str(lambda) ], '' , ...
% [ 'true '] )

set(gca,'fontsize',15)
return %pause


%% Illustration of iterates and the MSE criteria
% 
%plotlms = 0;  % == 1 plot LMS estimates
plotinterval = 1;

nop = 200;
TH1 = linspace(-.0 , 1.3,nop);
TH2 = linspace(-.9 , 0.3,nop);
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
title('Contour plot of MSE criterion function')

for k = 1:plotinterval:N/2
    hold on
    plot(ths(1,k),ths(2,k),'x','MarkerSize',12),
    if plotlms == 1, plot(thl(1,k),thl(2,k),'o'), end
    if k < 400, text(ths(1,k),ths(2,k),int2str(k)), end
    drawnow
end
hold off, 

pause
hold on
contour(TH1,TH2,MSE2',30,'--g','LineWidth',2)
pause
for k = N/2+1 :plotinterval: N
    hold on
    plot(ths(1,k),ths(2,k),'x','MarkerSize',12),
    if plotlms==1, plot(thl(1,k),thl(2,k),'o'), end 
    drawnow
   text(ths(1,k),ths(2,k),int2str(k))
end
hold off

