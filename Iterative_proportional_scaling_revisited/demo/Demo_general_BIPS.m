%opts           Stuct contatins  the following options
%TOLGRAD        Tolerance for termination
%MAX_ITERS      Maxmimum iterations allowed
%MAX_TIME       Maxmimum time allowed 
%BETA_INIT      p*1 column vector as the initial value for the
%               beta estimation (default value: all components are 0)
%Method         Newton or LBFGS

%% Parameter Specification
clear,clc;
addpath(genpath(pwd))
rng(101);
N = 2000;
p = 1000;
rho = 0.8;    
q = ones(N,1);  %offset
gk = 200;    %block size



%Set options
MAX_ITERS = 1e+2;
MAX_TIME = 20;
TOLGRAD = 1e-6;
opts=struct('MAX_ITERS',MAX_ITERS,'MAX_TIME',MAX_TIME,'TOLGRAD',TOLGRAD, 'Method','Newton');
%'Method','lbfgs');



%% Design matrix generation
corrMat = rho .^ abs((1:p)'*ones(1,p)-ones(p,1)*(1:p));
mu = zeros(N,p);
X = mvnrnd(mu,corrMat); 
fprintf('Design matrix obtained\n')
X = X/(100*max(max(abs(X))));

%% Beta_true generation
mu = [-10;10];
sigma = cat(3,[1],[1]);
qq = ones(1,2)/2;
obj = gmdistribution(mu,sigma,qq);
beta_true = random(obj,p+1);
beta_true(1) = 10;

% Observed counts generation
n = poissrnd(q.*exp([ones(N,1),X]*beta_true));


%% Perform BIPS
[beta,~,~,gradval,iter,time] = B_IPS(X,n,q,gk,opts,beta_true);

iter
time








