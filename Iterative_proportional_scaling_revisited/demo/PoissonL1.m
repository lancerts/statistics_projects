function [betaEst, totNum, lossval, fval] = PoissonL1(X, y, q, Lambda, errBnd, maxIT, betaInit)
% This is the main L1-penalized Poisson procedure. 
% Arguements:   X -- BINARY model matrix (=0/1) which may contain the intercept
%               y -- Poisson responses (>=0), 
%               q -- offset (>0),
%               Lambda -- the lambda parameter vector (same size as beta or a scalar), 
%               errBnd -- error tolerance, 
%               maxIT -- max iteration number allowed,
%               betaInit -- initial estimate, 
%%%%%
% Values: betaEst gives the coefficient estimates The
%           number of iterations performed is recorded in totNum.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We do not perform column-scaling on X 

if any(y < 0) || any(q <= 0) || any(any(X~= round(X))) || any (Lambda < 0) 
    error('Wrong input!')
end
if ~isempty(find(~sum(X), 1))
    error ('There are zero columns!')
end

[n, d] = size(X);  
if isempty(betaInit)
    betaInit = zeros(d, 1); 
else
    betaInit = reshape(betaInit, d, 1);
end
if numel(Lambda) == 1
    Lambda_new = Lambda * ones(d, 1);
else
    Lambda_new = reshape(Lambda, d, 1);
end

beta_cur =  betaInit;
t = 0;
debug = 1; %1;
%%%%%%%%%%% ITERATION starts %%%%%%%%%%%%%%%%%%%%
while 1
    t = t + 1;
    gamma = beta_cur;
    for j = 1:d
        q_new =  q .* exp( X * gamma - X(:,j) * gamma(j) );
%        q_new =  q .* exp( X * gamma );
        gamma(j) = PoissonSoftTh(X(:, j), y, q_new, Lambda_new(j));
    end
    beta_new = gamma;
    
    if debug == 1
        sprintf(['iter--', num2str(t), ', Func val: ', num2str(PoissonL1Obj(X, y, q, beta_new, Lambda))])
    end
    
    if  (max(abs(beta_new-beta_cur)) < errBnd ) || t >= maxIT 
        break;
    end
    
    beta_cur = beta_new;
end

betaEst = beta_new;
totNum = t;
[fval, lossval, ~] = PoissonL1Obj(X, y, q, betaEst, Lambda, 1);
end

function [gamma] = PoissonSoftTh(x, y, q_new, lambda);
innprod_xy = sum(y .* x);
innprod_xq = sum(q_new .* x);
if innprod_xq==0
     warning('Marginal total is zero!')
end
delta = innprod_xy - innprod_xq;

if abs(delta) >= lambda
    gamma = log( (innprod_xy - lambda * sign(delta)) / innprod_xq );
else
    gamma = 0;
end
end

function [fval, lossval, penval] = PoissonL1Obj(X, y, q, beta, Lambda, pen)

if nargin == 5 | ~exist('pen', 'var') | isempty(pen)
    pen = 1;
end

eta = X * beta;

lossval = [];
penval = [];

lossval = - sum(y .* eta) + sum(q .* exp(eta));

if pen == 1 %|| pen == 2
    penval = norm(Lambda .* beta, 1);
    fval = lossval + penval;
end
if pen == 0 %|| pen == 2
    fval =  lossval;
end
end


