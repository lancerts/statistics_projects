%The miminization via Random BCD is performed.

%Input:
%X       Design matrix
%n       Observed value
%q       Offset
%gk 	 Block size

%Output:
%beta       Estimated coefficients  
%mu         Estimated mean
%fval       Row vector records the value of the objective function
%           during the progress
%gradval    Row vector records the value of the infinity norm of the
%           gradient during the progress
%iter       Iteration number when BIPS is terminated
%time       Computational time of BIPS 


function [beta,mu,fval,gradval,iter,time] = B_IPS(X,n,q,gk,opts,beta_true)
%Get parameters, options
N=size(X,1);
p=size(X,2);
beta=zeros(p,1);
beta0=0;


%options for Newton (Quasi-Newton) subroutine
options={};
options.display='off';
%options.progTol= 1e-12;
%options.optTol= 1e-3;

if (~isempty(opts))
    if isfield(opts,'TOLGRAD');TOLGRAD = opts.TOLGRAD;end
    if isfield(opts,'MAX_ITERS');MAX_ITERS = opts.MAX_ITERS;end
    if isfield(opts,'BETA_INIT');beta = opts.BETA_INIT;end
    if isfield(opts,'MAX_TIME');MAX_TIME = opts.MAX_TIME;end
    if isfield(opts,'Method');Method=opts.Method;end
    if strcmp(Method,'newton')
        options.Method='newton';
    else
        options.Method='lbfgs';
    end
end




%Compute total number of blocks m
if mod(p,gk)==0
    m= p/gk;
else
    m= floor(p/gk)+1;
end


%Initialize
mu = q.*exp(X*beta);
t = 0; %the iteration number 

k = 1; 
fval = zeros(1,MAX_ITERS+1);
gradval = zeros(1,MAX_ITERS+1);
time = 0;
values.A = n.'*[ones(N,1),X];
values.q = q;
values.X = [ones(N,1),X];
values.B = sum(n);
values.n = n;
[yold,grad0] = f([beta0;beta],values);
grad0 = norm(grad0,Inf);
fval(1) = yold;
gradval(1) = grad0;
permutation = randperm(p);
tic

while 1
    
    startk = gk*(k-1)+1;
    if k<m
        endk = gk*k;
    else
        endk = p;
    end
    % Define the subproblem and use newton (quasi-netwon) to solve it
    index = permutation(startk:endk);
    betak_old = beta(index);
    values.Z = X(:,index);
    values.Ak = n.'*values.Z;
    values.Bk = mu.*exp(-values.Z*betak_old);
    g = @(x) f_BCD_re(x,values);
    betak_new = minFunc(g,betak_old,options);
    mu = mu.*exp(values.Z*(betak_new-betak_old));
    beta(index) = betak_new;
    
    k = k+1;
    if k>m
        t = t+1;
        k = 1;
        permutation = randperm(p);
        time = toc+time;
        beta0 = log(sum(n))-log(sum(mu));
        [ynew,grad] = f([beta0;beta],values);
        fval(t+1) = ynew;
        gradval(t+1) = norm(grad,Inf);
        
        %when each cycle completes, check relative function value change
        if norm(grad,Inf)/max(grad0,1)>TOLGRAD
            if mod(t,1) == 0
                beta0 = log(sum(n))-log(sum(mu));
                error_beta = norm([beta0;beta]-beta_true,2)^2/norm(beta_true,2)^2;
                fprintf('B_IPS2 iteration= %.i, time= %.2f \n', t, time);
                fprintf('B_IPS2 relative function value change = %.4e \n', (yold-ynew)/abs(ynew));
                fprintf('B_IPS2 relative grad = %.4e \n', norm(grad,Inf)/max(grad0,1));
                fprintf('Error of estimated beta compared to beta_true = %.8e .\n', error_beta);
            end
            yold=ynew;
        else
            break
        end
        if t == MAX_ITERS;
            fprintf('MAX_ITERS has reached but the BIPS has not yet converged\n')
            break
        end
        if time > MAX_TIME;
            fprintf('MAX_TIME has reached but the BIPS has not yet converged\n')
            break
        end
        tic
    end
end
iter = t;
beta = [beta0;beta];
mu = mu*exp(beta0);
end

function [y,g] = f(beta,values)
%Poisson negative log-liklihood 

%y = q.'*exp(X*beta)-n.'*X*beta;

%Since this function maybe called in the loop for many times, we should not
%evaluate n.'*X each time as it is fixed, therefore, define A=n.'*X, we
%have:

mu=values.q.*exp(values.X*beta);
y = sum(mu,1)-values.A*beta;
if nargout > 1 % gradient required
    g = -((values.n-mu).'*values.X).';
end
end

function [y,g,h] = f_BCD_re(beta,values)

%Reparametrized Poisson Loss

mu=values.Bk.*exp(values.Z*beta);
sum_mu=sum(mu);
y = values.B*log(sum_mu)-values.Ak*beta;
if nargout > 1 % gradient required
%   h21 = (mu.'*values.Z).'; 
    h21 = values.Z.'*mu;
    g = -values.Ak.'+values.B*h21/sum_mu;
    if nargout > 2 % Hessian required
        h1 = bsxfun(@times,values.Z,mu).'*values.Z;
        h = values.B*(h1/sum_mu-h21*h21.'/sum_mu^2);
    end
    
end
end

