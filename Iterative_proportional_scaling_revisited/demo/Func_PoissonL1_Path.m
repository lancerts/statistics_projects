function [beta_path, nz_Nums, nz_Supps, lossVals, beta_path_calib, lossVals_calib] = ...
    Func_PoissonL1_Path(X, y, q, lambdaGrid, lambdaWeights, betaInit, nz_Ubnd, convCtrl, calib, warmstarts,calib_options)
%==========================================================================
% This function calculates the solution path of the l1-penalized Poisson model
%------------------------ Input Variables ---------------------------------
%
%------------------------ Output Variables ---------------------------------
%
%==========================================================================
output = 1;

betaInit0 = betaInit;

lambdaGrid = sort(lambdaGrid, 'descend');
lambdaGridsize = length(lambdaGrid);    
nz_Nums = zeros(lambdaGridsize, 1);
nz_Supps = cell(lambdaGridsize, 1);
lossVals = zeros(lambdaGridsize, 1);
lossVals_calib = zeros(lambdaGridsize, 1);
beta_path = zeros(size(X,2), lambdaGridsize);
beta_path_calib = beta_path;


for j = 1:lambdaGridsize
    lambda = lambdaGrid(j);
    if output == 1; disp(['=========== Parameter index: ', num2str(j), ', lambda=', num2str(lambda), ' ===========']), end
    
    Lambda = lambdaWeights * lambda; % same size as beta
    
    
    [beta_path(:,j), ~, lossVals(j), ~] = PoissonL1(X, y, q, Lambda, convCtrl.ErrBnd, convCtrl.MaxIt, betaInit);
    
    tol = 1e-4;
    nz_Supps{j} = find(abs(beta_path(:,j)) > tol);
    nz_Nums(j) = numel(nz_Supps{j});
    
    if output == 1; disp(['***** # of nzs: ', num2str(nz_Nums(j)), ' *****']), end
        
    if calib == 1 % perform bias correction, by restricting a nonzeros to the learned support 
        disp('....Calibration....')

        % For tuning purposes, we usually want to calibrate the estimate. 
        % That is, the preferred estimate is an unbiased estimate on the selected
        % dimensions. This can be achieved by setting a proper lambda vector

        X_calib=X(:,nz_Supps{j});
        
        %Calibration with Newton if the number of the selected nonzero coeeficients
        %is less than 400. Use LBFGS otherwise.
        calib_options.display = 0;
        if size(X_calib,2)<400
            calib_options.Method ='newton';
        else
            calib_options.Method ='lbfgs';
        end
        g=@(x) PoissonObj(X_calib, y, q, x);
        beta_non_zero=beta_path(nz_Supps{j}, j);
        [beta_non_zero_calib,lossVals_calib(j),~,~]=minFunc(g,beta_non_zero,calib_options); 
        beta_path_calib(nz_Supps{j}, j)=beta_non_zero_calib;

    end
        
    if warmstarts == 1
    % Use warm starts 
        betaInit = beta_path(:, j);
    else
        betaInit = betaInit0;
    end
    
    % Early stopping
    if nz_Nums(j) > nz_Ubnd
        disp('Too many nonzeros')
        beta_path = beta_path(:, 1:j);
        nz_Nums = nz_Nums(1:j); 
        nz_Supps = nz_Supps(1:j);
        lossVals = lossVals(1:j);
        beta_path_calib = beta_path_calib(:, 1:j);
        lossVals_calib = lossVals_calib(1:j);
        
        lambdaGridsize = j;
        break;
    end
end
end



function [fval,g,h] = PoissonObj(X, y, q,beta)
%definition of the function used for minimization, negative log-liklihood:

%have:
eta = X * beta;
mu=q .* exp(eta);
fval = - sum(y .* eta) + sum(mu);

if nargout > 1 % gradient required
    g = -X.'*(y-mu);
    if nargout > 2 % Hessian required
        h = bsxfun(@times,X,mu).'*X;  
    end
end
end
