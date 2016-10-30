function [alpha_path, W_path, nz_Nums, nz_Supps, objVals, alpha_path_calib, W_path_calib, objVals_calib] = ...
    Func_iGGL_Path(Y, optType, intercept_exist, inits, grid, nz_Ubnd, convCtrl, rho, methodW, lossType, scale, calib, warmstarts,...
    phiVaryFactor,tol, mixed_index)
%==========================================================================
% This function calculates the solution path corresponding to lambda_W
%------------------------ Input Variables ---------------------------------
% Y               - Target data matrix
% optType         - which variables to be optimized
% intercept_exist - Whether intercept is considered in the optimization
% inits           - Provide some parameter values and some starting values of unknown parameters. 
% grid            - A vector of values for lambda_W
% nz_Ubnd         - Upperbound for number of nonzeros to stop computing the path
% convCtrl        - Provide control of the stopping criterion in iterations.
% rho             - fixed at 1 in this code to keep things simple
% methodW         - Technique used to optimize W
% lossType        - loss function
% scale           - Only used in robust loss functions, to provide a robust scale estimate
% calib           - Perform bias correction or not. If =1, refit a model restricted to the learned support.
%------------------------ Output Variables ---------------------------------
% alpha_path, W_path, alpha_path_calib, W_path_calib: solution paths for alpha, W, and calbrated A and W 
% nz_Nums:  vector of number of nonzeros in the off-diagonal part of W
% nz_Supps: support sets of the off-diagonal part of W estimates
%==========================================================================
output = 1;
inits0 = inits;
% warmstarts = 1;

grid = sort(grid, 'descend');
gridsize = length(grid);    
nz_Nums = zeros(gridsize, 1);
nz_Supps = cell(gridsize, 1);
objVals = cell(gridsize, 1);
objVals_calib = cell(gridsize, 1);
alpha_path = cell(gridsize, 1);
W_path = cell(gridsize, 1);
alpha_path_calib = cell(gridsize, 1);
W_path_calib = cell(gridsize, 1);

for j = 1:gridsize
    lambda_W = grid(j);
    if output == 1; disp(['=========== Parameter index: ', num2str(j), ', lambda=', num2str(lambda_W), ' ===========']), end
    
    Mask = ones(size(Y,2)); Mask = Mask - diag(diag(Mask)); 
        % The weighting matrix to be multiplied with lambda_W. 1: possibly nonzero; 0: no sparsity imposed
    Lambda_W_Mat = lambda_W .* Mask;
    
%     inits.alpha', inits.W(1:5,1:5), inits.Xi(1:5, 1:5)
    
    [alpha_path{j}, W_path{j}, Xi_Est, objVals{j}] = Func_iGGL_AOS(Y, optType, intercept_exist, ...
        inits, convCtrl, rho, Lambda_W_Mat, methodW, lossType, scale, phiVaryFactor, mixed_index);
    
%     alpha_path{j}', W_path{j}(1:5,1:5), Xi_Est(1:5, 1:5)
    

    [nz_Nums(j), nz_Supps{j}] = Func_CalcNz(W_path{j}, 1, tol);  % 1 means we do not want to count the diagonal entries
        % number of nonzeros on off-diagonal entries. Note that here, we
        % scale (back) the result. 
    if output == 1; disp(['***** # of nzs (off-diag): ', num2str(nz_Nums(j)), ' *****']), end
        
    if calib == 1 % perform bias correction, by restricting a nonzeros to the learned support 
        disp('....Calibration....')
        Mask2 = ones(size(Y,2)); Mask2 = Mask2 - diag(diag(Mask2)); 
        Mask2(nz_Supps{j}) = 0;
        Lambda_W_Mat2 =  (1e+20) * Mask2;
        inits2 = inits; inits2.alpha = alpha_path{j}; inits2.W = W_path{j}; init2.Xi = Xi_Est;
        convCtrl2 = convCtrl;
        convCtrl2.outMaxIt = convCtrl.outMaxIt / 10;
        [alpha_path_calib{j}, W_path_calib{j}, ~, objVals_calib{j}] = Func_iGGL_AOS(Y, optType, intercept_exist, ...
            inits2, convCtrl2, rho, Lambda_W_Mat2, methodW, lossType, scale, phiVaryFactor, mixed_index);
        if output == 1; disp(['***** # of nzs (off-diag): ', num2str(nz_Nums(j)), ' *****']), end
    
    end
        
    if warmstarts == 1
    % Use warm starts 
        inits.alpha = alpha_path{j};
        inits.W = W_path{j};
        inits.Xi = Xi_Est;
    else
        inits = inits0;
    end
    
    % Early stopping
    if nz_Nums(j) > nz_Ubnd
        disp('Too many nonzeros')
        alpha_path = alpha_path(1:j);
        W_path = W_path(1:j);
        nz_Nums = nz_Nums(1:j); 
        nz_Supps = nz_Supps(1:j);
        objVals = objVals(1:j);
        alpha_path_calib = alpha_path_calib(1:j);
        W_path_calib = W_path_calib(1:j);
        objVals_calib = objVals_calib(1:j);
        
        gridsize = j;
        break;
    end
end
