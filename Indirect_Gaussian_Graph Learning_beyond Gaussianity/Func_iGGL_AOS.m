function [alpha_Est, W_Est, Xi_Est, objVal] = Func_iGGL_AOS(Y, optType, intercept_exist, inits, convCtrl, rho, Lambda_W_Mat, methodW, mixed_lossType, scale, ...
    phiVaryFactor, mixed_index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function implements the iterative GGL algorithm.
% Input:
% inits: Initial values to be used in iteration--inits.Xi, inits.W, inits.alpha
% Output:
% alpha_Est, W_Est, Xi_Est: estimates.
% objVal: a column vector with the overall f-value, the loss value without penalty, and the phi part only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assume rho = 1
if rho ~= 1
    error('Not implemented yet.')
end
output = 1; %1;
debug = 1; %1;

% ============== Initialization ==============
[n,m] = size(Y);

alpha_cur = inits.alpha;
W_cur = inits.W;
Xi = inits.Xi; % changing over time (did not give record its value though)
phi = inits.phi; % initial

M_cur = ones(n,1) * alpha_cur';
Theta_cur = phi * (M_cur - Xi) * (W_cur) + Xi;


if debug
    disp(['norm and phi: ', num2str(norm(W_cur, 2)), ', ', num2str(phi)])
end

if debug == 1
    funcval_cur = Inf;
    %     [tempvals] = Func_Obj_L1(W_cur, M_cur, Lambda_W_Mat, phi, Xi, Y, Theta_cur, lossType, scale);
    %     funcval_cur = tempvals(1);
    %     % not useful as a good reference when lambda values changes a lot.
    %     % Because the W estimate may not correspond to the given Xi at the
    %     % correct threshold level
end

% ============== Outer Loop ==============
t = 0;
while 1
    t = t + 1;
    if output == 1
        disp(['Outer iter=', num2str(t)])
    end
    % Carry out the surrogate step
    gradval = zeros(n,m);
    for i = 1:length(mixed_lossType)
        lossType = mixed_lossType{i};
        [~, gradval(:,mixed_index{i})] = Func_Loss_Grad(Y(:,mixed_index{i}), Theta_cur(:,mixed_index{i}), lossType, 1, scale);
    end
    
    Xi = Theta_cur - gradval;
    if debug
        disp(['Max diff between Xi and Y--  ', num2str(max(max(abs(Xi - Y))))])
    end
    % ========Inner Loop ==============
    % Here, Xi is given. We can update the mean and W alternatively.
    % But if optType is NOT 'joint', the inner loop stops after 1 iteration.
    alpha_cur_inn = alpha_cur; W_cur_inn = W_cur;
    alpha_new_inn = zeros(m,1);
    tt = 0;
    if output == 1, fprintf('InnerIter='), end
    while 1
        tt = tt+1;
        if output == 1, fprintf([num2str(tt),',']), end
        % Update mean if required
        if strcmpi(optType, 'M') || strcmpi(optType, 'Joint')
            pseudoDesign = ones(n,1);
            for i = 1:length(mixed_lossType)
                if intercept_exist{i}
                    alpha_new_inn(mixed_index{i}) = ((pseudoDesign' * pseudoDesign) \ (pseudoDesign' * Xi(:, mixed_index{i})))';
                else
                    alpha_new_inn(mixed_index{i}) = alpha_cur_inn(mixed_index{i});
                end
            end
        else
            alpha_new_inn = alpha_cur_inn;
        end
        
        
        
        
        % Update W if required
        M_new_inn  = ones(n,1) * alpha_new_inn';
        S = (Xi - M_new_inn)' * (Xi - M_new_inn) / (n);
        
        if strcmpi(optType, 'W') || strcmpi(optType, 'Joint')
            [W_new_inn] = Func_MyGLasso(S, Lambda_W_Mat, W_cur_inn, convCtrl, methodW);
        else
            W_new_inn = W_cur_inn;
        end
        
        err_alpha_inn = max(max(abs(alpha_new_inn - alpha_cur_inn)));
        err_W_inn = max(max(abs(W_new_inn - W_cur_inn)));
        if tt >= convCtrl.inMaxIt || err_alpha_inn <= convCtrl.inErrBnd || err_W_inn <= convCtrl.inErrBnd
            break;
        else
            alpha_cur_inn = alpha_new_inn;
            W_cur_inn = W_new_inn;
        end
    end
    if output == 1, fprintf('\n'), end
    % ========Inner Loop ends==============
    
    % Update Theta
    
    alpha_new = alpha_new_inn;
    W_new = W_new_inn;
    
    if isempty(phiVaryFactor)
        phi = inits.phi;
    else
        phi = phiVaryFactor/norm(W_new, 2);
    end
    
    if debug
        %     disp(num2str(alpha_new'))
        disp(['norm and phi: ', num2str(norm(W_new, 2)), ', ', num2str(phi), ...
            '; num of direct nzs: ', num2str(sum(sum(abs(W_new - diag(diag(W_new)))> 0 )))])
    end
    
    %     diag(W_new)', max(diag(W_new))
    if norm(W_new, 2) * phi > 1
        warning('Spectral constraint (upper) not satisfied. If this occurs at a later stage, perhaps need to decrease phi further for the given value of lambda')
    end
    
    M_new  = ones(n,1) * alpha_new';
    Theta_new = phi * (M_new - Xi) * (W_new) + Xi;
    
    err_Theta = max(max(abs(Theta_new-Theta_cur)));
    err_W = max(max(abs(W_new-W_cur)));
    
    if debug == 1
        [tempvals] = Func_Obj_L1(W_new, M_new, Lambda_W_Mat, phi, Xi, Y, Theta_new, mixed_lossType, scale, mixed_index);
        funcval_new = tempvals(1);
        if output == 1
            disp(['func val=', num2str(funcval_new), ', decreasing = ', num2str(funcval_new - funcval_cur), ', err_W=', num2str(err_W), ', err_Theta=', num2str(err_Theta)])
            if isempty(phiVaryFactor) % not check this when using varying phi
                if funcval_new>funcval_cur + 1e-3
                    error('func-val increasing')
                end
            end
        end
        funcval_cur = funcval_new;
    end
    
    if t >= convCtrl.outMaxIt || err_Theta < convCtrl.outErrBnd || err_W < convCtrl.outErrBnd
        [objVal] = Func_Obj_L1(W_new, M_new, Lambda_W_Mat, phi, Xi, Y, Theta_new, mixed_lossType, scale, mixed_index);
        break;
    else
        Theta_cur = Theta_new;
        alpha_cur = alpha_new;
        W_cur = W_new;
    end
end
alpha_Est = alpha_new;
W_Est = W_new;
Xi_Est = Xi;

function [fvalvec] = Func_Obj_L1(W, M, Lambda_W_Mat, phi, Xi, Y, Theta, mixed_lossType, scale, mixed_index)
% fval: the overall, fval1: the loss without penalty, fval0: the phi part only
fval0 = 0;
for i = 1:length(mixed_lossType)
    lossType = mixed_lossType{i};
    [fval0_partial, ~] = Func_Loss_Grad(Y(:,mixed_index{i}), Theta(:,mixed_index{i}), lossType, 0, scale);
    fval0 = fval0 + fval0_partial;
end


[n,m] = size(Y);
temp = W - phi * W * W;
if n <= m
    fval1 = fval0 / phi + trace((Xi - M) * temp * (Xi - M)') / 2 -  logdet(W, 'chol') * n /2; %log(det(W)) * n /2; %
else
    fval1 = fval0 / phi + trace((Xi - M)' *(Xi - M)* temp  ) / 2 -  logdet(W, 'chol') * n /2; %log(det(W)) * n /2; %
end
fval = fval1 + (n /2) * sum(sum( abs(Lambda_W_Mat .* W) ));
fvalvec = [fval, fval1, fval0]';
