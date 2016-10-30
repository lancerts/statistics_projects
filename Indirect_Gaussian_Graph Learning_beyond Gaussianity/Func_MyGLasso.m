function [W_est] = Func_MyGLasso(S, Lambda_W_Mat, initW, convCtrl, methodW)

W_cur = initW;
m = size(S, 1);
switch lower(methodW)
    case 'cvx'
        cvx_quiet(1) ;
        warning off;
        cvx_begin sdp
            variable W(m, m) semidefinite
            minimize( - log_det(W) + trace(S*W) + sum(sum( abs(Lambda_W_Mat .* W) )) )  
%             subject to 
%                 W <= K * eye(m);
        cvx_end
        warning on;
        W_est = W;
    case 'l1_general'
        options = struct('verbose', 0, 'maxIter', convCtrl.wMaxIt, 'optTol', convCtrl.wErrBnd);
        
        indices = (1:m*m)';
        funObj = @(x)sparsePrecisionObj(x, m, indices, S);
        
        Lamda_W_vec = Lambda_W_Mat(:);
        
        wEst = L1General2_PSSgb(funObj, initW(:), Lamda_W_vec, options);
        W_est = reshape(wEst, m, m);

    otherwise
        disp('not implemented');
end

