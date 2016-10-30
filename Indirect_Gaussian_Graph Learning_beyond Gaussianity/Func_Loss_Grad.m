function [F, G] = Func_Loss_Grad(Y, Theta, lossType, gradDesired, scale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the loss function value and its gradient function value
% This is the loss of the overall entries under column-dependency.
% Y denotes the oberved data matrix. Theta is the estimate which may NOT on
% the same scale of Y.
% The penalty part is NOT included.
% If gradDesired=0: only the loss function; 1: only the grad; 2: both (default). 
% The scale parameter is used in the robust loss function (only); by default, it is 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('scale', 'var') || isempty('scale')
    scale = 1;
end
if ~exist('gradDesired', 'var') || isempty('gradDesired')
    gradDesired = 2;
end

debug = 0;%1;
F = []; G = [];
switch lower(lossType)
    case 'gaussian'
        if gradDesired ~= 1 
            F = norm(Y-Theta, 'fro')^2/2;
        end
        if gradDesired ~= 0
            G = Theta - Y;
        end
    case {'bernoulli'}
        % Don't directly calculate log(1+exp(Theta)).
        % Theta may have very large values, in which case the quantity should be essentially 1
        Lfactor = 1/4;
        if gradDesired ~= 1 
            F = (-sum(sum(Y.*Theta)) + sum(sum( Theta(Theta>500) )) + sum(sum( log(1+exp( Theta(Theta<=500) )) ))) /Lfactor;
        end
        if gradDesired ~= 0
            G = ( - Y +  1./(1+exp(-Theta)) ) /Lfactor;
            % better than exp(Theta) ./ (exp(Theta) + 1);  due to occuring of INF! 
        end
    case 'poisson' % re-parametrized poisson loss  
        Lfactor = norm(Y, 1)/2;
        if gradDesired ~= 1 
%         	F = ( -sum(sum(Y.*Theta)) + sum( sum(Y) .* log(sum(exp(Theta))) ) ) / Lfactor; 
                % This is easy but may not be scalable when say >800
            colMax = max(Theta);
            Theta2 = bsxfun(@minus, Theta, colMax); %Theta - repmat(colMax, size(Theta,1), 1);
            F = ( -sum(sum(Y.*Theta)) + sum( sum(Y) .* (colMax + log(sum(exp(Theta2)))) ) ) / Lfactor; 
        end
        if gradDesired ~= 0
            colMax = max(Theta);
            Theta2 = bsxfun(@minus, Theta, colMax); %Theta - repmat(colMax, size(Theta,1), 1);
            G = ( - Y +  bsxfun(@times, exp(Theta2), sum(Y)./sum(exp(Theta2)))  ) /Lfactor;
            %G = ( - Y +  repmat(sum(Y)./sum(exp(Theta2)), size(Theta,1), 1) .* exp(Theta2) ) /Lfactor;
        end
    case 'huber'
        lambda = 1.345 * scale;
        Lfactor = 1;
        R = Theta - Y;
        ttInds1 = find(abs(R) <= lambda);
        ttInds2 = find(abs(R) > lambda);
    
        if gradDesired ~= 1 
            F = zeros(size(Theta));
            F(ttInds1) = R(ttInds1).^2 / 2;
            F(ttInds2) = lambda * abs(R(ttInds2)) - lambda^2/2;
            F = sum(sum(F)) /Lfactor;
        end
        if gradDesired ~= 0
            G = zeros(size(Theta));
            G(ttInds1) = R(ttInds1);
            G(ttInds2) = lambda * sign(R(ttInds2));
            G = G /Lfactor;
        end
        
    case 'tukey'
        lambda = 4.685 * scale;
        Lfactor = 1;
        R = Theta - Y;
        ttInds1 = find(abs(R) <= lambda);
        ttInds2 = find(abs(R) > lambda);
        if debug
            disp(['Current # of gross outlying entries: ' num2str(length(ttInds2)), ...
                ', # of at least moderate ones: ', num2str(sum(sum(abs(R)>lambda/sqrt(5)))), ...
                ', range of R: ', num2str(max(range(R))) ])
        end
        if gradDesired ~= 1 
            F = zeros(size(Theta));
            F(ttInds1) = 1 - (1 - (R(ttInds1)/lambda).^2).^3;
            F(ttInds2) = 1;
            F = sum(sum(F)) * lambda^2 /6 /Lfactor;
        end      
        
        if gradDesired ~= 0
            G = zeros(size(Theta));
            G(ttInds1) = R(ttInds1) .* ((1 - (R(ttInds1)/lambda).^2 ).^2);
                %or R(ttInds1) - 2 * R(ttInds1).^3 ./ lambda^2 + R(ttInds1).^5 ./ lambda.^4; 
            G = G / Lfactor;
        end
    case 'lorenz'
        if any(any(abs(Y)~=1)), error('Wrong input values for Lorenz loss'), end
        Lfactor = 2;
        R = Theta .* Y;
        ttInds1 = find((R) <= 1);
        ttInds2 = find((R) > 1);
        
        if gradDesired ~= 1 
            F = zeros(size(Theta));
            F(ttInds1) = log(1 + (R(ttInds1) - 1).^2);
            F = sum(sum(F)) /Lfactor;
        end      
        
        if gradDesired ~= 0
            G = zeros(size(Theta));
            G(ttInds1) = ( 2*(Theta(ttInds1) - Y(ttInds1)) ) ./ (1 + (R(ttInds1) - 1).^2);
            G = G / Lfactor;
        end
        
    case 'hampel'
        error('To be implemented')
    otherwise
        error('Not implemented yet')
end

end