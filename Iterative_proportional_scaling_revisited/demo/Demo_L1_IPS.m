clear,  clc
addpath(genpath(pwd))
rng('default')
rng(0)

%% Generate some test data
n = 100;
p = 500;

betaTrue = [[1, 5], zeros(1, p-4),[5, 3]]';
X = [ones(n, 1), round(rand(n, p-1))];
q = rand(n, 1);
tempMean = exp(X * betaTrue + log(q));
y = poissrnd(tempMean); % Poisson counts


%% Specify parameter valus
% Convergence control
convCtrl.ErrBnd = 5e-02;
convCtrl.MaxIt = 1e+03;

%% Initial estimate for the whole sol path
betaInit = zeros(p, 1); 

%% Specify the grid of values for the l1 regularization parameter  
lambdaWeights = sqrt(sum(X.^2))';
lambdaWeights(1) = 0; % assume the first is intercept

lambdaGridsize = 50;
lambdaUBnd = log(max(abs((y - q)' * X(:,2:end)) ./ lambdaWeights(2:end)')); % guaranteed to get a zero estimate
lambdaGrid = exp(linspace(lambdaUBnd, 0, lambdaGridsize));
% Construct the weights for lambda

%% Other parameters in controlling the path run function
nz_Ubnd  = min(size(X,1) * 0.8, size(X, 2)); 
warmstarts = 1; %  0 % 
calib = 1;

%% Run the function to get the l1 solution path
[beta_path, nz_Nums, nz_Supps, lossVals, beta_path_calib, lossVals_calib] = ...
    Func_PoissonL1_Path(X, y, q, lambdaGrid, lambdaWeights, betaInit, nz_Ubnd, convCtrl, calib, warmstarts);
nz_Nums'


grid = log(lambdaGrid);
Coefficients = beta_path(2:end,:).'; % excluding the intercept
n_grid = size(Coefficients);
figure;
plot(grid(1:n_grid),Coefficients,'LineWidth',1);
hold on
% line([grid(optInd) grid(optInd)], ylim, 'LineStyle',':','LineWidth',1,...
%     'Color',[0.5,0.5,0.5])
xlabel('log(\lambda)');
ylabel('Coefficients');
set(gca,'FontName','Helvetica','FontWeight','normal','FontSize',16);
fig = gcf;
set(gca, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 2        );


%% Tuning via Predictive Information Criterion
trnerrs = lossVals_calib'; %  lossVals; % 
    % Below, we assume intercept exists. So we exclude the intercept in the
    % inflation calculation
pens = 1*(1 * nz_Nums + 1 * (nz_Nums-1) .* (1 + log( (size(X, 2)-1)./ max(1, nz_Nums-1) )))'; % the max makes zero as zero eventually
criScores = trnerrs + pens;
optInd = find(criScores  == min(criScores), 1, 'first')

beta_opt = beta_path_calib(:, optInd);
support_opt = nz_Supps{optInd};
numNz_opt = nz_Nums(optInd)

[betaTrue, beta_opt]