clear,  clc
addpath(genpath(pwd))
rng('default')
rng(0)


%% Input the design matrix from
mypath0 = 'C:\Users\student\Dropbox\project1_IPFP\IPFP\code\real data';
%mypath0 = 'C:\Users\shao.tang\Dropbox\project1_IPFP\IPFP\code\real data';
cd(mypath0)
option='twoway';
path='C:\Users\student\Dropbox\project1_IPFP\IPFP\code\real data\bank marketing\aggregate_zero_margin_removed';
 
mydata=sparse(csvread(fullfile(path,sprintf(('%s_design_matrix.dat'),option))));

X=mydata(:,1:end-1);
y=mydata(:,end);

nz_index1=find(sum(X)); %index of non-zero columns
nz_index2=find(X.'*y~=0); %index of non-zero mariginal total columns

nz_index=nz_index2; %intercetion of two sets is nz_index2

fprintf('Total number of predictors including itercept is %d \n.',length(nz_index))


% X=X0(:,any(X0,1)); %delete zero columns
% 
% X(:,find(X.'*y==0))=[];


n = size(X,1);
p = size(X,2);


q = ones(n, 1);




delete(fullfile(path,sprintf(('bank_%s_logfile_modified_full_path.txt'),option)));
diary(fullfile(path,sprintf(('bank_%s_logfile_modified_full_path.txt'),option)));

%% Specify parameter valus
% Convergence control
convCtrl.ErrBnd = 1e-03;
convCtrl.MaxIt = 1e+03;
convCtrl

calib_options.MaxIter=1e+4;
calib_options.MaxFunEvals=1e+4;
        
%% Initial estimate for the whole sol path
betaInit = zeros(p, 1); 

%% Specify the grid of values for the l1 regularization parameter  
lambdaWeights = sqrt(sum(X.^2))';
lambdaWeights(1) = 0; % assume the first is intercept

lambdaGridsize = 100;
lambdaUBnd = log(max(abs((y - q)' * X(:,2:end)) ./ lambdaWeights(2:end)')); % guaranteed to get a zero estimate
lambdaGrid = exp(linspace(lambdaUBnd, -4, lambdaGridsize));
% Construct the weights for lambda

%% Other parameters in controlling the path run function
nz_Ubnd  = min(size(X,1) * 0.8, size(X, 2)); 
warmstarts = 1; %  0 % 
calib = 1;

%% Run the function to get the l1 solution path

tic
[beta_path, nz_Nums, nz_Supps, lossVals, beta_path_calib, lossVals_calib] = ...
    Func_PoissonL1_Path_modified(X, y, q, lambdaGrid, lambdaWeights, betaInit, nz_Ubnd, convCtrl, calib, warmstarts,[]);

fprintf('The running time is %.2f seconds \n', toc)
nz_Nums'
subplot(1,2,1)
plot(beta_path(2:end,:)') % excluding the intercept
subplot(1,2,2)
plot(beta_path_calib(2:end,:)') % excluding the intercept

%% Tuning
trnerrs = lossVals_calib'; %  lossVals; % 
    % Below, we assume intercept exists. So we exclude the intercept in the
    % inflation calculation
pens = 1*(1 * nz_Nums + 1 * (nz_Nums-1) .* (1 + log( (size(X, 2)-1)./ max(1, nz_Nums-1) )))'; % the max makes zero as zero eventually
criScores = trnerrs + pens;
optInd = find(criScores  == min(criScores), 1, 'first')

beta_opt = beta_path_calib(:, optInd);
support_opt = nz_Supps{optInd};
numNz_opt = nz_Nums(optInd)


%[betaTrue, beta_opt]

%% Extract non-zero features on the original index before the zero columns are removed

beta_opt(beta_opt~=0);
nz_index(beta_opt~=0);

fid = fopen( fullfile(path,sprintf(('L1_feature_%s_full_path.txt'),option)), 'wt');
fprintf(fid,'%d\n',nz_index(beta_opt~=0).');
fclose(fid);

save(fullfile(path,sprintf(('%s_result_full_path.mat'),option)));
diary off