clear,  clc

rng('default')
rng(0)
addpath(genpath(pwd))
%% Pre-process
%mypath0 = 'C:\work\Dropbox\main\work\iGGL\Code\'
mypath0 = 'C:\Users\student\Dropbox\iGGL\code\refs\mixed_data\';
cd(mypath0)
dataName =   'newsgroup'


runRScript = 1; % 0 %
resoutput = 1; % 0 %



if resoutput, diary(['matlaboutput_', dataName, '.txt']), end

mypath = [mypath0, dataName, '\'];
data = dlmread([mypath, 'data.txt'], ' ');


k = 50;
mixed_lossType = {'tukey', 'poisson','lorenz'} %{'bernoulli'} %{'poisson'}%{'bernoulli'} %;
mixed_index = {1:k, (k+1):2*k, (2*k+1):3*k};

if strcmpi(dataName,'newsgroup')
    if isempty(setxor(mixed_lossType,{'gaussian'})) ||  isempty(setxor(mixed_lossType,{'tukey'}));
        data = data(:,mixed_index{1});
        mixed_index={1:k};
        
    elseif isempty(setxor(mixed_lossType,{'poisson'}));
        data = data(:,mixed_index{2});
        mixed_index={1:k};
        
    elseif isempty(setxor(mixed_lossType,{'bernoulli'})) ||  isempty(setxor(mixed_lossType,{'lorenz'}));
        data = data(:,mixed_index{3});
        mixed_index={1:k};
    end
end

if resoutput
    fid = fopen([mypath0, 'mypath.txt'], 'wt');
    fprintf(fid, '%s\n', mypath);
    fclose(fid);
    fid = fopen([mypath, 'mypath.txt'], 'wt');
    fprintf(fid, '%s\n', mypath);
    fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,m] = size(data);
alpha = zeros (m,1);
Xi = zeros(n,m);



for i = 1:length(mixed_lossType)
    Y = data(:,mixed_index{i});
    lossType = mixed_lossType{i};
    [n, m] = size(Y);
    switch(lower(lossType))
        case {'bernoulli'}
            Y = double(abs(Y)>1e-6);
        case {'lorenz'}
            Y = double(abs(Y)>1e-6);
            Y = Y * 2 - 1;
        case {'poisson'}
            if ~all(all(mod(Y, 1)==0))
                warning('Data matrix not interger-valued')
                Y = round(Y);
            end
            if min(min(Y)) < 0
                warning('DataS matrix not nonnegative')
                Y = Y - min(min(Y));
            end
        case {'gaussian', 'tukey'}
            % Since this is Gaussian-type data, we can scale it.
            % The purpose is to make W properly smalle (better with largest eigenvalue =1
            
            % % way 1
            % Y = Y / min(sqrt(var(Y)));
            % Way 2
            Y = bsxfun(@rdivide, Y, sqrt(var(Y))); % Use the correlation matrix as suggested by Rothman et al. (2008)
            % % way 3
            % Y_centered = bsxfun(@minus, Y, mean(Y));
            % scale = range(quantile(Y_centered(:), [0.25 1-0.25]))/range(norminv([0.25 1-0.25])); % more robust than sqrt(trace(cov(Y))/m);
            % Y = Y / scale;
            
        otherwise
    end
    data(:,mixed_index{i}) = Y;
    
    if strcmpi(lossType, 'poisson')
        intercept_exist{i} = 0;  % 0; %% whether to account for baseline effect.
    else
        intercept_exist{i} = 1;  % 0; %% whether to account for baseline effect.
    end
    
    
    scale = 1;
    phiVaryFactor = []; %  0.001; %
    switch(lower(lossType))
        case {'bernoulli'}
            inits.alpha = (log(mean(Y)./(1-mean(Y))))'; % zeros(m,1); %  mean(Y)'; %
            inits.Xi = Y; %ones(n,1) * inits.alpha'; %   %zeros(size(Y)); -bad
            inits.W =  diag(1./var(Y));
            
        case {'lorenz'}
            inits.alpha =  mean(Y)'; % (log(mean(Y)./(1-mean(Y))))'; % zeros(m,1);
            inits.Xi = Y; %ones(n,1) * inits.alpha'; %   %zeros(size(Y)); -bad
            inits.W =  diag(1./var(Y));
            
        case {'poisson'}
            inits.alpha = zeros(m,1); %  mean(Y)'; %
            inits.Xi = Y; %bsxfun(@rdivide, Y, sum(Y)); %  ones(n,1) * inits.alpha'; %   %zeros(size(Y)); -bad
            inits.W =  diag(1./var(Y)); %diag(1./var(bsxfun(@rdivide, Y, sum(Y)))); %
            
        case {'gaussian', 'tukey'}
            %         phiVaryFactor = [] ; %0.001; %0.5; % 0.9; %  0.1; % 1;
            if strcmpi(lossType, 'gaussian')
                inits.alpha = mean(Y)';
                inits.Xi = Y; % ones(n,1) * inits.alpha'; %  %zeros(size(Y)); -bad
                inits.W =  diag(1./var(Y));
            else
                inits.alpha = median(Y)';
                %     inits.Xi = ones(n,1) * inits.alpha'; %  Y; %
                
                Y_centered = bsxfun(@minus, Y, inits.alpha'); % Y - inits.Xi; %
                scale = range(quantile(Y_centered(:), [0.25 1-0.25]))/range(norminv([0.25 1-0.25])) % more robust than sqrt(trace(cov(Y))/m);
                
                Y_centered_rob = Y_centered;
                scales = sqrt(var(Y));
                for colInd = 1:size(Y,2)
                    scales(colInd) = range(quantile(Y_centered(:, colInd), [0.25 1-0.25]))/range(norminv([0.25 1-0.25]));
                    Y_centered_rob(:, colInd) = max(-sqrt(2*log(n))*scales(colInd), min(sqrt(2*log(n))* scales(colInd), Y_centered(:, colInd)) );
                end
                inits.Xi = Y_centered_rob + (Y - Y_centered);
                if min(abs(scales))>1e-4
                    inits.W =  diag(1./(scales.^2));
                else
                    inits.alpha = mean(Y)';
                    inits.Xi = Y;
                    inits.W =  diag(1./var(Y));
                end
            end
            
        otherwise
            warning('Which specific case?')
            inits.alpha = mean(Y)';
            inits.Xi = Y; % ones(n,1) * inits.alpha'; %  %zeros(size(Y)); -bad
            inits.W =  diag(1./var(Y));
    end
    alpha(mixed_index{i}) = inits.alpha;
    Xi(:,mixed_index{i}) =  inits.Xi;
end





if strcmpi(dataName,'newsgroup')
    if isempty(setxor(mixed_lossType,{'gaussian'}))  
        inits.phi = 1e-3/norm(inits.W, 2);
        gridsize = 50; grid = exp(linspace(-1, -3, gridsize));
    elseif  isempty(setxor(mixed_lossType,{'tukey'}));
        inits.phi = 1e-4/norm(inits.W, 2);
        gridsize = 50; grid = exp(linspace(-1.5, -3, gridsize));
        
    elseif isempty(setxor(mixed_lossType,{'poisson'}));
        inits.phi = 1e-3/norm(inits.W, 2);
        gridsize = 50; grid = exp(linspace(-2, -4, gridsize));
        
    elseif isempty(setxor(mixed_lossType,{'bernoulli'}));
        inits.phi = 1e-3/norm(inits.W, 2);
        gridsize = 50; grid = exp(linspace(1, -1, gridsize));
        
    elseif isempty(setxor(mixed_lossType,{'lorenz'}));
        inits.phi = 1e-3/norm(inits.W, 2);
        gridsize = 50; grid = exp(linspace(-2, -5, gridsize));
        
        
    elseif isempty(setxor(mixed_lossType,{'gaussian', 'poisson','lorenz'}));
        inits.phi = 1e-3/norm(inits.W, 2);
        gridsize = 50; grid = exp(linspace(-1.25, -4, gridsize));
        inits.alpha = alpha;
        inits.Xi = Xi;
        inits.W = diag(1./var(data));
        
    elseif isempty(setxor(mixed_lossType,{'tukey', 'poisson','lorenz'}));
        inits.phi = 1e-3/norm(inits.W, 2);
        gridsize = 50; grid = exp(linspace(-1.25, -4, gridsize));
        inits.alpha = alpha;
        inits.Xi = Xi;
        inits.W = diag(1./var(data));
    else error('mixed_lossType specified??')
    end
end

%inits.phi = 1e-10/norm(inits.W, 2); % small enough


rho = 1; % fixed. We assume that the loss function is Lip(1)
methodW = 'l1_general'; % 'cvx'; %

convCtrl.outMaxIt = 1e+03;
convCtrl.outErrBnd = 1e-03;

convCtrl.inMaxIt = 1e+1;
convCtrl.inErrBnd = 1e-02;

convCtrl.wMaxIt = 1000;
convCtrl.wErrBnd = 1e-5;
convCtrl

tol = 1e-4

warmstarts = 0; %  0 %

optType =  'Joint' %'W' %'Joint' % 'W' %;


%gridsize = 50; grid =  exp(linspace(-1, -3, gridsize));


phiVaryFactor, scale, warmstarts
Y = data;
[n, m] = size(Y);

% No need of scaling. This parameter grid is for regularization on the inverse covariance directly.
nz_Ubnd = min(m*m*.75, n*m/2); % At most 80% connected.


calib = 1;

% warning off;
[alpha_path, W_path, nz_Nums, nz_Supps, objVals, alpha_path_calib, W_path_calib, objVals_calib] = ...
    Func_iGGL_Path(Y, optType, intercept_exist, inits, grid, nz_Ubnd, convCtrl, rho, methodW, mixed_lossType,...
    scale, calib, warmstarts, phiVaryFactor, tol, mixed_index);
% warning on;

nz_Nums'

objVals_ref = objVals_calib; %  objVals; %
trnerrs = reshape(cell2mat(objVals_ref), [size(objVals_ref{1}), numel(objVals_ref)]);
trnerrs = reshape(trnerrs(2,:,:), 1, numel(objVals_ref));
pens = (1 * nz_Nums + 1 * nz_Nums .* (1 + log( (size(Y, 2)*(size(Y, 2)-1))./ max(1, nz_Nums) )))'; % the max makes zero as zero eventually
criScores = trnerrs + pens;
optInd = find(criScores  == min(criScores), 1, 'first')

alpha_opt = alpha_path_calib{optInd};
W_opt = W_path_calib{optInd};
support_opt = nz_Supps{optInd};
numNz_opt = nz_Nums(optInd)



sum(sum(abs(W_opt-W_opt.')))


if strcmpi(dataName,'newsgroup')
    if isempty(setxor(mixed_lossType,{'gaussian'}))
        figure;
        subplot(1,2,1);
        imagesc(abs(W_opt)>tol);
        colormap('Gray')
        subplot(1,2,2);
        imagesc(abs(inv(cov(Y)))>850*tol);
        colormap('Gray')
    end
end

if resoutput
    dlmwrite([mypath, sprintf('WEst_%s.txt',strjoin(mixed_lossType,'_'))], W_opt, ',');
    %     varnames = cellstr(num2str([1:m]'));
    %     fid = fopen([mypath, 'varnames.txt'], 'wt');
    %     fprintf(fid, '%s\n', varnames{:})
    %     fclose(fid);
end

if resoutput
    fid = fopen([mypath, 'plotfile.txt'], 'wt');
    fprintf(fid, '%s\n', sprintf('WEst_%s.txt',strjoin(mixed_lossType,'_')));
    fclose(fid);
end

if runRScript
    cmdline = ['R CMD BATCH ', mypath, 'plotgraph.R ', mypath, 'tempout_R.txt'];
    status  = system(cmdline, '-echo')
    if status ~= 0, warning('R program not run successfully. Check tempout_R for more details'), end
end

if resoutput, diary off, end
if resoutput
    save([mypath, dataName, datestr(clock, '-yyyy-mm-dd-HH.MM.SS-'),strjoin(mixed_lossType,'_'), '.mat'])
end



phiVaryFactor, warmstarts, inits.phi*norm(inits.W,2), grid
