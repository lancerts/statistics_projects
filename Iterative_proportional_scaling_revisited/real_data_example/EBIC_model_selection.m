clear,  clc
addpath(genpath(pwd))
rng('default')
rng(0)

mypath0 = 'C:\Users\student\Dropbox\project1_IPFP\IPFP\code\real data';
%mypath0 = 'C:\Users\shao.tang\Dropbox\project1_IPFP\IPFP\code\real data';
cd(mypath0)
%option='twoway';
option='threeway';
path='C:\Users\student\Dropbox\project1_IPFP\IPFP\code\real data\bank marketing\aggregate_zero_margin_removed\EBIC_analysis';

load(fullfile(path,sprintf(('%s_result.mat'),option)));
%load(fullfile(path,sprintf(('%s_result_full_path.mat'),option)));

%% Tuning
trnerrs = 2*lossVals_calib'; %  lossVals; % 
gamma = 1;
pens =  nz_Nums.'*(log(size(X, 1)) + 2*gamma*log(size(X, 2))); 
criScores = trnerrs + pens;
optInd = find(criScores  == min(criScores), 1, 'first')

beta_opt = beta_path_calib(:, optInd);
support_opt = nz_Supps{optInd};
numNz_opt = nz_Nums(optInd)


path='C:\Users\student\Dropbox\project1_IPFP\IPFP\code\real data\bank marketing\aggregate_zero_margin_removed\EBIC_analysis';
fid = fopen( fullfile(path,sprintf(('L1_EBIC_feature_%s.txt'),option)), 'wt');
fprintf(fid,'%d\n',nz_index(beta_opt~=0).');
fclose(fid);



