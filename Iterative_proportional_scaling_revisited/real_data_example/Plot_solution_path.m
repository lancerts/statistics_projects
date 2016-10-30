clear,  clc
addpath(genpath(pwd))
rng('default')
rng(0)

mypath0 = 'C:\Users\student\Dropbox\project1_IPFP\IPFP\code\real data';
%mypath0 = 'C:\Users\shao.tang\Dropbox\project1_IPFP\IPFP\code\real data';
cd(mypath0)
option='threeway';
path='C:\Users\student\Dropbox\project1_IPFP\IPFP\code\real data\bank marketing\aggregate_zero_margin_removed';

%load(fullfile(path,sprintf(('%s_result.mat'),option)));
load(fullfile(path,sprintf(('%s_result_time_1.000000e-01.mat'),option)));
grid = log(lambdaGrid);
Coefficients = beta_path(2:end,:).';




figure;
plot(grid,Coefficients,'LineWidth',1);
hold on
% line([grid(optInd) grid(optInd)], ylim, 'LineStyle',':','LineWidth',1,...
%     'Color',[0.5,0.5,0.5])
xlabel('log(\lambda)');
ylabel('Coefficients');
set(gca,'FontName','Helvetica','FontWeight','normal','FontSize',16);
xlim([-2 3.5])
fig = gcf;
set(gca, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 2        );
set(gcf, 'PaperPosition', [0 0 6 5.5]);
set(gcf, 'PaperSize', [5.7 5.4]);
save_address=fullfile(path,sprintf('%s_solution_path_full_path.pdf',option));
saveas(fig, save_address)



figure;
plot(grid(1:75),Coefficients(1:75,:),'LineWidth',1);
hold on
% line([grid(optInd) grid(optInd)], ylim, 'LineStyle',':','LineWidth',1,...
%     'Color',[0.5,0.5,0.5])
xlabel('log(\lambda)');
ylabel('Coefficients');
set(gca,'FontName','Helvetica','FontWeight','normal','FontSize',16);
xlim([0.5 3.5])
fig = gcf;
set(gca, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 2        );
set(gcf, 'PaperPosition', [0 0 6 5.5]);
set(gcf, 'PaperSize', [5.7 5.4]);
save_address=fullfile(path,sprintf('%s_solution_path_zoom_full_path.eps',option));
print(fig, '-depsc', save_address)
 








option='threeway';

load(fullfile(path,sprintf(('%s_result.mat'),option)));
grid = log(lambdaGrid);
Coefficients = beta_path(2:end,:).';




figure;
plot(grid,Coefficients,'LineWidth',1);
hold on
% line([grid(optInd) grid(optInd)], [-1,3], 'LineStyle',':','LineWidth',1,...
%     'Color',[0.5,0.5,0.5])
xlabel('log(\lambda)');
ylabel('Coefficients');
set(gca,'FontName','Helvetica','FontWeight','normal','FontSize',16);
xlim([0.5 3])
fig = gcf;
set(gca, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 2        );
set(gcf, 'PaperPosition', [0 0 6 5.5]);
set(gcf, 'PaperSize', [5.7 5.4]);
set(gcf,'GraphicsSmoothing','on')
save_address=fullfile(path,sprintf('%s_solution_path.pdf',option));
print(fig, '-dpdf', '-r300', save_address)



figure;
plot(grid,Coefficients,'LineWidth',1);
hold on
% line([grid(optInd) grid(optInd)], [-1,1], 'LineStyle',':','LineWidth',1,...
%     'Color',[0.5,0.5,0.5])
xlabel('log(\lambda)');
ylabel('Coefficients');
set(gca,'FontName','Helvetica','FontWeight','normal','FontSize',16);
xlim([0.5 3])
ylim([-1 1])
fig = gcf;
set(gca, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 2        );
set(gcf, 'PaperPosition', [0 0 6 5.5]);
set(gcf, 'PaperSize', [5.7 5.4]);
set(gcf,'GraphicsSmoothing','on')
save_address=fullfile(path,sprintf('%s_solution_path_zoom.pdf',option));
print(fig, '-dpdf', '-r300', save_address)
print(fig, '-depsc', save_address)
 
 