% plot a logistic psychometric function for all subjects
% Anne Urai, 2015

clear; clc; close all;
figpath = '~/Dropbox/Figures/uncertainty';

subjects = 1:27;
% all the coherence level
allcohs = [-0.3 -0.2 -0.1 -0.05 -0.025 -0.0125 -0.0063 0.0063 0.0125 0.025 0.05 0.1 0.2 0.3];
newx = linspace(min(allcohs), max(allcohs), 100);
cols = linspecer(5); cols = cols([1 4], :);

% preallocate
grandavg.betas = nan(length(subjects),  2);
grandavg.yfit = nan(length(subjects), 121);
grandavg.xpts = nan(length(subjects), length(allcohs));
grandavg.ypts = nan(length(subjects), length(allcohs));

for sj = unique(subjects),
    
    disp(sj);
    % get all the data
    data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
    
    % make it such that we can plot the GA of datapoints
    LogisticFit = Psychometric_Logistic(data, 0, 0);
    for j = 1:length(allcohs),
        idx = find( abs(LogisticFit.coherences - allcohs(j)) < 0.0001);
        try
            grandavg.xpts(find(sj==subjects),  j) = LogisticFit.coherences(idx);
            grandavg.ypts(find(sj==subjects),  j) = LogisticFit.meanresp(idx);
        end
    end

    % do a quicker logistic fit without the lapse rate 
    x = data.coherence .* data.stim;
    y = data.resp; y(y==-1) = 0;
    
    % sort both
    [x, idx] = sort(x);
    y = y(idx);
    
    [beta,dev,stats] = glmfit(x, y, ...
        'binomial','link','logit');
    
    yfit = glmval(beta, newx,'logit');
    grandavg.betas(sj, :) = beta;
    
    % also save the curve itself
    grandavg.yfit(sj, :) = LogisticFit.y;
    
    
    if 0,
        % also get the individual d prime
        % use 0 (for weaker motion answer) instead of -1
        resp = data.resp; resp(resp==-1) = 0;
        stim = data.stim; stim(stim==-1) = 0;
        
        hit  = length(find(resp(find(stim == 1))==1)) / length(resp(find(stim == 1)));
        fa   = length(find(resp(find(stim == 0))==1)) / length(resp(find(stim == 0)));
        
        grandavg.dprime(sj) = norminv(hit) - norminv(fa);
        
        % see if there is a difference in the actual stimulus as well
        grandavg.difficulty(sj, :) = mean(abs(data.motionstrength));
    end
end

subplot(441)
% plot the whole thing
plot(LogisticFit.x*100, nanmean(grandavg.yfit), 'color', [0.3 0.3 0.3]); % average function fit
hold on;

% then the datapoints, with standard deviation
p = ploterr(100*squeeze(nanmean(grandavg.xpts)), squeeze(nanmean(grandavg.ypts)), ...
    100*squeeze(nanstd(grandavg.xpts)), ...
    squeeze(nanstd(grandavg.ypts)) , 'k.', 'hhxy', 0.00000000000000000000000001);
p(1).MarkerSize = 8; p(1).Color = 'w'
box off; axis tight; set(gca, 'ytick', [0 0.5 1], 'xtick', [-30 -15 0 15 30]);
axis square;
xlabel('delta motion coherence'); ylabel('% response ''stronger''');
set(gca, 'xcolor', 'k', 'ycolor', 'k');
print(gcf, '-dpdf', sprintf('%s/fig1a_psychometricFunc.pdf', figpath));

