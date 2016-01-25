function [grandavg] = fig3de_Uncertainty_Accuracy(nbins)
% plots uncertainty by accuracy both for the modelfits and the pupil

global mypath;

if ~exist('nbins', 'var'), nbins = 20; end
subjects      = 1:27;
fitIndividual = false;

for sj = unique(subjects),
    
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));

    % divide into bins
    [ grandavg.pup(sj, :), grandavg.acc(sj, :), stdx, stdy] = ...
        divideintobins(data.decision_pupil, data.correct, nbins);
    
    % do a logistic regression
    x = [zscore(abs(data.motionstrength)) zscore(data.decision_pupil) ];
    x(data.trialnr < 2, :) = nan;
    
    [b] = glmfit(x, data.correct, ...
        'binomial','link','logit');
    
    % save for later
    grandavg.b(sj, :) = b;
    
    if fitIndividual,
        subplot(5,6,sj);
        yfit = glmval(b, x,'logit');
        hold on;
        plot(grandavg.pup(sj, :), grandavg.acc(sj,:), 'b.');
        axis tight; ylim([0.5 1]); xlims = get(gca, 'xlim');
        title(sprintf('P%02d', sj));
    end
    
end

% now, plot
grandavg.acc = 100 * grandavg.acc;
boundedline(1:nbins, nanmean(grandavg.acc), ...
    nanstd(grandavg.acc) / sqrt(length(subjects)), 'cmap', [0 0 0]);

xlabel('Pupil response');
ylabel({'Percent correct'});
axis tight; axis square;
set(gca, 'ytick', [60:5:80]);
set(gca, 'xtick', [1 nbins/2 nbins], 'xticklabel', {'low', 'medium', 'high'});
%offsetAxes;
box off; ylim([65 80]);

% do stats on these logistic betas
[h, p, ci, stats] = ttest(grandavg.b(:, end));
text(2, 69, sprintf('b = %.3f', mean(grandavg.b(:, end))));
text(2, 67, sprintf('p < %.3f', p));
set(gca, 'xcolor', 'k', 'ycolor', 'k');

%% fit logistic slope on high vs low pupil bins
subjects = 1:27;
grandavg.betas = nan(length(subjects), 2, 2);

for sj = unique(subjects),
    
    disp(sj);
    % get all the data
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', sj));
    
    % divide into low and high pupil bins
    puptrls{1} = find(data.decision_pupil < quantile(data.decision_pupil, 0.5));
    puptrls{2} = find(data.decision_pupil > quantile(data.decision_pupil, 0.5));
    
    % clear puptrls; % puptrls{1} = 1:height(data);
    
    for b = 1:length(puptrls),
        thisdat = data(puptrls{b}, :);
        
        % fit on all the trials, dont average beforehand!
        x = thisdat.motionstrength;
        y = thisdat.resp; y(y==-1) = 0;
        
        % sort both
        [x, idx] = sort(x);
        y = y(idx);
        
        [beta,dev,stats] = glmfit(x, [y ones(size(y))], ...
            'binomial','link','logit');
        newx = -3:0.01:3;
        yfit = glmval(beta, newx,'logit');
        % plot(newx, yfit,'-', 'color', cols(b, :));
        
        % also save the curve itself
        grandavg.betas(sj, b, :) = beta;
        grandavg.yfit(sj, b, :) = yfit;
        
        % also get the individual d prime
        % use 0 (for weaker motion answer) instead of -1
        resp = thisdat.resp; resp(resp==-1) = 0;
        stim = thisdat.stim; stim(stim==-1) = 0;
        
        hit  = length(find(resp(find(stim == 1))==1)) / length(resp(find(stim == 1)));
        fa   = length(find(resp(find(stim == 0))==1)) / length(resp(find(stim == 0)));
        
        grandavg.dprime(sj, b) = norminv(hit) - norminv(fa);
        
        % see if there is a difference in the actual stimulus as well
        grandavg.difficulty(sj, b) = mean(abs(thisdat.motionstrength));
    end
end

% make sure I actually plot the bar graph here to reproduce figure 3!
subplot(4,8,16);
hold on;
bar(1, mean(grandavg.betas(:, 1, 2)), 'facecolor', [0.6 0.6 0.6], 'edgecolor', 'none', 'barwidth', 0.4);
bar(2, mean(grandavg.betas(:, 2, 2)), 'facecolor', [0.4 0.4 0.4], 'edgecolor', 'none', 'barwidth', 0.4);
errorbar(1:2, squeeze(nanmean(grandavg.betas(:, :, 2))),  ...
    squeeze(std(grandavg.betas(:, :, 2))) / sqrt(length(subjects)), '.k');

% do ttest on regression coefficients
[~, pval(3), ~, stat] = ttest(grandavg.betas(:, 1, 1), grandavg.betas(:, 1, 2));
xlim([0.5 2.5]); box off;
mysigstar([1 2], [1 1], pval(3)); 
ylim([0.4 1.1]); 
set(gca, 'ytick', [0.5 1]);
xlabel('Pupil response'); set(gca, 'xtick', 1:2, 'xticklabel', {'low', 'high'});
ylabel('Current trial slope');
set(gca, 'xcolor', 'k', 'ycolor', 'k');

end
