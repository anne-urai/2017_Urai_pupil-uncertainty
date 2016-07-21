function [] = MedianSplit(whichmodulator, field)
% plot correlation between subjects

global mypath;
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));
mycolmap = cbrewer('div', 'PuOr', 3);

% split between repeaters and alternators
mat = dat.([field '_' whichmodulator])(:, 1);

for i = 1:2,
    
    switch i
        case 1
            subjects{i} = find(dat.response(:, 1) < 0);
            colors = mycolmap(1,:);
        case 2
            subjects{i} = find(dat.response(:, 1) > 0);
            colors = mycolmap(3,:);
            
    end
    plotBetasSwarmUnpaired(i, mat(subjects{i}, :), colors)
end

% do stats between them
xlim([0.5 2.5]);
set(gca, 'xtick', 1:2, 'xticklabel', {'alternators', 'repeaters'}, ...
   'xticklabelrotation', -30);

ylims = get(gca, 'ylim');
yrange = range(ylims);
ylim([ylims(1) - yrange*0.2 ylims(2) + yrange*0.1]);

[h, p, ci, stats] = ttest2(mat(subjects{1}), mat(subjects{2}));
mysigstar(gca, [1 2], min(get(gca, 'ylim')), p);

end

function [] = plotBetasSwarmUnpaired(i, beta, colors)
% using the timewindow that is indicated in the regression timecourse plot,
% show the shape of the pupil vs motionstrength pattern

if ~exist('colors', 'var'); colors = cbrewer('qual', 'Set1', 8); 
colors = [0 0 0; colors]; end % start with black
hold on; % paired

bar(i, squeeze(mean(beta)), 'edgecolor', 'none', ...
    'facecolor', [0.8 0.8 0.8], 'barwidth', 0.5);

% scatter all the points
scatter(i * ones(1, size(beta, 1)), beta, ...
    10, colors, 'o', 'linewidth', 0.5, 'jitter', 'on', 'jitteramount', 0);
set(gca, 'xtick', [1 2], 'xminortick', 'off');
ylabel('Beta weights (a.u.)'); 

xlim([0.5 size(beta,2) + 0.5]);
axis tight;
box off;
[~, pval, ~, stat] = ttest(beta, 0, 'tail', 'both');
mysigstar(gca, i, min(get(gca, 'ylim')), pval);

end

