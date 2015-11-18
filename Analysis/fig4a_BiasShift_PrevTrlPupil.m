%% split trials by their previous response and  
clear; clc; close all;

subjects = 1:27;
for sj = unique(subjects),
    
    data = readtable(sprintf('~/Data/UvA_pupil/CSV/2ifc_data_sj%02d.csv', sj));
    
    % outcome vector need to be 0 1 for logistic regression
    data.resp(data.resp == -1) = 0;
    % remove the trials at the end of each block
    endOfBlockTrls = find(data.trialnr == 50);
    
    % previous response
    resps = [0 1];
    for r = 1:2,
        
        trls = find(data.resp == resps(r) & data.correct == 1);
        
        % with this selection, take the trials after that
        % exclude trials at the end of the block
        thisdat = data(setdiff(trls, endOfBlockTrls)+1, :);
        
        [b, dev, stats] = glmfit(thisdat.motionstrength, thisdat.resp, ...
            'binomial','link','logit');
        x = -5:0.01:5;   yfit = glmval(b, x, 'logit');
        %  plot(x, yfit);
        
        % save betas
        grandavg.logistic(sj, r,  :) = b;
        grandavg.curve(sj, r,  :)    = yfit; % curve
        
    end % resp
end % sj


figure;
colors = linspecer(4); colors = colors([1 4], :);

subplot(4,4,1); hold on;
for i = 1:2,
    b = bar(i, mean(grandavg.logistic(:, i, 1)), ...
        'edgecolor', 'none', 'barwidth', 0.3, 'facecolor', colors(i, :));
end
h = ploterr(1:2, squeeze(mean(grandavg.logistic(:, :, 1))), ...
    [], squeeze(std(grandavg.logistic(:, :, 1))) / sqrt(length(subjects)), ...
    'k.', 'hhxy', 0.0000000000000001);
set(h(1), 'marker', 'none');
axis tight;
axis square;

ylabel({'Response bias'});
set(gca, 'ytick', [-.15 0  0.15], 'ylim', [-0.15 0.15]);
xlabel('Previous response');

[~, p ] = permtest(grandavg.logistic(:, 1, 1), grandavg.logistic(:, 2, 1))

set(gca, 'xtick', [1 2], 'xticklabel', {'s1 < s2', 's1 > s2'});
print(gcf, '-dpdf', '~/Dropbox/Figures/uncertainty/Fig4_valueShift_allTrials.pdf');


%% 2. dependence on pupil-based uncertainty, all trials
% does this pattern increase with higher uncertainty on previous trial?
clear; clc; 
nbins = 3;
subjects = 1:27;
for sj = unique(subjects),
    
    data = readtable(sprintf('~/Data/UvA_pupil/CSV/2ifc_data_sj%02d.csv', sj));
    
    % outcome vector need to be 0 1 for logistic regression
    data.resp(data.resp == -1) = 0;
    % remove the trials at the end of each block
    endOfBlockTrls = find(data.trialnr == 50);
    
    % previous response
    resps = [0 1];
    for r = 1:2,
        
        if nbins > 2,
            uncQs = quantile(data.decision_pupil(data.resp == resps(r)), nbins - 1);
        elseif nbins == 2,
            uncQs = median(data.decision_pupil(data.resp == resps(r)));
        end
        
        % uncertainty bins
        for u = 1:nbins,
            
            switch u
                case 1
                    trls = find(data.resp == resps(r) & data.decision_pupil < uncQs(u));
                case nbins
                    trls = find(data.resp == resps(r) & data.decision_pupil > uncQs(u-1));
                otherwise
                    trls = find(data.resp == resps(r) & ...
                        data.decision_pupil > uncQs(u-1) & data.decision_pupil < uncQs(u));
            end
            
            % only use correct trials!
            trls = intersect(trls, find(data.correct == 1));
            
            % with this selection, take the trials after that
            % exclude trials at the end of the block
            thisdat = data(setdiff(trls, endOfBlockTrls)+1, :);
            
            [b, dev, stats] = glmfit(thisdat.motionstrength, thisdat.resp, ...
                'binomial','link','logit');
            x = -5:0.01:5;   yfit = glmval(b, x, 'logit');
            %  plot(x, yfit);
            
            % save betas
            grandavg.logistic(sj, r, u, :) = b;
            grandavg.curve(sj, r, u, :)    = yfit; % curve
            
        end % uncertainty bin
    end % resp
end % sj


cnt = 3;
colors = linspecer(4); colors = colors([1 4], :);
x = 1:u;

subplot(4,4,cnt); hold on;
for r = 1:2,
    % b = bar(x, squeeze(nanmean(grandavg.logistic(:, r, :, 1))), ...
    %     'FaceColor', colors(r, :), 'EdgeColor', 'none', 'BarWidth', 0.4);
    h = ploterr(x, ...
        squeeze(nanmean(grandavg.logistic(:,r, :, 1))),  [], ...
        squeeze(nanstd(grandavg.logistic(:, r, :, 1))) ./ sqrt(length(subjects)), ...
        '-',  'hhxy', 0.001);
    set(h(1), 'color', colors(r, :), ...
        'marker', '.', 'markerfacecolor', colors(r, :), 'markersize', 12);
    set(h(2), 'color', colors(r, :));
    x = x + 0.2;
end
plot([1 max(x)], [0 0], '-', 'color', [0.4 0.4 0.4], 'LineWidth', 0.1);

ylabel({'Response bias'});
set(gca, 'ytick', [-.15 0  0.15]);
xlabel('Previous trial pupil');

axis square;
%set(gca, 'xtick', [1 2.1 3.2], 'xticklabel', {'low',  'medium', 'high'});
%box off; offsetAxes;
%set(gca, 'xtick', [1.1 2.1 3.1], 'xticklabel', {'low', 'medium', 'high'});
print(gcf, '-dpdf', '~/Dropbox/Figures/uncertainty/Fig4_valueShift_allTrials.pdf');


%% same plot, but then split by pupil * response  (since the SJ dont have access to the difficulty level, only their internal state
