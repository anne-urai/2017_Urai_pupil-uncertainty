%% Weibull curves
clearvars -except mypath
close all; clc;
nbins    = 3;
grandavg = postPupilBehaviour('pupil', nbins, []);

% get individual choice weights
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));
[w, idx] = sort(dat.response(:, 1));

if 0,
    %% visually inspect
    figure;
    colors = cbrewer('qual', 'Set1', nbins); %colors = colors(2:end, :);
    for sj = 1:27,
        subplot(5,6,sj); hold on;
        for b = 1:nbins,
            p(b) = plot(linspace(0, 6, 100), squeeze(grandavg.weibullCurve(idx(sj), b, :)), 'color', colors(b,:));
            plot(squeeze(grandavg.weibullPtsX(idx(sj), b, :))', squeeze(grandavg.weibullPtsY(idx(sj), b, :))',...
                'o', 'color', colors(b,:), 'markersize', 2);
        end
        xlim([0 7]); ylim([0.45 1]); box off; axis square;
        title(sprintf('P%02d', idx(sj)));
    end
    l  = legend(p, {'low', 'medium', 'high'}); legend boxoff;
    sp = subplot(5,6,sj+1); spos = get(gca, 'position');
    set(l, 'position', spos); axis off;
    suplabel('Subjects sorted from alternators to repeaters', 't');
    suplabel('Evidence strength', 'x');
    suplabel('Accuracy', 'y');
    print(gcf, '-dpdf', sprintf('%s/Figures/WeibullByLastTrialPupil.pdf', mypath));
end

%% overview for in rebuttal letter
clearvars -except dat grandavg nbins mypath; clf;
groups{1} = find(dat.response(:, 1) < 0); % alternators
groups{2} = find(dat.response(:, 1) > 0); % repeaters

mycolmap = cbrewer('div', 'PuOr', 3);
colors{1} = mycolmap(1,:);
colors{2} = mycolmap(3,:);

for g = 1:length(groups),
    cnt = 0;
    plotFields = {'accuracy', 'weibull'};
    
    for s = 1:length(plotFields),
        cnt = cnt + 1;
        if length(plotFields) > 4,
            nrsubplots = length(plotFields);
        else nrsubplots = 4; 
        end
        subplot(nrsubplots, nrsubplots, cnt);
        hold on;
        
        x   = [1:nbins] + (g-1)*0.1;
        y   = grandavg.(plotFields{s})(groups{g}, :); % only those subjects
        thiscolor = colors{g};
        
        % errorbar
        h = ploterr(x, nanmean(y), [], nanstd(y) ./sqrt(length(y)), 'k-',  'abshhxy', 0);
        set(h(1), 'color', thiscolor, 'markersize', 10, 'marker', '.');
        set(h(2), 'color', thiscolor); % line color
        
        % layout
        xticklabs       = repmat({' '}, 1, nbins);
        xticklabs{1}    = 'low';
        xticklabs{end}  = 'high';
        if nbins == 3, xticklabs{2} = 'med'; end
        set(gca, 'xlim', [0.5 nbins+0.5], 'xtick', 1:nbins,  'xticklabel', xticklabs, ...
            'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5, 'box', 'off', 'xminortick', 'off', 'yminortick', 'off');
        axis square; xlim([0.5 nbins+0.5]);
        
        % determine y label and limits
        switch plotFields{s}
            case 'accuracy'
                ylabel('Accuracy (%)');
                ylim([0.7 0.78]);
            case 'weibull'
                ylabel('Threshold (a.u.)');
                ylim([1.5 2.4]);
            otherwise
                ylabel(plotFields{s});
        end
        xlabel('Previous trial pupil');
        
        if 0,
            % do statistics, repeated measures anova across bins
            if size(y, 1) > 1,
                
                sj      = repmat(1:size(y, 1), nbins, 1)';
                ft      = repmat(transpose(1:nbins), 1, size(y, 1))';
                f{1}    = ft(:);
                stats   = rm_anova(y(:), sj(:), f); % from Valentin
                
                %                 % do Bayesian ANOVA to get Bayes Factors
                %                 statdat         = table;
                %                 statdat.DV      = y(:);
                %                 statdat.subjnr  = sj(:);
                %                 statdat.prevPupilBins = ft(:);
                %                 writetable(statdat, sprintf('%s/Data/CSV/ANOVAdat.csv', mypath));
                %                 system('/Library/Frameworks/R.framework/Resources/bin/R < BayesFactorANOVA.R --no-save');
                %                 statres{s} = readtable(sprintf('%s/Data/CSV/ANOVAresults.csv', mypath)); % fetch results
                %
                if stats.f1.pvalue < 0.05, % only show if significant
                    yval = max(get(gca, 'ylim'));
                    mysigstar(gca, [1 nbins], [yval yval], stats.f1.pvalue, 'k', 'down');
                end
            else
                axis tight;
            end
        end
        
    end
end
print(gcf, '-dpdf', sprintf('%s/Figures/PostPupilAccuracy_perGroup.pdf', mypath));

