% figure 4 overview
clearvars -except mypath; clc; close all;
global mypath;
mods = {'fbpupil', 'fb+decpupil'};
nbins = 3;
close; figure;

for m = 1:length(mods),
    
    grandavg{m} = postPupilBehaviour(mods{m}, nbins, []);
    subplot(4,4,(m-1)*4+1);
    
    x   = 1:nbins;
    y   = grandavg{m}.repetition;
    
    thismarker = '.';
    thismarkersize = 14;
    thiscolor = [0 0 0];
    
    % line to indicate 50 % repetition
    plot([1 nbins], [0.5 0.5], 'k:', 'linewidth', 0.5); hold on;
    h = ploterr(x, nanmean(y), [], nanstd(y) ./sqrt(27), 'k-',  'abshhxy', 0);
    set(h(1), 'color', thiscolor, 'markersize', thismarkersize, 'marker',thismarker);
    set(h(2), 'color', thiscolor); % line color
    
    xticklabs       = repmat({' '}, 1, nbins);
    xticklabs{1}    = 'low';
    xticklabs{end}  = 'high';
    if nbins == 3, xticklabs{2} = 'med'; end
    set(gca, 'xlim', [0.5 nbins+0.5], 'xtick', 1:nbins,  'xticklabel', xticklabs, ...
        'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5, 'box', 'off', 'xminortick', 'off', 'yminortick', 'off');
    axis square; xlim([0.5 nbins+0.5]);
    
    % determine y label and limits
    set(gca, 'ylim', [0.48 0.54], 'ytick', [0.48:0.02:0.54]);
    ylabel('P(repeat)');
    
    % do statistics, repeated measures anova across bins
    sj      = repmat(1:27, nbins, 1)';
    ft      = repmat(transpose(1:nbins), 1, 27)';
    f{1}    = ft(:);
    stats   = rm_anova(y(:), sj(:), f); % from Valentin
    
    % do Bayesian ANOVA to get Bayes Factors
    statdat         = table;
    statdat.DV      = y(:);
    statdat.subjnr  = sj(:);
    statdat.prevPupilBins = ft(:);
    writetable(statdat, sprintf('%s/Data/CSV/ANOVAdat.csv', mypath));
    system('/Library/Frameworks/R.framework/Resources/bin/R < BayesFactorANOVA.R --no-save');
    statres = readtable(sprintf('%s/Data/CSV/ANOVAresults.csv', mypath)); % fetch results
    
    yval    = max(get(gca, 'ylim'));
    if stats.f1.pvalue < 0.05, % only show if significant
        mysigstar(gca, [1 nbins], [yval yval], statres.pvalue, 'k', 'down');
    end
    fprintf('\n %s, ANOVA F(%d,%d) = %.3f, p = %.3f, bf10 = %.3f \n', ...
        mods{m}, statres.df1, statres.df2, statres.F, statres.pvalue, statres.bf10);
    
    switch mods{m}
        case 'pupil'
            xlabel('Previous trial pupil');
        case 'rt'
            xlabel('Previous trial RT');
        otherwise
            xlabel(sprintf('Previous trial %s', mods{m}));
    end
    
    subplot(4,8,(m-1)*8+3);
    if m < 3,
        load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, mods{m}));
        plotBetas([dat.response_pupil(:, 1) ...
            dat.stimulus_pupil(:, 1)], ...
            0.5*ones(3,3));
        set(gca, 'xtick', 1:2, 'xticklabel', []);
    else
        load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, mods{m-1}));
        plotBetas([dat.response_rt(:, 1) ...
            dat.stimulus_rt(:, 1)], ...
            0.5*ones(3,3));
        set(gca, 'xtick', 1:2, 'xticklabel', ...
            {'Pupil x choice', 'Pupil x stimulus'}, ...
            'xticklabelrotation', -30);
    end
    ylim([-0.09 0.06]);
end
set(gca, 'xtick', 1:2, 'xticklabel', ...
    {'Choice', 'Stimulus'}, ...
    'xticklabelrotation', -30);

print(gcf, '-dpdf', sprintf('%s/Figures/FigureS9.pdf', mypath));
