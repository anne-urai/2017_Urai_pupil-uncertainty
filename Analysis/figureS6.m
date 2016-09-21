
% reproduces
clearvars -except mypath;
global mypath;
close all; 

%% use nice shades of red and green
cols = cbrewer('qual', 'Set1', 8);
cols = cols([1 2], :);
nbins = 3;
mods = {'pupil', 'rt'};
cors = [1 0];
cnt = 1;
subjects = 1:27;

for m = 1:length(mods),
    for c = 1:length(cors),
        grandavg = postPupilBehaviour(mods{m}, nbins, cors(c));
        subplot(4,4,cnt); cnt = cnt + 1;
        
        % Use HOLD and ERRORBAR, passing axes handles to the functions.
        colors = cbrewer('qual', 'Set1', 8);
        switch cors(c)
            case 1
                thismarker = '.';
                thismarkersize = 14;
                thiscolor = colors(2,:);
            case 0
                thismarker = '^';
                thismarkersize = 4;
                thiscolor = colors(1,:);
        end
        
        % line to indicate 50% repetition
        y = grandavg.repetition;
        plot([1 nbins], [0.5 0.5], 'k:', 'linewidth', 0.5); hold on;
        
        % errorbar
        h = ploterr(1:nbins, nanmean(y), [], nanstd(y) ./sqrt(27), 'k-',  'abshhxy', 0);
        set(h(1), 'color', thiscolor, 'markersize', thismarkersize, 'marker',thismarker);
        
        % only for error markers, triangle
        if cors(c) == 0, set(h(1), 'markerfacecolor', 'w', 'markeredgecolor', thiscolor); end
        set(h(2), 'color', thiscolor); % line color
        
        xticklabs       = repmat({' '}, 1, nbins);
        xticklabs{1}    = 'low';
        xticklabs{end}  = 'high';
        if nbins == 3, xticklabs{2} = 'med'; end
        
        set(gca, 'xlim', [0.5 nbins+0.5], 'xtick', 1:nbins,  'xticklabel', xticklabs, ...
            'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5, 'box', 'off', 'xminortick', 'off', 'yminortick', 'off');
        axis square; 
        
        % determine y label and limits
        set(gca, 'ylim', [0.46 0.58], 'ytick', 0.48:0.05:0.58);
        ylabel('P(repeat)');
        xlim([0.5 nbins+0.5]);
        
        % do Bayesian ANOVA to get Bayes Factors
        statdat             = table;
        statdat.DV          = y(:);
        sj                  = repmat(1:27, nbins, 1)';
        statdat.subjnr      = sj(:);
        ft                  = repmat(transpose(1:nbins), 1, 27)';
        statdat.prevPupilBins = ft(:);
        writetable(statdat, sprintf('%s/Data/CSV/ANOVAdat.csv', mypath));
        system('/Library/Frameworks/R.framework/Resources/bin/R < BayesFactorANOVA.R --no-save');
        statres{cnt} = readtable(sprintf('%s/Data/CSV/ANOVAresults.csv', mypath)); % fetch results
        
        yval = max(get(gca, 'ylim'));
        if statres{cnt}.pvalue < 0.05, % only show if significant
            mysigstar(gca, [1 nbins], [yval yval], statres{cnt}.pvalue, 'k', 'down');
        end
        fprintf('\n %s, correct %d, ANOVA F(%d,%d) = %.3f, p = %.3f, bf10 = %.3f \n', ...
            mods{m}, cors(c), statres{cnt}.df1, statres{cnt}.df2, statres{cnt}.F, statres{cnt}.pvalue, statres{cnt}.bf10);
        
        switch mods{m}
            case 'pupil'
                xlabel('Previous trial pupil');
            case 'rt'
                xlabel('Previous trial RT');
            otherwise
                xlabel(sprintf('Previous trial %s', mods{m}));
        end
    end
end

%% also show the bar graphs of the error and correct fruend-derived model
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));

colors = cbrewer('qual', 'Set1', 8);
subplot(4,4,5); plotBetas([dat.correct(:, 1) dat.incorrect(:, 1)], colors([2 1], :));
set(gca, 'xtick', 1:2, 'xticklabel', {'Correct', 'Error'}, 'xticklabelrotation', -30); %ylim([-0.35 0.3]);
axis square;

colors = cbrewer('qual', 'Set1', 8);
subplot(4,4,6); plotBetas([dat.correct_pupil(:, 1) dat.incorrect_pupil(:, 1)], colors([2 1], :));
set(gca, 'xtick', 1:2, 'xticklabel', {'Pupil x correct', 'Pupil x error'}, 'xticklabelrotation', -30); %ylim([-0.35 0.3]);
axis square;

subplot(4,4,7); plotBetas([dat.correct_rt(:, 1) dat.incorrect_rt(:, 1)], colors([2 1], :));
set(gca, 'xtick', 1:2, 'xticklabel', {'RT x correct', 'RT x error'}, 'xticklabelrotation', -30); %ylim([-0.35 0.3]);
axis square;

print(gcf, '-dpdf', sprintf('%s/Figures/FigureS6.pdf', mypath));
