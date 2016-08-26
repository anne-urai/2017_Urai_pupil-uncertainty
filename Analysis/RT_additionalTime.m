
clearvars -except mypath; close all; clc;
global mypath;

% it was not clear to me whether the reduced
% sequential dependencies for long RT trials
% could simply be due to the additional passage of time.

load(sprintf('%s/Data/GrandAverage/pupilgrandaverage.mat', mypath));

% get the data that still include sample timings!
data = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
data.latencies = nan(size(data.rt));

for sj = unique(data.subjnr)',
    trl       = pupilgrandavg.timelock{sj}(1).lock.trialinfo;
    latencies = circshift(trl(:, 6), -1) - trl(:, 6); % time between this stimulus and the next
    latencies(find(diff(trl(:, 12) ~= 1))) = NaN;
    latencies(latencies < 0)        = NaN;
    lat = latencies(~isnan(latencies));
    
    % little function that properly rank transforms
    ranked = ranktransform(lat);

    % put back into the array with nans
    latencies(~isnan(latencies)) = ranked;
    data.latencies(data.subjnr == sj) = latencies;
end

writetable(data, sprintf('%s/Data/CSV/2ifc_data_allsj_withlatencies.csv', mypath));

%% get data
data = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj_withlatencies.csv', mypath));

for sj = unique(data.subjnr)',
    subplot(5,6,sj);
    plot(data.rt(data.subjnr == sj), data.latencies(data.subjnr == sj), '.');
    [rho(sj), pval(sj)] = corr(data.rtNorm(data.subjnr == sj), ...
        data.latencies(data.subjnr == sj) ./ 100, ...
        'type', 'spearman', 'rows', 'pairwise');
    title(sprintf('rho = %.3f, p = %.3f', rho(sj), pval(sj)));
    
    [ binnedex(sj, :), binnedy(sj, :)] = divideintobins(data.rtNorm(data.subjnr == sj), ...
        data.latencies(data.subjnr == sj), 5, [], @nanmedian);
    axis square; axis tight; box off;
end

% now see if latencies between this trial and the next change as a function of RT
fprintf('mean rho: %.3f, range %.3f to %.3f, significant in %d out of 27 participants \n', ...
    mean(rho), min(rho), max(rho), length(find(pval < 0.05)));

figure;
subplot(331);
errorbar(nanmean(binnedex), nanmean(binnedy), nanstd(binnedy) ./ sqrt(27), 'linewidth', 1);
xlabel('Response time (z)');
ylabel('Time between stim (rank)');
axisNotSoTight;
box off; axis square;

%% next: check whether regressing out these latencies reduces the effect of RT

correctness = []; % empty; both correct and error trials will be used
nbins = 3;
mods = {'rt_withlatencies', 'latencies'}; cnt = 1;
for m = 1:length(mods),
    
    % get all the measures as a function of previous trial pupil
    % how about post-error slowing and signed choice bias?
    
    grandavg{m} = postPupilBehaviour(mods{m}, nbins, correctness);
    plotFields = {'repetition'};
    
    for s = 1:length(plotFields),
        
        cnt = cnt + 1;
        subplot(3,3,cnt);
        
        x   = 1:nbins;
        y   = grandavg{m}.(plotFields{s});
        
        % Use HOLD and ERRORBAR, passing axes handles to the functions.
        colors = cbrewer('qual', 'Set1', 8);
        if isempty(correctness),
            thismarker = '.';
            thismarkersize = 14;
            thiscolor = [0 0 0];
        else
            switch correctness
                case 1
                    thismarker = '.';
                    thismarkersize = 14;
                    thiscolor = colors(2,:);
                case 0
                    thismarker = '^';
                    thismarkersize = 4;
                    thiscolor = colors(1,:);
            end
        end
        
        switch plotFields{s},
            case 'repetition'
                % line to indicate 50 % repetition
                plot([1 nbins], [0.5 0.5], 'k:', 'linewidth', 0.5); hold on;
            case {'pes', 'pesMatched', 'pesRegressedout'}
                plot([1 nbins], [0.0 0.0], 'k:', 'linewidth', 0.5); hold on;
        end
        
        % errorbar
        h = ploterr(x, nanmean(y), [], nanstd(y) ./sqrt(27), 'k-',  'abshhxy', 0);
        set(h(1), 'color', thiscolor, 'markersize', thismarkersize, 'marker',thismarker);
        
        % only for error markers, triangle
        if correctness == 0, set(h(1), 'markerfacecolor', 'w', 'markeredgecolor', thiscolor); end
        set(h(2), 'color', thiscolor); % line color
        
        xticklabs       = repmat({' '}, 1, nbins);
        xticklabs{1}    = 'low';
        xticklabs{end}  = 'high';
        if nbins == 3, xticklabs{2} = 'med'; end
        
        set(gca, 'xlim', [0.5 nbins+0.5], 'xtick', 1:nbins,  'xticklabel', xticklabs, ...
            'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5, 'box', 'off', 'xminortick', 'off', 'yminortick', 'off');
        axis square; xlim([0.5 nbins+0.5]);
        
        % determine y label and limits
        switch plotFields{s}
            case 'sensitivity'
                set(gca, 'ylim', [0.7 0.9], 'ytick', 0.7:0.1:0.9);
                ylabel('Slope');
            case 'RT'
                ylabel('Response time (s)');
                set(gca, 'ylim', [0.25 0.45]);
            case {'pes', 'pesMatched', 'pesRegressedout'}
                ylabel('Post-error slowing (s)');
                % ylabel(plotFields{s});
                set(gca, 'ylim', [-0.02 0.02], 'ytick', [-0.02 0 0.02]); % in s, so 40 ms
            case 'absoluteBias'
                set(gca, 'ylim', [0.2 0.4], 'ytick', [0.2 0.3 0.4]);
                ylabel('Absolute bias');
            case 'repetition'
                set(gca, 'ylim', [0.48 0.54], 'ytick', 0.48:0.02:0.54);
                ylabel('P(repeat)');
            otherwise
                ylabel(plotFields{s});
        end
        
        % do statistics, repeated measures anova across bins
        if size(grandavg{m}.(plotFields{s}), 1) > 1,
            
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
            statres{m} = readtable(sprintf('%s/Data/CSV/ANOVAresults.csv', mypath)); % fetch results
            
            yval    = max(get(gca, 'ylim'));
            if stats.f1.pvalue < 0.05, % only show if significant
                mysigstar(gca, [1 nbins], [yval yval], statres{m}.pvalue, 'k', 'down');
            end
        else
            axis tight;
        end
        
        switch mods{m}
            case 'pupil'
                xlabel('Previous trial pupil');
            case 'rt'
                xlabel('Previous trial RT');
            case 'rt_withlatencies'
                xlabel({'Previous trial RT'; 'latencies removed'});
            otherwise
                xlabel(sprintf('Previous trial %s', mods{m}));
        end
    end
end

for s = 1:2,
    fprintf('\n %s, ANOVA F(%d,%d) = %.3f, p = %.3f, bf10 = %.3f \n', ...
        mods{s}, statres{m}.df1, statres{s}.df2, statres{s}.F, statres{s}.pvalue, statres{s}.bf10);
end
% save
print(gcf, '-dpdf', sprintf('%s/Figures/RTadditionalTime.pdf', mypath));

