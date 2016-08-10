% reproduces
global mypath;
close; figure;

correctness = []; % empty; both correct and error trials will be used
nbins = 3;

subplot(441); psychFuncShift_Bias_byResp('pupil', nbins, correctness);
subplot(442); psychFuncShift_Bias_Slope('pupil', nbins, correctness);

subplot(443); psychFuncShift_Bias_byResp('rt', nbins, correctness);
subplot(444); psychFuncShift_Bias_Slope('rt', nbins, correctness);

print(gcf, '-dpdf', sprintf('%s/Figures/figure5.pdf', mypath));

%% ========================================================= %
% new figure 5 (will become figure 4)
% ========================================================= %

clf;
correctness = []; % empty; both correct and error trials will be used
nbins = 3;
mods = {'pupil', 'rt'}; cnt = 0;
for m = 1:length(mods),
    
    % get all the measures as a function of previous trial pupil
    % how about post-error slowing and signed choice bias?
    
    grandavg = postPupilBehaviour(mods{m}, nbins, correctness);
    
    plotFields = {'sensitivity', 'RT', 'absoluteBias', 'repetition'};
    for s = 1:length(plotFields),
        
        cnt = cnt + 1;
        subplot(length(plotFields), length(plotFields), cnt);
        
        x   = 1:nbins;
        y   = grandavg.(plotFields{s});
        
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
        
        if strcmp(plotFields{s}, 'repetition'),
            % line to indicate 50 % repetition
            plot([1 nbins], [0.5 0.5], 'k:', 'linewidth', 0.5); hold on;
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
            'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5, 'box', 'off', 'xminortick', 'off', 'yminortick', 'on');
        axis square; xlim([0.5 nbins+0.5]);
        
        % determine y label and limits
        switch plotFields{s}
            case 'sensitivity'
                set(gca, 'ylim', [0.7 0.9]);
                ylabel('Slope');
            case 'RT'
                ylabel('Response time (s)');
                set(gca, 'ylim', [0.25 0.45]);
            case 'pes'
                ylabel('Post-error slowing (s)');
                set(gca, 'ylim', [-0.04 0.04]); % in s, so 40 ms
            case 'absoluteBias'
                set(gca, 'ylim', [0.2 0.4], 'ytick', [0.2 0.3 0.4]);
                ylabel('|Bias|');
            case 'repetition'
                set(gca, 'ylim', [0.48 0.55], 'ytick', [0.48:0.02:0.56]);
                ylabel('P(repeat)');
            otherwise
                ylabel(plotFields{s});
        end
        
        % do statistics, repeated measures anova across bins
        if size(grandavg.(plotFields{s}), 1) > 1,
            sj      = repmat(1:27, nbins, 1)';
            ft      = repmat(transpose(1:nbins), 1, 27)';
            f{1}    = ft(:);
            stats   = rm_anova(y(:), sj(:), f);
            yval    = max(get(gca, 'ylim'));
            if stats.f1.pvalue < 0.05, % only show if significant
                mysigstar(gca, [1 nbins], [yval yval], stats.f1.pvalue, 'k', 'down');
            end
        else
            axis tight; 
        end
        
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

% save
print(gcf, '-dpdf', sprintf('%s/Figures/figure4_zscoreRTperblock.pdf', mypath));


