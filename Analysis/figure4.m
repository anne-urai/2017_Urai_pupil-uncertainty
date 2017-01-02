% This code reproduces the analyses in the paper
% Urai AE, Braun A, Donner THD (2016) Pupil-linked arousal is driven
% by decision uncertainty and alters serial choice bias.
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% If you use the Software for your own research, cite the paper.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%
% Anne Urai, 2016
% anne.urai@gmail.com

clearvars -except mypath;
clc;
close all;
correctness = []; % empty; both correct and error trials will be used
nbins = 3;
mods = {'pupil', 'rt'}; cnt = 0;
for m = 1:length(mods),
    
    % get all the measures as a function of previous trial pupil
    grandavg{m} = postPupilBehaviour(mods{m}, nbins, correctness);
    disp(mods{m});
    
    plotFields = {'repetition', 'lapse', 'absoluteBias', ...
        'sensitivity', 'lapse','pesRegressedout'};
    
    for s = 1:length(plotFields),
        
        cnt = cnt + 1;
        subplot(length(plotFields), length(plotFields), cnt);
        
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
        
        switch mods{m}
            case 'pupil'
                xticklabs       = repmat({' '}, 1, nbins);
                xticklabs{1}    = 'low';
                xticklabs{end}  = 'high';
                % if nbins == 3, xticklabs{2} = 'med'; end
            case 'rt'
                xticklabs       = repmat({' '}, 1, nbins);
                xticklabs{1}    = 'fast';
                xticklabs{end}  = 'slow';
                % if nbins == 3, xticklabs{2} = 'med'; end
        end
        
        set(gca, 'xlim', [0.5 nbins+0.5], 'xtick', 1:nbins,  'xticklabel', xticklabs, ...
            'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5, 'box', 'off', 'xminortick', 'off', 'yminortick', 'off');
        axis square; xlim([0.5 nbins+0.5]);
        
        % determine y label and limits
        switch plotFields{s}
            case 'sensitivity'
                set(gca, 'ylim', [0.75 0.95], 'ytick', 0.75:0.1:1);
                ylabel('Slope');
            case 'RT'
                ylabel('Response time (s)');
                set(gca, 'ylim', [0.25 0.45]);
            case {'pes', 'pesMatched', 'pesRegressedout'}
                ylabel('Post-error slowing (s)');
                set(gca, 'ylim', [-0.02 0.02], 'ytick', [-0.02 0 0.02]); % in s, so 40 ms
            case 'absoluteBias'
                set(gca, 'ylim', [0.2 0.4], 'ytick', [0.2:0.1:0.6]);
                ylabel('Absolute bias');
            case 'repetition'
                set(gca, 'ylim', [0.48 0.54], 'ytick', [0.48:0.03:0.54]);
                ylabel('P(repeat)');
            case 'lapse'
                set(gca, 'ylim', [0 0.02], 'ytick', [0:0.01:0.02]);
                ylabel('P(lapse)');
            otherwise
                ylabel(plotFields{s});
        end
        
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
        
        % point to your version of R here
        [status, cmdout] = system('/Library/Frameworks/R.framework/Resources/bin/R < BayesFactorANOVA.R --no-save');
        statres{s} = readtable(sprintf('%s/Data/CSV/ANOVAresults.csv', mypath)); % fetch results
        
        yval    = max(get(gca, 'ylim'));
        if stats.f1.pvalue < 0.05, % only show stars if significant
            mysigstar(gca, [1 nbins], [yval yval], statres{s}.pvalue, 'k', 'down');
        else
            ylims = get(gca, 'ylim');
            text(1, ylims(2), sprintf('Bf_{10} = %.3f', statres{s}.bf10));
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
    
    % ========================================================= %
    % just do stats
    % ========================================================= %
    
    statsFields = {'sensitivity', 'signedBias', ...
        'absoluteBias','pesRegressedout', 'repetition', 'RT', 'lapse'};
    
    for s = 1:length(statsFields),
        y       = grandavg{m}.(statsFields{s});
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
        [status, cmdout] = system('/Library/Frameworks/R.framework/Resources/bin/R < BayesFactorANOVA.R --no-save');
        statres{s} = readtable(sprintf('%s/Data/CSV/ANOVAresults.csv', mypath)); % fetch results
    end
    
    % print results
    for s = 1:length(statres),
        fprintf('\n %s, ANOVA F(%d,%d) = %.3f, p = %.3f, bf10 = %.3f \n', ...
            statsFields{s}, statres{s}.df1, statres{s}.df2, statres{s}.F, statres{s}.pvalue, statres{s}.bf10);
    end
    cnt = cnt + length(plotFields);
end

% save
print(gcf, '-dpdf', sprintf('%s/Figures/Figure4.pdf', mypath));
