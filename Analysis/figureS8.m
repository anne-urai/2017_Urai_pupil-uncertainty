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

clearvars -except mypath; close all; clc;
global mypath;
cd(sprintf('%s/Code/Analysis', mypath));

% it was not clear to me whether the reduced
% sequential dependencies for long RT trials
% could simply be due to the additional passage of time.

% 1. median time between responses
load(sprintf('%s/Data/GrandAverage/pupilgrandaverage.mat', mypath));
timebewteenResp = {};
for sj = 1:length(pupilgrandavg.timelock),
    respdiff = diff(pupilgrandavg.timelock{sj}(1).lock.trialinfo(:, 9)) ./ 100;
    trldiff = diff(pupilgrandavg.timelock{sj}(1).lock.trialinfo(:, 12));
    respdiff(trldiff ~= 1) = []; % only use the difference between subsequent trials
    timebetweenResp{sj} = respdiff;
end
timebetweenResp = cat(1, timebetweenResp{:});
median(timebetweenResp); % long-tailed distribution, so mean is biased

% preallocate
latencies = {'latency_fixation', 'latency_interval', 'latency_rtfeedback', ...
    'latency_feedbackisi', 'latency_total'};
nanzscore = @(x) (x - nanmean(x)) ./ nanstd(x);

if ~exist(sprintf('%s/Data/CSV/2ifc_data_allsj_withlatencies.csv', mypath), 'file'),
    %% compute and add latencies to csv file
    load(sprintf('%s/Data/GrandAverage/pupilgrandaverage.mat', mypath));
    
    % get the data that still include sample timings!
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
    
    latidx = [2 1; 6 2; 11 9; 19 11; 20 6];
    for l = 1:length(latencies),
        data.(latencies{l}) = nan(size(data.rt));
    end
    
    % compute all the latencies in each trial and across trials
    for sj = unique(data.subjnr)',
        
        trl         = pupilgrandavg.timelock{sj}(1).lock.trialinfo;
        trl(:, 19)  = circshift(trl(:, 1), -1);     % shift start of next fixation
        trl(:, 20)  = circshift(trl(:, 6), -1);     % shift start of next test stim
        badtrls     = find(diff(trl(:, 12)) ~= 1);  % only use trial with subsequent nrs
        
        % blocknrs for z-scoring
        blockchange = find(diff(trl(:, 12)) < 0);
        blocknrs = zeros(size(trl, 1), 1);
        for b = 1:length(blockchange)-1,
            blocknrs(blockchange(b)+1:blockchange(b+1)) = blocknrs(blockchange(b))+1;
        end
        blocknrs(blockchange(end)+1:end) = blocknrs(blockchange(end))+1;
        
        for l = 1:length(latencies),
            % latency between these events in seconds
            thislatency = (trl(:, latidx(l,1)) - trl(:, latidx(l,2))) ./ 100;
            
            if l > 3, % remove last trial when comparing to next one
                thislatency(badtrls) = NaN;
                thislatency(end) = NaN;
                assert(~any(thislatency < 0)); % if this happens something is wrong!
                data.(latencies{l})(data.subjnr == sj) = thislatency;
            else
                assert(~any(thislatency < 0)); % if this happens something is wrong!
                % normalize within each block, just as pupil and normRT
                blocks = unique(blocknrs)';
                thislatencyNorm = thislatency;
                for b = blocks,
                    thislatencyNorm(b==blocknrs) = nanzscore(thislatency(b==blocknrs));
                end
                data.(latencies{l})(data.subjnr == sj) = thislatencyNorm;
            end
            
        end
    end
    
    writetable(data, sprintf('%s/Data/CSV/2ifc_data_allsj_withlatencies.csv', mypath));
else
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj_withlatencies.csv', mypath));
end

% check
for l = 1:length(latencies),
    subplot(3,3,l);
    plot(data.(latencies{l}), '.');
    axisNotSoTight; box off;
    set(gca, 'xtick', find(diff(data.subjnr) ~= 0), 'xgrid', 'on', 'xticklabel', 1:27);
    title(latencies{l});
end

%% see if RTs correlate with overall latencies
for sj = unique(data.subjnr)',
    subplot(5,6,sj);
    plot(data.rt(data.subjnr == sj), data.latency_total(data.subjnr == sj), '.');
    [rho(sj), pval(sj)] = corr(data.rtNorm(data.subjnr == sj), ...
        data.latency_total(data.subjnr == sj) ./ 100, ...
        'type', 'spearman', 'rows', 'pairwise');
    title(sprintf('rho = %.3f, p = %.3f', rho(sj), pval(sj)));
    
    [ binnedex(sj, :), binnedy(sj, :)] = divideintobins(data.rtNorm(data.subjnr == sj), ...
        data.latency_total(data.subjnr == sj), 5, [], @nanmedian);
    axis square; axis tight; box off;
end

% now see if latencies between this trial and the next change as a function of RT
fprintf('mean rho: %.3f, range %.3f to %.3f, significant in %d out of 27 participants \n', ...
    mean(rho), min(rho), max(rho), length(find(pval < 0.05)));
close all;

%% next: check whether regressing out these latencies reduces the effect of RT
correctness = []; % empty; both correct and error trials will be used
nbins = 3;
mods = {'rt_withlatencies', 'latency_total'};
figure; cnt = 0;
for m = 1:length(mods),
    
    % get all the measures as a function of previous trial pupil
    % how about post-error slowing and signed choice bias?
    
    grandavg{m} = postPupilBehaviour(mods{m}, nbins, correctness);
    plotFields = {'repetition'};
    
    for s = 1:length(plotFields),
        
        cnt = cnt + 1;
        subplot(4,4,cnt);
        
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
                set(gca, 'ylim', [0.46 0.58], 'ytick', 0.48:0.05:0.58);
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
            case 'latency_total'
                xlabel('Current - previous trial latency');
            case 'rt'
                xlabel('Previous trial RT');
            case 'rt_withlatencies'
                xlabel('Previous trial RT');
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

%%
mods{1} = 'cleanpupil+rt';
if ~exist(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, mods{1}), 'file'),
    
    alldata = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj_withlatencies.csv', mypath));
    % corrplot(alldata, {'decision_pupil', 'latency_fixation', 'latency_interval', 'latency_rtfeedback', 'latency_feedbackisi'})
    cd(sprintf('%s/Code/serial-dependencies/data', mypath));
    %figure;
    for sj = unique(alldata.subjnr)',
        data = alldata(alldata.subjnr == sj, :);
        
        blockchange = find(diff(data.trialnr) < 0);
        blocknrs = zeros(size(data.rt));
        for b = 1:length(blockchange)-1,
            blocknrs(blockchange(b)+1:blockchange(b+1)) = blocknrs(blockchange(b))+1;
        end
        blocknrs(blockchange(end)+1:end) = blocknrs(blockchange(end))+1;
        
        % put all data in
        cleanpupil = projectout(data.decision_pupil, data.latency_fixation);
        cleanpupil = projectout(cleanpupil, data.latency_interval);
        cleanpupil = projectout(cleanpupil, data.latency_rtfeedback);
        assert(~any(isnan(cleanpupil)));
        
        newdat = [blocknrs data.sessionnr abs(data.motionstrength) (data.stim > 0) (data.resp > 0) ...
            nanzscore(cleanpupil) nanzscore(data.rtNorm)];
        dlmwrite(sprintf('2ifc_cleanpupil+rt_sj%02d.txt', sj), ...
            newdat,'delimiter','\t','precision',4);
    end
    
    % run model and retrieve data
    cd(sprintf('%s/Code/serial-dependencies', mypath));
    system([sprintf('for sj in {1..27}; do filename=$(printf "data/2ifc_%s_sj%%02d.txt" $sj);', mods{1}), ...
        sprintf('echo $filename; python2.7 analysis.py -fr -n10 -p "%s/Data/serialmodel/" $filename; sleep 5; done', mypath)]);
    cd(sprintf('%s/Code/Analysis', mypath));
    a6_retrieveDataFromPython(mods{1});
end

% print the regression weights
mods = {'cleanpupil+rt'};
for m = 1:length(mods),
    load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, mods{m}));
    subplot(4,4,m+2);
    plotBetas([dat.response_pupil(:, 1) ...
        dat.stimulus_pupil(:, 1)  dat.response_rt(:, 1) dat.stimulus_rt(:, 1)], 0.5);
    [~, pval, ~, stat] = ttest(dat.response_pupil(:, 1), dat.stimulus_pupil(:, 1));
    mysigstar(gca, [1 2], -0.08, pval);
    [~, pval, ~, stat] = ttest(dat.response_rt(:, 1), dat.stimulus_rt(:, 1));
    mysigstar(gca, [3 4], -0.08, pval);
    xlim([0.5 4.5]);
    set(gca, 'xtick', 1:4, 'xticklabel', ...
        {'Pupil x choice', 'Pupil x stimulus', 'RT x choice', 'RT x stimulus'}, ...
        'xticklabelrotation', -30);
end

alldata = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj_withlatencies.csv', mypath));
for sj = unique(alldata.subjnr)',
    [rho1(sj), pval1(sj)] = corr(alldata.decision_pupil(alldata.subjnr == sj), ...
        alldata.latency_interval(alldata.subjnr == sj), 'type', 'spearman');
    [rho2(sj), pval2(sj)] = corr(alldata.decision_pupil(alldata.subjnr == sj), ...
        alldata.latency_rtfeedback(alldata.subjnr == sj), 'type', 'spearman');
end

fprintf('mean Spearman''s rho %.3f, range %.3f to %.3f, significant in %d out of 27 observers \n', ...
    mean(rho1), min(rho1), max(rho1), length(find(pval1 < 0.05)));
fprintf('mean Spearman''s rho %.3f, range %.3f to %.3f, significant in %d out of 27 observers\n', ...
    mean(rho2), min(rho2), max(rho2), length(find(pval2 < 0.05)));

print(gcf, '-dpdf', sprintf('%s/Figures/FigureS8.pdf', mypath));

