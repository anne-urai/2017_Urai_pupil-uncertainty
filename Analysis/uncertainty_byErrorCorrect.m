function [b, bint] = uncertainty_byErrorCorrect(field, nbins)
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

global mypath;

% get all data
alldata = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
warning('error', 'stats:LinearModel:RankDefDesignMat'); % stop if this happens
subjects = unique(alldata.subjnr)'; 
if ~exist('nbins', 'var'); nbins = 6; end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE OVERVIEW OF THE PUPIL UNCERTAINTY CORRELATION FOR ALL THESE DIFFERENT FIELDS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grandavg.(field).data = nan(length(subjects), 2, nbins);

for sj = subjects,
    % get data from this subject
    data = alldata(alldata.subjnr == sj, :);
    
    % if looking at the baseline, get the one from the next trial
    switch field
        case 'baseline_pupil'
            data.baseline_pupil = circshift(zscore(data.baseline_pupil), -1);
    end
    
    % normalization etc
    data.motionstrength = (abs(data.motionstrength));
    
    % loop over error and correct
    cors = [0 1];
    for corr = cors,
        
        % RATHER THAN DISCRETE CATEGORIES, BIN BY motionenergy
        clear trls;
        trls = find(data.subjnr==sj & data.correct==corr);
        
        switch field
            case 'rt'
                summaryFunc = @nanmedian;
                distFunc    = @nanstd;
            otherwise
                summaryFunc = @nanmean;
                distFunc    = @naniqr;
        end
        
        % get the mean pupil dilation out
        [grandavg.xMean(find(sj==subjects), find(corr==cors), :), ...
            grandavg.(field).data(find(sj==subjects), find(corr==cors), :), ...
            grandavg.xStd(find(sj==subjects), find(corr==cors), :), ...
            grandavg.(field).wgt(find(sj==subjects), find(corr==cors), :)] = ...
            divideintobins(data.motionstrength(trls), data.(field)(trls), nbins, [], summaryFunc);
    end
    
    % project RT out of the pupil and vice versa
    switch field
        case 'rt'
            data.(field) = projectout(nanzscore(data.rtNorm), nanzscore(data.decision_pupil));
        case 'decision_pupil'
            data.(field) = projectout(data.decision_pupil, nanzscore(data.rtNorm));
    end
    
    % loop over error and correct
    cors = [0 1];
    for corr = cors,
        
        % FIT BETAS ON THE FULL MODEL, NOT BINNED
        trls = find(data.subjnr == sj & data.correct == corr);
        
        % include RT as a regressor
        try
            mdl = fitlm(nanzscore(data.motionstrength(trls)),  ...
                nanzscore(data.(field)(trls)));
            
            % SAVE BETAS FOR THIS PARTICIPANT
            grandavg.(field).regline(find(sj==subjects), find(corr==cors), :) = ...
                mdl.Coefficients.Estimate;
            bint = mdl.coefCI; coef = mdl.Coefficients.Estimate;
            grandavg.(field).bint(find(sj==subjects), find(corr==cors), :) = bint(2,2) - coef(2);
        catch
            grandavg.(field).regline(find(sj==subjects), find(corr==cors), :) = ...
                nan(size(mdl.Coefficients.Estimate));
        end
    end
end

% PLOT
% use nice shades of red and green
colors = cbrewer('qual', 'Set1', 8);
cols = colors([1 2], :);

hold on;
markers = {'^', '.'}; markersizes = [4 14];
for co = 1:2,
    % use double error bars
    h = ploterr(squeeze(nanmean(grandavg.xMean(:, co, :))), ...
        squeeze(nanmean(grandavg.(field).data(:, co, :))), ...
        squeeze(nanstd(grandavg.xMean(:, co, :))) / sqrt(length(subjects)), ...
        squeeze(nanstd(grandavg.(field).data(:, co, :))) / sqrt(length(subjects)), ...
        'k-',  'abshhxy', 0);
    
    set(h(1), 'color', cols(co, :), ...
        'markersize', markersizes(co), ...
        'marker', markers{co});
    if co == 1,
        set(h(1), 'markerfacecolor', 'w', 'markeredgecolor', cols(co, :));
    end
    set(h(2), 'color', cols(co, :));
    set(h(3), 'color', cols(co, :));
    handles{co} = h(1);
end

xlabel('Evidence');
axis square;
xlim([-0.2 5.6]); set(gca, 'xtick', 0:2.75:5.5, ...
    'xticklabel',  {'weak', 'medium', 'strong'}, 'xminortick', 'off');

switch field
    case 'decision_pupil'
        ylim([0.2 0.61]); set(gca, 'ytick', [0.2 0.4 0.6]);
        ylabel('Pupil response (z)');
        savefast(sprintf('%s/Data/GrandAverage/grandavg_pupil_uncertainty.mat', mypath), 'grandavg');
    case 'baseline_pupil'
        ylim([-0.1 0.2]); set(gca, 'ytick', -1:0.1:1);
        ylabel('Next trial baseline pupil (z)');
    case 'rt'
        ylim([0.23 0.56]); set(gca, 'ytick', 0.25:0.15:0.55);
        ylabel('Reaction time (s)');
end

legend([handles{:}], {'error', 'correct'}, 'location', 'southwest'); legend boxoff;
set(gca, 'xcolor', 'k', 'ycolor', 'k');

% output
b = grandavg.(field).regline(:, :, 2);
bint = grandavg.(field).bint;
end
