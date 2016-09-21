function [] = figure1
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
close all;

% grey shades for different levels of difficulty
gr1 = [0.2 0.2 0.2];
gr2 = [0.6 0.6 0.6];

% ==================================================================
% 1. decision variable
% ==================================================================

% sigma from data
data = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
b = glmfit(data.motionstrength, (data.resp > 0), 'binomial', 'link', 'probit');

% standard deviation at these values is the inverse
sigma = 1/b(2);
stim = -15:0.01:15; % x axis for distributions

% simulate two motionstrength items
strength1 = 4;
strength2 = 1;

if 0,
    hold on;
    [binsM, binsR] = divideintobins(data.motionstrength, (data.resp > 0), 3);
    plot(binsM, binsR, '.'); hold on;
    y = glmval(b, sort(data.motionstrength),  'probit');
    % this is the code used by glmval
    % eta = x*beta + offset;
    % lowerBnd = norminv(eps(dataClass)); upperBnd = -lowerBnd;
    % ilink = @(eta) normcdf(constrain(eta,lowerBnd,upperBnd));
    
    plot(sort(data.motionstrength), y, '.');
    
    % how to transfer probit measures into normcdf params?
    plot(sort(data.motionstrength), normcdf(sort(data.motionstrength), b(1), 1/b(2)), '.');
    plot(sort(data.motionstrength), normcdf(sort(data.motionstrength), -b(1), 1/b(2)), '.');
    
    grid on; xlim([-0.05 0.05]);
    legend({'data', 'glmval', 'normpdf1', 'normpdf2'});
    
end

s1 = normpdf(stim, -strength1, sigma);
s2 = normpdf(stim,  strength1, sigma);
s3 = normpdf(stim, -strength2, sigma);
s4 = normpdf(stim,  strength2, sigma);

subplot(449);
p = plot(stim, s1, stim, s2, stim, s3, stim, s4, zeros(1, 2), [0 max(s2)], 'k');
p(1).Color = gr1; p(2).Color = gr1;
p(3).Color = gr2; p(4).Color = gr2;

hold on;
pt = 6;
plot(pt, s2(dsearchn(stim', pt)), '.', ...
    'MarkerFaceColor', gr1, 'MarkerEdgeColor', gr1, 'MarkerSize', 10);
plot([pt pt], [0 s2(dsearchn(stim', pt))], 'color', gr1, 'linestyle', ':');

axis tight; ylim([0 0.2]);
ylabel('Probability density'); xlabel('Decision variable');
text(-12, 0.18, 'Stimulus A');
text(2, 0.18, 'Stimulus B');

% add legend
l = legend([p(3), p(1)], {'weak', 'strong'});
lpos = l.Position;
lpos(1) = lpos(1) + 0.12; % right
lpos(2) = lpos(2) - 0.05; % down
l.Position = lpos;
l.Box = 'off';
text(15, 0.16, 'Evidence');

%offsetAxes(gca, 0.1, 1);
set(gca, 'xtick', [0 pt], 'xticklabel', {'c', 'dv_i'}, ...
    'box', 'off', 'tickdir', 'out', 'ytick', [], 'yticklabel', []);

% ==================================================================
% average confidence
% ==================================================================

subplot(551);
colors = cbrewer('qual', 'Set1', 8);
cols = colors([1 2], :);

% first the difficult stimulus
% indicate one correct and one error
hold on;
a = area(unique(abs(stim)), normpdf(unique(abs(stim)), strength2, sigma));
a.FaceColor = cols(2,:); a.EdgeColor = 'none';

a = area(unique(-abs(stim)), normpdf(unique(-abs(stim)), strength2, sigma));
a.FaceColor = cols(1,:); a.EdgeColor = 'none';

% now the distributions on top
s3 = normpdf(stim, -strength2, sigma);
s4 = normpdf(stim,  strength2, sigma);
p = plot(stim, s3, ':', stim, s4,  '-', zeros(1, 2), [0 max(s4)], 'k');
p(1).Color = gr2; p(2).Color = gr2;

axis tight; ylim([0 0.2]);
ylabel('Probability density');
title('Weak evidence');
set(gca, 'ytick', [], 'xtick', [], 'box', 'off');

% then the easy stimulus
subplot(556);
% indicate one correct and one error
hold on;
a = area(unique(abs(stim)), normpdf(unique(abs(stim)), strength1, sigma));
a.FaceColor = cols(2,:); a.EdgeColor = 'none';

a = area(unique(-abs(stim)), normpdf(unique(-abs(stim)), strength1, sigma));
a.FaceColor = cols(1,:); a.EdgeColor = 'none';

% now the distributions on top
s3 = normpdf(stim, -strength1, sigma);
s4 = normpdf(stim,  strength1, sigma);
p = plot(stim, s3, ':', stim, s4,  '-', zeros(1, 2), [0 max(s4)], 'k');
p(1).Color = gr1; p(2).Color = gr1;

axis tight; ylim([0 0.2]);
xlabel('Decision variable');
set(gca, 'xtick', 0, 'xticklabel', {'c'}, ...
    'box', 'off', 'tickdir', 'out', 'ytick', []);
title('Strong evidence');

% ==================================================================
% generate a measure of uncertainty for each sample
% ==================================================================

% sigma = 1;
model.evs           = unifrnd(min(data.motionstrength), max(data.motionstrength), 1, 1e7);
model.dvs           = model.evs + normrnd(0, sigma, size(model.evs));
model.choice        = sign(model.dvs);
model.correct       = (model.choice == sign(model.evs));

% normal cdf on absolute dv, see Lak et al. equation 7
dv2conf             = @(x, sigma) 0.5 * (1 + erf(abs(x) ./ (sqrt(2)*sigma)));
model.confidence    = dv2conf(model.dvs, sigma);
model.uncertainty   = 1 - model.confidence;

% ==================================================================
% plot mean confidence
% ==================================================================

nbins = 100;
subplot(552); hold on;
[stimlevs, confC] = divideintobins(abs(model.evs(model.correct == 1)), model.confidence(model.correct == 1), nbins);
plot(stimlevs, confC, '-',  'color', cols(2, :));
[stimlevs, confE] = divideintobins(abs(model.evs(model.correct == 0)), model.confidence(model.correct == 0), nbins);
plot(stimlevs, confE, '-', 'color', cols(1, :));

% also add some points
stimpts = [strength2 strength1];
confEpt = nan(size(stimpts)); confCpt = nan(size(stimpts));
for s = 1:length(stimpts),
    confEpt(s) = mean(model.confidence(model.evs < stimpts(s)*1.1 & model.evs > stimpts(s) * 0.9 & model.correct == 0));
    confCpt(s) = mean(model.confidence(model.evs < stimpts(s)*1.1 & model.evs > stimpts(s) * 0.9 & model.correct == 1));
end
plot(stimpts(1), confEpt(1), '^', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', gr2, 'MarkerSize', 4);
plot(stimpts(1), confCpt(1), '.', 'MarkerFaceColor', gr2, 'MarkerEdgeColor', gr2, 'MarkerSize', 14);
ploth(1) = plot(stimpts(2), confEpt(2), '^', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', gr1, 'MarkerSize', 4);
ploth(2) = plot(stimpts(2), confCpt(2), '.', 'MarkerFaceColor', gr1, 'MarkerEdgeColor', gr1, 'MarkerSize', 14);

axis tight; axis square;
set(gca,  'box', 'off', 'tickdir', 'out', 'xlim', [min(stimlevs) max(stimlevs)]);
ylabel('Confidence');
ylim([0.5 1]);
xlim([-0.5 max(stimlevs)]);
set(gca, 'ytick', [0:0.5:1]);
set(gca, 'xcolor', 'w');

l = legend(ploth, {'error', 'correct'});
lpos = get(l, 'position'); lpos(2) = lpos(2) - 0.10;
set(l, 'box', 'off', 'position', lpos);

% ==================================================================
% plot mean uncertainty below that
% ==================================================================

subplot(557); hold on;
[stimlevs, confC] = divideintobins(abs(model.evs(model.correct == 1)), model.uncertainty(model.correct == 1), nbins);
plot(stimlevs, confC, '-',  'color', cols(2, :));
[stimlevs, confE] = divideintobins(abs(model.evs(model.correct == 0)), model.uncertainty(model.correct == 0), nbins);
plot(stimlevs, confE, '-', 'color', cols(1, :));

% also add some points
stimpts = [strength2 strength1];
confEpt = nan(size(stimpts)); confCpt = nan(size(stimpts));
for s = 1:length(stimpts),
    confEpt(s) = mean(model.uncertainty(model.evs < stimpts(s)*1.1 & model.evs > stimpts(s) * 0.9 & model.correct == 0));
    confCpt(s) = mean(model.uncertainty(model.evs < stimpts(s)*1.1 & model.evs > stimpts(s) * 0.9 & model.correct == 1));
end
plot(stimpts(1), confEpt(1), '^', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', gr2, 'MarkerSize', 4);
plot(stimpts(1), confCpt(1), '.', 'MarkerFaceColor', gr2, 'MarkerEdgeColor', gr2, 'MarkerSize', 14);
ploth(1) = plot(stimpts(2), confEpt(2), '^', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', gr1, 'MarkerSize', 4);
ploth(2) = plot(stimpts(2), confCpt(2), '.', 'MarkerFaceColor', gr1, 'MarkerEdgeColor', gr1, 'MarkerSize', 14);

axis square;
set(gca,'box', 'off', 'tickdir', 'out', 'xlim', [min(stimlevs) max(stimlevs)]);
xlabel('Evidence');
ylabel('Uncertainty');
ylim([-0.05 0.5]);
xlim([-0.5 max(stimlevs)]);
set(gca, 'ytick', [0:0.5:1]);
set(gca, 'xtick', stimpts, 'xticklabel', {'weak', 'strong'});

% ==================================================================
% from Sanders, show accuracy vs confidence
% ==================================================================

[unc, correct] = divideintobins(model.uncertainty, model.correct, nbins);

subplot(553); cla; hold on;
plot(unc, correct*100, 'k-'); % show
ylim([45 100]); set(gca, 'ytick', [50 75 100]); ylabel('Accuracy (%)');
xlim([-0.02 max(unc)]); xlabel('Uncertainty');
plot([min(unc) max(unc)], [50 50], 'color', [0.5 0.5 0.5], 'linewidth', 0.5, 'linestyle', ':');
set(gca, 'xtick', [min(unc) max(unc)], 'xticklabel', {'low', 'high'});
box off;  axis square;

% ==================================================================
% from Sanders, psychFuncs by uncertainty
% ==================================================================

uncMed = median(model.uncertainty); % median split
% PLOT
subplot(558);
hold on;
colors(1,:) = [0.5 0.5 0.5];
colors(2,:) = [0.2 0.2 0.2];
for r = 1:2,
    switch r
        case 1
            trls = find(model.uncertainty < uncMed);
        case 2
            trls = find(model.uncertainty > uncMed);
    end

    % make weibull fit faster
    [binnedx, binnedy] = divideintobins(abs(model.evs(trls)), model.correct(trls), 100);
    h = plot(binnedx, 100*binnedy, 'color', colors(r, :)) ;
    handles{r} = h;
end

set(gca, 'xlim', [-0.5 5.5], 'ylim', [45 100], 'ytick', 50:25:100);
xlim([-0.5 5.5]); set(gca, 'xtick', [0 2.75 5.5], 'xticklabel', {'weak', 'medium', 'strong'}, 'xminortick', 'off');
ylabel('Accuracy (%)'); xlabel('Evidence');
box off;  axis square;
l = legend([handles{:}], {'low', 'high'}, 'location', 'southeast');
legend boxoff;
lpos = get(l, 'position');
lpos(1) = lpos(1) + 0.05;
set(l, 'position', lpos);

% save
print(gcf, '-dpdf', sprintf('%s/Figures/Figure1.pdf', mypath));



