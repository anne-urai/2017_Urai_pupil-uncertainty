function [] = figure1_withLapses
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
% Anne Urai, 2017
% anne.urai@gmail.com

global mypath;
close all;

colors = cbrewer('qual', 'Set1', 8);
cols   = colors([1 2], :);
nbins  = 10;

% sigma from data
data    = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
b       = glmfit(data.motionstrength, (data.resp > 0), 'binomial', 'link', 'probit');

% standard deviation at these values is the inverse
sigma   = 1/b(2);

%% ==================================================================
% generate a measure of uncertainty for each sample
% ==================================================================

model.evs           = unifrnd(min(data.motionstrength), max(data.motionstrength), 1, 1e7);
% dv is evidence corrupted by noise
model.dvs           = model.evs + normrnd(0, sigma, size(model.evs)); 

% add the possibility of a lapse on each trial
lapserate                   = 0.05; % 10%
model.lapse                 = logical(binornd(1, lapserate, size(model.dvs)));

% on no-lapse trials, base the choice on the noisy evidence
model.choice                = sign(model.dvs);
model.choice(model.choice == -1) = 0;

% on a lapse trial, randomly guess by flipping a coin
model.choice(model.lapse)   = binornd(1, 0.5, 1, sum(model.lapse)); 

% recode into [-1,1]
model.choice(model.choice == 0) = -1;

% did the model make the objective correct decision?
model.correct               = (model.choice == sign(model.evs));

% normal cdf on absolute dv, see Lak et al. equation 7
dv2conf             = @(x, sigma) 0.5 * (1 + erf(abs(x) ./ (sqrt(2)*sigma)));
model.confidence    = dv2conf(model.dvs, sigma);

% if the trial was a lapse, confidence must be 0.5 
model.confidence(model.lapse) = 0.5;
model.uncertainty   = 1 - model.confidence;

%% ==================================================================
% SHOW PSYCHOMETRIC FUNCTION WITH LAPSE RATE
% ==================================================================

figure; subplot(441);
[stimlevels, pchoice] = divideintobins(model.evs, (model.choice > 0), nbins);
plot(stimlevels, pchoice, 'k');
xlabel('Sensory evidence');
ylabel('P(choice == 1)');
axis tight; axis square; box off;
ylim([0 1]); xlim([-6 6]); set(gca, 'xtick', [-6 0 6]);

%% ==================================================================
% plot mean confidence
% ==================================================================

subplot(442); hold on;
[stimlevs, confC] = divideintobins(abs(model.evs(model.correct == 1)), model.confidence(model.correct == 1), nbins);
plot(stimlevs, confC, '-',  'color', cols(2, :));
[stimlevs, confE] = divideintobins(abs(model.evs(model.correct == 0)), model.confidence(model.correct == 0), nbins);
plot(stimlevs, confE, '-', 'color', cols(1, :));

axis tight; axis square;
set(gca,  'box', 'off', 'tickdir', 'out', 'xlim', [min(stimlevs) max(stimlevs)]);
ylabel('Confidence');
ylim([0.5 1]);
xlim([-0.5 max(stimlevs)]);
set(gca, 'ytick', [0:0.5:1]);
set(gca, 'xtick', [-0.5 max(stimlevs)], 'xticklabel', {'weak', 'strong'});
%offsetAxes;

%% ==================================================================
% plot mean uncertainty below that
% ==================================================================

subplot(443); hold on;
[stimlevs, confC] = divideintobins(abs(model.evs(model.correct == 1)), model.uncertainty(model.correct == 1), nbins);
plot(stimlevs, confC, '-',  'color', cols(2, :));
[stimlevs, confE] = divideintobins(abs(model.evs(model.correct == 0)), model.uncertainty(model.correct == 0), nbins);
plot(stimlevs, confE, '-', 'color', cols(1, :));

axis square;
set(gca,'box', 'off', 'tickdir', 'out', 'xlim', [min(stimlevs) max(stimlevs)]);
xlabel('Evidence');
ylabel('Uncertainty');
ylim([0 0.5]);
xlim([-0.5 max(stimlevs)]);
set(gca, 'ytick', [0:0.5:1]);
set(gca, 'xtick', [-0.5 max(stimlevs)], 'xticklabel', {'weak', 'strong'});
%offsetAxes;

%% ==================================================================
% how does the lapse rate scale?
% ==================================================================

subplot(444); hold on;
[stimlevs, confC] = divideintobins(abs(model.evs(model.correct == 1)), model.lapse(model.correct == 1), nbins);
plot(stimlevs, confC, '-',  'color', cols(2, :));
[stimlevs, confE] = divideintobins(abs(model.evs(model.correct == 0)), model.lapse(model.correct == 0), nbins);
plot(stimlevs, confE, '-', 'color', cols(1, :));

axis square;
set(gca,'box', 'off', 'tickdir', 'out', 'xlim', [min(stimlevs) max(stimlevs)]);
xlabel('Evidence');
ylabel('P(lapse)');
ylim([0 0.5]); set(gca, 'ytick', [0 0.5]);
xlim([min(stimlevs) max(stimlevs)]);
set(gca, 'ytick', [0:0.5:1]);
set(gca, 'xtick', [min(stimlevs) max(stimlevs)], 'xticklabel', {'weak', 'strong'});
%offsetAxes;

%% ==================================================================
% from Sanders, show accuracy vs confidence
% ==================================================================

[unc, correct] = divideintobins(model.confidence, double(model.correct), nbins);

subplot(446); hold on;
plot(unc, correct*100, 'k-'); % show
ylim([50 100]); set(gca, 'ytick', [50 75 100]); ylabel('Accuracy (%)');
xlim([min(unc) max(unc)]); xlabel('Confidence');
set(gca, 'xtick', [min(unc) max(unc)], 'xticklabel', {'low', 'high'});
box off;  axis square;

[unc, correct] = divideintobins(model.uncertainty, double(model.correct), nbins);

subplot(447); hold on;
plot(unc, correct*100, 'k-'); % show
ylim([50 100]); set(gca, 'ytick', [50 75 100]); ylabel('Accuracy (%)');
xlim([min(unc) max(unc)]); xlabel('Uncertainty');
set(gca, 'xtick', [min(unc) max(unc)], 'xticklabel', {'low', 'high'});
box off;  axis square;

% also for lapse rate
[unc, correct] = divideintobins(double(model.lapse), double(model.correct), 2);

subplot(448); hold on;
bar(unc, correct*100, 'edgecolor', 'none', 'facecolor', [0.5 0.5 0.5]); % show
ylim([45 100]); set(gca, 'ytick', [50 100]); ylabel('Accuracy (%)');
xlim([-0.5 1.5]); xlabel('Lapse');
set(gca, 'xtick', [min(unc) max(unc)], 'xticklabel', {'no', 'yes'});
box off;  axis square;

% save
print(gcf, '-dpdf', sprintf('%s/Figures/Figure1_withLapses.pdf', mypath));

