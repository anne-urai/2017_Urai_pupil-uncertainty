function [] = figure1b
% generate figures to show the uncertainty model from Kepecs et al. 2008
% in this version, rather than fitting on my data, use Kepecs parameters
%
% Anne Urai, 2015

global mypath;

% grey shades for different levels of difficulty
gr1 = [0.2 0.2 0.2];
gr2 = [0.6 0.6 0.6];

% ==================================================================
% 1. decision variable
% ==================================================================

% sigma from data - as in Kepecs
data    = table;
a       = linspace(0.1, 99.1, 100000);
b       = fliplr(a);
data.motionstrength = log(a ./ b);

% standard deviation at these values is the inverse!
sigma = 0.5;

% our stimuli range from -6 to 6
stim = -6:0.01:6;

% simulate two motionstrength items
strength1 = 1;
strength2 = 0.3;

s1 = normpdf(stim, -strength1, sigma);
s2 = normpdf(stim,  strength1, sigma);
s3 = normpdf(stim, -strength2, sigma);
s4 = normpdf(stim,  strength2, sigma);

subplot(331);
p = plot(stim, s1, stim, s2, stim, s3, stim, s4, zeros(1, 2), [0 max(s2)], 'k');
p(1).Color = gr1; p(2).Color = gr1;
p(3).Color = gr2; p(4).Color = gr2;

hold on;
pt = 1.8; % illustrate a single DV
plot(pt, s2(dsearchn(stim', pt)), '.', ...
    'MarkerFaceColor', gr1, 'MarkerEdgeColor', gr1, 'MarkerSize', 10);
plot([pt pt], [0 s2(dsearchn(stim', pt))], 'color', gr1, 'linestyle', ':');

axis tight; xlim([-3 3]);
ylabel('Probability density'); xlabel('Decision variable');
text(-2, 0.8, 'Stimulus A');
text(2, 0.8, 'Stimulus B');

% add legend
l = legend([p(3), p(1)], {'weak', 'strong'});
lpos = l.Position;
lpos(1) = lpos(1) + 0.12; % right
lpos(2) = lpos(2) - 0.05; % down
l.Position = lpos;
l.Box = 'off';
text(3, 0.16, 'Evidence');

%offsetAxes(gca, 0.1, 1);
set(gca, 'xtick', [0 pt], 'xticklabel', {'c', 'dv_i'}, ...
    'box', 'off', 'tickdir', 'out', 'ytick', [], 'yticklabel', []);

% ==================================================================
% average confidence
% ==================================================================

subplot(443);
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

axis tight; xlim([-3 3]);
ylabel('Probability density');
title('Weak evidence');
set(gca, 'ytick', [], 'xtick', [], 'box', 'off');

% then the easy stimulus
subplot(4,4,7);
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

axis tight; xlim([-3 3]);
xlabel('Decision variable');
set(gca, 'xtick', 0, 'xticklabel', {'c'}, ...
    'box', 'off', 'tickdir', 'out', 'ytick', []);
title('Strong evidence');

% ==================================================================
% generate a measure of uncertainty for each sample
% ==================================================================

model.evs           = unifrnd(min(data.motionstrength), max(abs(data.motionstrength)), 1, 1e7);
model.dvs           = model.evs + normrnd(0, sigma, size(model.evs));
model.choice        = sign(model.dvs);
model.correct       = (model.choice == sign(model.evs));
model.confidence    = tanh(abs(model.dvs));
model.uncertainty   = 1 - model.confidence;

% ==================================================================
% plot mean confidence
% ==================================================================

subplot(4,4,4); hold on;
[stimlevs, confC] = divideintobins(abs(model.evs(model.correct == 1)), model.confidence(model.correct == 1), 4);
plot(stimlevs, confC, '-',  'color', cols(2, :));
[stimlevs, confE] = divideintobins(abs(model.evs(model.correct == 0)), model.confidence(model.correct == 0), 4);
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
lpos = get(l, 'position'); lpos(2) = lpos(2) - 0.15;
set(l, 'box', 'off', 'position', lpos);

% ==================================================================
% plot mean uncertainty below that
% ==================================================================

subplot(4,4,8); hold on;
[stimlevs, confC] = divideintobins(abs(model.evs(model.correct == 1)), model.uncertainty(model.correct == 1), 20);
plot(stimlevs, confC, '-',  'color', cols(2, :));
[stimlevs, confE] = divideintobins(abs(model.evs(model.correct == 0)), model.uncertainty(model.correct == 0), 20);
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

subplot(4,4,9); hold on;
[unc, correct] = divideintobins(model.uncertainty, model.correct, 100);
plot(unc, correct*100, 'k-');
ylim([45 100]); set(gca, 'ytick', [50 75 100]); ylabel('Accuracy (%)');
xlim([-0.05 1.02]); xlabel('Uncertainty');
plot([0 1], [50 50], 'color', [0.5 0.5 0.5], 'linewidth', 0.5, 'linestyle', ':'); hold on;
box off; 

% ==================================================================
% from Sanders, psychFuncs by uncertainty
% ==================================================================

uncMed = quantile(model.uncertainty, 2); % median
% PLOT
subplot(4,4,10);
hold on;
colors(1,:) = [0.5 0.5 0.5];
colors(2,:) = [0.2 0.2 0.2];
newx = linspace(min(abs(model.evs)), max(abs(model.evs)), 100);

for r = 1:2,
    switch r
        case 1
            trls = find(model.uncertainty < uncMed(1));
        case 2
            trls = find(model.uncertainty > uncMed(2));
    end
    % [ev, acc] = divideintobins(abs(model.evs(trls)), model.correct(trls), nbins);
    % fit cumulative weibull to this psychfunc
    
    % make weibull fit faster
    [binnedx, binnedy] = divideintobins(abs(model.evs(trls)), model.correct(trls), 100);
    [slope, threshold, lapse] = fitWeibull(binnedx, binnedy);
    h = plot(newx, 100* Weibull([slope, threshold, lapse], newx), 'color', colors(r, :)) ;
    
    % datapoints on top
   % h = plot(ev, acc*100, '.', 'markersize', 12', 'color', colors(r, :));
   handles{r} = h;
end

set(gca, 'xlim', [-0.5 5.5], 'ylim', [45 100], 'ytick', 50:25:100);
xlim([-0.5 5.5]); set(gca, 'xtick', [0 2.75 5.5], 'xticklabel', {'weak', 'medium', 'strong'}, 'xminortick', 'off');
ylabel('Accuracy (%)'); xlabel('Evidence');
 box off;

l = legend([handles{:}], {'low', 'high'}, 'location', 'southeast');
legend boxoff;
lpos = get(l, 'position');
lpos(1) = lpos(1) + 0.05;
set(l, 'position', lpos);

% save

cd(mypath); if ~exist('Figures', 'dir'); mkdir Figures; end
print(gcf, '-dpdf', sprintf('%s/Figures/figure1.pdf', mypath));



