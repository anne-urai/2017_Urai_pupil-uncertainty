function [] = figure1
% generate figures to show the uncertainty model from Kepecs et al. 2008
%
% Anne Urai, 2015

global mypath;

% ==================================================================
% 1. decision variable
% ==================================================================

% grey shades
gr1 = [0.2 0.2 0.2];
gr2 = [0.6 0.6 0.6];

% sigma from datas
data = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
b = glmfit(data.motionstrength, (data.resp > 0), 'binomial', 'link', 'probit');

if 0,
    hold on;
    [binsM, binsR] = divideintobins(data.motionstrength, (data.resp > 0), 20);
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

% standard deviation at these values is the inverse!
sigma = 1/b(2);

subplot(331);
% our stimuli range from -6 to 6
stim = -15:0.01:15;

% simulate two motionstrength items
strength1 = 4;
strength2 = 1;

s1 = normpdf(stim, -strength1, sigma);
s2 = normpdf(stim,  strength1, sigma);
s3 = normpdf(stim, -strength2, sigma);
s4 = normpdf(stim,  strength2, sigma);
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

axis tight; ylim([0 0.2]); 
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

axis tight; ylim([0 0.2]);
xlabel('Decision variable');
set(gca, 'xtick', 0, 'xticklabel', {'c'}, ...
    'box', 'off', 'tickdir', 'out', 'ytick', []);
title('Strong evidence');

% ==================================================================
% confidence split by error vs correct
% ==================================================================

% make a little function that computes this quantity
stimlevs = 0:0.01:6;
confE = nan(size(stimlevs)); confC = nan(size(stimlevs));
for s = 1:length(stimlevs),
    [confE(s), confC(s)] = simulateConf(stimlevs(s), sigma);
end

subplot(4,4,4); hold on;
plot(stimlevs, confC, '-',  'color', cols(2, :));
plot(stimlevs, confE, '-', 'color', cols(1, :));

% also add some points
stimpts = [strength2 strength1];
confEpt = nan(size(stimpts)); confCpt = nan(size(stimpts));
for s = 1:length(stimpts),
    [confEpt(s), confCpt(s)] = simulateConf(stimpts(s), sigma);
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
% uncertainty = 1 - confidence
% ==================================================================

% make a little function that computes this quantity
confE = nan(size(stimlevs)); confC = nan(size(stimlevs));
for s = 1:length(stimlevs),
    [confE(s), confC(s)] = simulateUnc(stimlevs(s), sigma);
end

subplot(4,4,8); hold on;
plot(stimlevs, confC, '-',  'color', cols(2, :));
plot(stimlevs, confE, '-', 'color', cols(1, :));

% also add some points
stimpts = [strength2 strength1];
confEpt = nan(size(stimpts)); confCpt = nan(size(stimpts));
for s = 1:length(stimpts),
    [confEpt(s), confCpt(s)] = simulateUnc(stimpts(s), sigma);
end
plot(stimpts(1), confEpt(1), '^', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', gr2, 'MarkerSize', 4);
plot(stimpts(1), confCpt(1), '.', 'MarkerFaceColor', gr2, 'MarkerEdgeColor', gr2, 'MarkerSize', 14);
plot(stimpts(2), confEpt(2), '^', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', gr1, 'MarkerSize', 4);
plot(stimpts(2), confCpt(2), '.', 'MarkerFaceColor', gr1, 'MarkerEdgeColor', gr1, 'MarkerSize', 14);

hold on; axis square;
set(gca,'box', 'off', 'tickdir', 'out', 'xlim', [min(stimlevs) max(stimlevs)]);
xlabel('Evidence');
ylabel('Uncertainty');
ylim([-0.05 0.5]);
xlim([-0.5 max(stimlevs)]);
set(gca, 'ytick', [0:0.5:1]);
set(gca, 'xtick', stimpts, 'xticklabel', {'weak', 'strong'});

cd(mypath); if ~exist('Figures', 'dir'); mkdir Figures; end
print(gcf, '-dpdf', sprintf('%s/Figures/figure1.pdf', mypath));

end


function [confE, confC] = simulateConf(x, sigma)
% compute confidence for correct and error trials based on stimulus strength

% assume no bias
bound = 0;

% simulate decision variables for this level of stimulus strength
dvs = x - bound + sigma * randn(100000, 1);

% find the mean distance to bound for the correct samples
confC = mean(tanh(abs(dvs(dvs > bound) - bound)));

% find the mean distance to bound for error samples
confE = mean(tanh(abs(dvs(dvs < bound) - bound)));

end


function [uncE, uncC] = simulateUnc(x, sigma)
% compute confidence for correct and error trials based on stimulus strength

% assume no bias
bound = 0;

% simulate decision variables for this level of stimulus strength
dvs = x - bound + sigma * randn(100000, 1);

% find the mean distance to bound for the correct samples
uncC = mean(1 - tanh(abs(dvs(dvs > bound) - bound)));

% find the mean distance to bound for error samples
uncE = mean(1 - tanh(abs(dvs(dvs < bound) - bound)));

end
