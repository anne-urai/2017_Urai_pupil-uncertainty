function [] = fig2b_MotionEnergy_Probability()
% after filtering motionenergy for all subjects, show the distribution of
% values we get (as a function of nominal coherence level)
close all;
figpath = '~/Dropbox/Figures/uncertainty';

if 0,
    % for each sj, check if this worked
    for sj = 1:27,
        t = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
        subplot(5,6,sj);
        plot(t.stim .* t.coherence, t.motionstrength, '.');
        ylim([-1 1]); xlim([-0.35 0.35]);
        title(sprintf('P%02d', sj)); axis tight;
    end
    close;
end

%% group level results
t = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_allsj.csv'));

stim = t.coherence;
stim = round((stim) * 100000)/100000;
stim((abs(stim-0.01)) < 0.000001) = 0.0125;
stim = t.stim .* stim;
strength = t.motionstrength;

% !!! add zero stimlevels into image MAT to ensure linearity
stepsize = min(diff(unique(stim)));
%stepsize = 2*min(diff(unique(stim)));
stimlevels = min(unique(stim)):stepsize:max(unique(stim));
stimlevels = [-0.31 stimlevels 0.31];
n = zeros(length(stimlevels), 99);

% find
for s = 1:length(stimlevels),
    trls = find(abs(stim - stimlevels(s)) < 0.000001);
    edges = linspace(min(strength), max(strength), 100);
    [n(s, :), edges] = histcounts(strength(trls), edges, ...
        'Normalization', 'probability');
end

% check
assert(isequal(length(unique(stim)), length(find(nansum(n, 2) > 0))), 'spacing not quite right');

subplot(441);
colormap hot;
imagesc(stimlevels, edges, n');
set(gca, 'ydir', 'normal');
xlabel({'Stimulus difficulty (%)'});
ylabel({'Motion energy (a.u.)'});

%offsetAxes
set(gca, 'tickdir', 'out', 'box', 'off');
set(gca, 'xcolor', 'k', 'ycolor', 'k');
xlabs = unique(stim) * 100;
xlabs = {'-30', ' ', ' ', ' ',  ' ', ' ', ' ', '0',  ' ', ' ', ' ',  ' ', ' ', '30'};
set(gca, 'xtick', unique(stim), 'xticklabel', xlabs, 'xticklabelrotation', 0);
set(gca, 'ytick', -3:3:3);
set(gca, 'XAxisLocation','top')

% make the colorbar prettier, move to the side
c = colorbar('Location', 'SouthOutside');
cpos = c.Position;
cpos(2) = cpos(2)* 0.9; % up
cpos(4) = 0.6*cpos(4); % thinner
cpos(1) = cpos(1)*1.15; % to the rigth
cpos(3) = cpos(3) *0.75; % shorter
c.Position = cpos;
c.Box = 'off';

axis square
% put a string on top of the colorbar
c.Label.String = 'Probability';
c.TickDirection = 'out';
print(gcf, '-dpdf', sprintf('%s/Fig1b_motionenergyResult.pdf', figpath));

end

