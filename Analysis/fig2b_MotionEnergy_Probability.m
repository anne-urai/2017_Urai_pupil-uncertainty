function [] = fig2b_MotionEnergy_Probability()
% after filtering motionenergy for all subjects, show the distribution of
% values we get (as a function of nominal coherence level)

if 0,
    % for each sj, check if this worked
    for sj = 1:27,
        t = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
        subplot(5,6,sj);
        plot(t.stim .* t.coherence, t.motionstrength, '.');
        title(sprintf('P%02d', sj)); axis tight;
        box off;
        ylim([-7 7]);
    end
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
colormap hot;
imagesc(stimlevels, edges, n');
set(gca, 'ydir', 'normal');
xlabel({'Stimulus strength (%)'});
ylabel({'Motion energy (a.u.)'});

%offsetAxes
set(gca, 'tickdir', 'out', 'box', 'off');
set(gca, 'xcolor', 'k', 'ycolor', 'k');
xlabs = unique(stim) * 100;
xlabs = {'-30', ' ', ' ', ' ',  ' ', ' ', ' ', '0',  ' ', ' ', ' ',  ' ', ' ', '30'};
set(gca, 'xtick', unique(stim), 'xticklabel', xlabs, 'xticklabelrotation', 0);
set(gca, 'ytick', -6:3:6);
set(gca, 'XAxisLocation','top');

% put a string on top of the colorbar
% make the colorbar prettier, move to the side
c = colorbar('Location', 'SouthOutside');
cpos = c.Position;
cpos(2) = cpos(2)*0.85; % down
cpos(4) = 0.6*cpos(4); % thinner
cpos(1) = cpos(1)*1.06; % to the rigth
cpos(3) = cpos(3) *0.6; % shorter
c.Position = cpos;
c.Box = 'off';

c.Label.String = 'Probability';
c.TickDirection = 'out';
axis square;

end

