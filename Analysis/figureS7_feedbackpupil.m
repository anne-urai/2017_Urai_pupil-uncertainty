% figure 4 overview
global mypath;

lagGroups = 1;
mods = {'fbpupil', 'fb+decpupil'};
nbins = 3;
close; figure;

for m = 1:length(mods),
    
    subplot(4,4,(m-1)*4+1); psychFuncShift_Bias(mods{m}, nbins);
    ylim([0.48 0.54]);
    
    %subplot(4,4,(m-1)*4+2); psychFuncShift_Bias_Slope(mods{m}, nbins, []);
    subplot(4,8,(m-1)*8+3);
    
    if m < 3,
        load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, mods{m}));
        plotBetas([dat.response_pupil(:, 1) ...
            dat.stimulus_pupil(:, 1)], ...
            0.5*ones(3,3));
        set(gca, 'xtick', 1:2, 'xticklabel', []);
    else
        load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, mods{m-1}));
        plotBetas([dat.response_rt(:, 1) ...
            dat.stimulus_rt(:, 1)], ...
            0.5*ones(3,3));
        set(gca, 'xtick', 1:2, 'xticklabel', ...
            {'Pupil x choice', 'Pupil x stimulus'}, ...
            'xticklabelrotation', -30);
    end
    ylim([-0.09 0.06]);
end

set(gca, 'xtick', 1:2, 'xticklabel', ...
    {'Pupil x choice', 'Pupil x stimulus'}, ...
    'xticklabelrotation', -30);

print(gcf, '-dpdf', sprintf('%s/Figures/figureS4.pdf', mypath));
