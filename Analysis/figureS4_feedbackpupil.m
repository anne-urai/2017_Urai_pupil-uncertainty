% figure 4 overview
global mypath;

lagGroups = 1;
mods = {'fbpupil', 'fb+decpupil', 'dec+fbpupil'};
nbins = 3;
close;
figure;

for m = 1:length(mods),
    
    subplot(4,4,(m-1)*4+1); psychFuncShift_Bias_byResp(mods{m}, nbins);
    subplot(4,4,(m-1)*4+2); psychFuncShift_Bias_Slope(mods{m}, nbins, []);
    
    subplot(4,8,(m-1)*8+6);
    
    if m  < 3,
        load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, mods{m}));
        plotBetasSwarm([dat.response_pupil(:, 1) ...
            dat.stimulus_pupil(:, 1)], ...
            [0 0 0; 0 0 0; 0 0 0; 0 0 0]);
        set(gca, 'xtick', 1:2, 'xticklabel', []);
    else
        load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, mods{m-1}));
        plotBetasSwarm([dat.response_rt(:, 1) ...
            dat.stimulus_rt(:, 1)], ...
            [0 0 0; 0 0 0; 0 0 0; 0 0 0]);
        set(gca, 'xtick', 1:2, 'xticklabel', ...
            {'Pupil * choice', 'Pupil * stimulus'}, ...
            'xticklabelrotation', -30);
    end
end

print(gcf, '-dpdf', sprintf('%s/Figures/figureS4.pdf', mypath));
