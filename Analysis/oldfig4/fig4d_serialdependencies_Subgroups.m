function  fig4d_serialdependencies_Subgroups

% split subjects based on their plain history weights
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'plain'));
colors = linspecer(8);
clf;

ylim1st = [-.3 .3];
ylim2nd = [-0.1 0.1];

for gr = 1:2,
    
    % ============================================ %
    % select subgroup
    % ============================================ %
    load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'plain'));
    switch gr
        case 1
            theseSj = find(dat.response(:, 1) > 0);
            titcolor = colors(2,:);
            tit = 'Repeaters';
        case 2
            theseSj = find(dat.response(:, 1) < 0);
            titcolor = colors(5,:);
            tit = 'Alternators';
    end
    
    whichmodulator = 'pupil';
    load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));
    
    fldnms = fieldnames(dat);
    for f = 1:length(fldnms),
        dat.(fldnms{f}) = dat.(fldnms{f})(theseSj, :);
    end
    
    % ============================================ %
    % grand average history kernels
    % ============================================ %
    
    % stim resp
    subplot(4,4,(gr-1)*4+1); hold on;
    bar(1, nanmean(dat.stimulus(:, 1)), 'facecolor', colors(2, :), 'edgecolor', 'none');
    bar(2, nanmean(dat.response(:, 1)), 'facecolor', colors(4, :), 'edgecolor', 'none');
    errorbar(1:2, [nanmean(dat.stimulus(:, 1)) nanmean(dat.response(:, 1))], ...
        [nanstd(dat.stimulus(:, 1)) nanstd(dat.response(:, 1))] ./ sqrt(length(theseSj)), '.');
    set(gca, 'xtick', 1:2, 'xticklabel', {'stim', 'resp'});
    xlim([0.5 2+0.5]); ylim(ylim1st);
    ylims = get(gca, 'ylim'); set(gca, 'ylim', [-max(abs(ylims)*1.1) max(abs(ylims)*1.1)]);
    
    ylabel(tit, 'color', titcolor);
    
    % err correct
    subplot(4,4,(gr-1)*4+2); hold on;
    bar(1, nanmean(dat.correct(:, 1)), 'facecolor', colors(3, :), 'edgecolor', 'none');
    bar(2, nanmean(dat.incorrect(:, 1)), 'facecolor', colors(1, :), 'edgecolor', 'none');
    errorbar(1:2, [nanmean(dat.correct(:, 1)) nanmean(dat.incorrect(:, 1))], ...
        [nanstd(dat.correct(:, 1)) nanstd(dat.incorrect(:, 1))] ./ sqrt(length(theseSj)), '.');
    set(gca, 'xtick', 1:2, 'xticklabel', {'correct', 'error'});
    xlim([0.5 2+0.5]); ylim(ylim1st);
    ylims = get(gca, 'ylim'); set(gca, 'ylim', [-max(abs(ylims)*1.1) max(abs(ylims)*1.1)]);
    
    % now the pupil
    subplot(4,4,(gr-1)*4+3);
    hold on;
    bar(1, nanmean(dat.stimulus_pupil(:, 1)), 'facecolor', colors(2, :), 'edgecolor', 'none');
    bar(2, nanmean(dat.response_pupil(:, 1)), 'facecolor', colors(4, :), 'edgecolor', 'none');
    errorbar(1:2, [nanmean(dat.stimulus_pupil(:, 1)) nanmean(dat.response_pupil(:, 1))], ...
        [nanstd(dat.stimulus_pupil(:, 1)) nanstd(dat.response_pupil(:, 1))] ./ sqrt(length(theseSj)), '.');
    set(gca, 'xtick', 1:2, 'xticklabel', {'stim*pup', 'resp*pup'});
    xlim([0.5 2+0.5]); ylim(ylim2nd);
    ylims = get(gca, 'ylim'); set(gca, 'ylim', [-max(abs(ylims)*1.1) max(abs(ylims)*1.1)]);
    
    subplot(4,4,(gr-1)*4+4); hold on;
    bar(1, nanmean(dat.correct_pupil(:, 1)), 'facecolor', colors(3, :), 'edgecolor', 'none');
    bar(2, nanmean(dat.incorrect_pupil(:, 1)), 'facecolor', colors(1, :), 'edgecolor', 'none');
    errorbar(1:2, [nanmean(dat.correct_pupil(:, 1)) nanmean(dat.incorrect_pupil(:, 1))], ...
        [nanstd(dat.correct_pupil(:, 1)) nanstd(dat.incorrect_pupil(:, 1))] ./ sqrt(length(theseSj)), '.');
    set(gca, 'xtick', 1:2, 'xticklabel', {'correct', 'error'});
    xlim([0.5 2+0.5]); ylim(ylim2nd);
    ylims = get(gca, 'ylim'); set(gca, 'ylim', [-max(abs(ylims)*1.1) max(abs(ylims)*1.1)]);
    
end

%suplabel('Grand Average', 't');
print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/serial/historykernels_%s_subgroups.pdf', whichmodulator));

end
