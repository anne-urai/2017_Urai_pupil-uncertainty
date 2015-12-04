% ============================================ %
% correlation between history effect and pupil interaction
% ============================================ %
% group lags 4-7

whichmodulator = 'pupilR';
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));
colors = linspecer(8);

names = fieldnames(dat);
names = names(1:end-2); % remove the two pvals
for i = 1:length(names),
    olddat = dat.(names{i});
    newdat(:, 1:3) = olddat(:, 1:3);
    newdat(:, 4) = mean(olddat(:, 4:end), 2);
    datG.(names{i}) = newdat;
end

flds = {'response', 'stimulus', 'correct', 'incorrect'};
colorder = [2 4 3 1];
for l = 1:4,
    for f = 1:length(flds),
        subplot(4,4,f+(l-1)*4);
        
        plot(datG.(flds{f})(:, l), datG.([flds{f} '_pupil'])(:, l), ...
            'o', 'MarkerFaceColor', colors(colorder(f), :), 'MarkerEdgeColor', 'w');
        box off;
        
        % do stats
        [rho, pval] = corr(datG.(flds{f})(:, l), datG.([flds{f} '_pupil'])(:, l), 'type', 'spearman');
        if pval < 0.05,
            ls = lsline; set(ls, 'color', [0.7 0.7 0.7]);
          %  assert(1==0)
            disp([rho pval]);
            text(-0.3, -0.15, sprintf('rho = %.3f, p = %.3f', rho, pval'));
        end
        
        % layout
        if l == 1,
            title(flds{f});
        end
        if f == 1,
            if l == 4,
                ylabel('Lag 4-7');
            else
                ylabel(sprintf('Lag %d', l));
            end
        end
    end
end

suplabel('Pupil interaction weight', 'y');
suplabel('History weight', 'x');

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/serial/historyPupil_kernels_scatter.pdf'));

%%
for l = 1:4,
    for f = 1:length(flds),

        mainW = datG.(flds{f})(:, l);
        intW = datG.([flds{f} '_pupil'])(:, l);
        [h,p,ci,stats] = ttest2(intW(mainW<0),intW(mainW>0));
        if p < 0.05,
         %   assert(1==0);
        end
       
    end
end