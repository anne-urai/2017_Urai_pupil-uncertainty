function [] = fig4b_decisionStrategies_pupilSplit(whichmodulator)

if ~exist('whichmodulator', 'var'); whichmodulator = 'rt'; end

cnt = 1;
for m = 1:2,
    switch m
        case 1
            whichMod = 'decision_pupil'
        case 2
            whichMod = 'rt';
    end
    
    nbins = 3;
    % ========================================================= %
    % split in 3 pupil bins and show the response and stimulus history weights
    % ========================================================= %
    
    for sj = 1:27,
        data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
        
        % make design matrix
        designM = [data.motionstrength ...
            circshift(data.resp, 1) circshift(data.stim, 1)];
        %   data.decision_pupil = circshift(data.decision_pupil, 1);
        data.(whichMod) = circshift(data.(whichMod), 1);
        
        data.resp(data.resp < 0) = 0;
        
        % split by those pupil values
        if nbins > 2,
            uncQs = quantile(data.(whichMod), nbins - 1);
        elseif nbins == 2,
            uncQs = median(data.(whichMod));
        end
        
        for u = 1:nbins,
            
            switch u
                case 1
                    trls = find(data.(whichMod) < uncQs(u));
                case nbins
                    trls = find(data.(whichMod) > uncQs(u-1));
                otherwise
                    trls = find(data.(whichMod) > uncQs(u-1) & data.(whichMod) < uncQs(u));
            end
            
            [b] = glmfit(designM(trls, :), data.resp(trls), ...
                'binomial','link','logit');
            
            % save betas
            grandavg.logistic(sj, u, :) = b;
            
        end
    end
    
    for u = 1:nbins,
        
        subplot(4,4,cnt); cnt = cnt + 1;
        hold on;
        plot([-1 1], [-1 1], 'color', 'k', 'linewidth', 0.5);
        plot([-1 1], [1 -1], 'color', 'k', 'linewidth', 0.5);
        
        load('~/Data/pupilUncertainty/GrandAverage/sjcolormap.mat');
        load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'plain'));
        
        for sj = 1:27,
            h = ploterr(grandavg.logistic(sj, u, 3), grandavg.logistic(sj, u, 4), ...
                [], ...
                [], '.', 'abshhxy', 0);
            set(h(1), 'color', mycolmap(sj, :), 'markerfacecolor', mycolmap(sj, :));
            %  set(h(2), 'color', mycolmap(sj, :), 'linewidth', 0.5);
            % set(h(3), 'color', mycolmap(sj, :), 'linewidth', 0.5);
        end
        
        % also show the mean
        h = ploterr(mean(grandavg.logistic(:, u, 3)), mean(grandavg.logistic(:, u, 4)), ...
            std(grandavg.logistic(:, u, 3)), ...
            std(grandavg.logistic(:, u, 4)), ...
            'o', 'abshhxy', 0);
        set(h(1), 'markeredgecolor', 'k', 'markerfacecolor', 'w', 'markersize', 4);
        set(h(2), 'color', 'k');
        set(h(3), 'color', 'k');
        
        % layout
        maxlim = 0.4;
        xlim([-maxlim maxlim]); ylim([-maxlim maxlim]);
        maxlim = 0.4;
        set(gca, 'xtick', -maxlim:maxlim:maxlim, 'ytick', -maxlim:maxlim:maxlim);
        
        if m == 2, xlabel('Choice weight'); end
        if u == 1, ylabel('Stimulus weight'); end
        box on; axis square;
        
        switch m
            case 1
                whichm = 'pupil';
            case 2
                whichm = 'RT';
        end
        
        switch u
            case 1
                title(['Low ' whichm]);
            case 2
                title(['Medium ' whichm]);
            case 3
                title(['High ' whichm]);
        end
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
    end
    cnt = cnt + 1;
end

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4c_decisionStrategies_byPupil.pdf'));
