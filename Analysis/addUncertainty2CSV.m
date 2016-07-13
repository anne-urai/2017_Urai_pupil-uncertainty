
%%  mediation analysis: add model-based uncertainty to csv files
plotme = 0;

for sj = 1:27,
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % fit a probit slope so we can get the sigma
    b       = glmfit(data.motionstrength, (data.resp > 0), 'binomial', 'link', 'probit');
    sigma   = 1/b(2);  % standard deviation at these values is the inverse!
    bound   = -b(1); % the bound is negative if people say
    
    % for each trial, compute the average level of uncertainty
    data.uncertainty = arrayfun(@simulateUncertainty, abs(data.motionstrength), ...
        data.correct, sigma*ones(length(data.correct), 1), bound*ones(length(data.correct), 1));
    
    if plotme,
        figure(1);
        subplot(5,6,sj); hold on;
        plot(data.uncertainty(data.correct == 1), ...
            data.decision_pupil(data.correct == 1), 'g.');
        plot(data.uncertainty(data.correct == 0), ...
            data.decision_pupil(data.correct == 0), 'r.');
        axis tight; axis square; box off;
        l = lsline;
        l(1).Color = l(1).Color * 0.8;
        l(2).Color = l(2).Color * 0.8;
        
        rho1 = corr(data.uncertainty(data.correct == 1), ...
            data.decision_pupil(data.correct == 1), 'type', 'spearman');
        rho2 = corr(data.uncertainty(data.correct == 0), ...
            data.decision_pupil(data.correct == 0), 'type', 'spearman');
        title(sprintf('%.2f / %.2f', rho1, rho2));
        
        figure(2)
        subplot(5,6,sj); hold on;
        plot(data.uncertainty(data.correct == 1), ...
            data.rt(data.correct == 1), 'g.');
        plot(data.uncertainty(data.correct == 0), ...
            data.rt(data.correct == 0), 'r.');
        axis tight; axis square; box off;
        l = lsline;
        l(1).Color = l(1).Color * 0.8;
        l(2).Color = l(2).Color * 0.8;
        
        rho1 = corr(data.uncertainty(data.correct == 1), ...
            data.rt(data.correct == 1), 'type', 'spearman');
        rho2 = corr(data.uncertainty(data.correct == 0), ...
            data.rt(data.correct == 0), 'type', 'spearman');
        title(sprintf('%.2f / %.2f', rho1, rho2));
    end
    
    % save for later analyses
    writetable(data, sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
end

if plotme,
    figure(1);
    suplabel('Uncertainty', 'x');
    suplabel('Pupil responses', 'y');
    print(gcf, '-dpdf', sprintf('%s/Figures/uncertaintyPupil.pdf', mypath));
    
    figure(2);
    suplabel('Uncertainty', 'x');
    suplabel('Reaction time', 'y');
    print(gcf, '-dpdf', sprintf('%s/Figures/uncertaintyRT.pdf', mypath));
end
