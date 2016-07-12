global mypath;
mypath = '/Users/anne/Data/pupilUncertainty_FigShare';

% add the code folder
addpath([mypath '/Code/Analysis/']);
cd([mypath '/Code/Analysis/']);
close all;

%%  mediation analysis: add model-based uncertainty to csv files
for sj = 1:27,
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % fit a probit slope so we can get the sigma
    b       = glmfit(data.motionstrength, (data.resp > 0), 'binomial', 'link', 'probit');
    sigma   = 1/b(2);  % standard deviation at these values is the inverse!
    bound   = -b(1); % the bound is negative if people say
    
    % for each trial, compute the average level of uncertainty
    data.uncertainty = arrayfun(@simulateUncertainty, abs(data.motionstrength), ...
        data.correct, sigma*ones(length(data.correct), 1), bound*ones(length(data.correct), 1));
    
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
    
    % save for later analyses
    writetable(data, sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
end

figure(1);
suplabel('Uncertainty', 'x');
suplabel('Pupil responses', 'y');
print(gcf, '-dpdf', sprintf('%s/Figures/uncertaintyPupil.pdf', mypath));

figure(2);
suplabel('Uncertainty', 'x');
suplabel('Reaction time', 'y');
print(gcf, '-dpdf', sprintf('%s/Figures/uncertaintyRT.pdf', mypath));

%% check that this looks the same as RT and pupil for uncertainty encoding
figure(3);
subplot(441); b = Uncertainty_byErrorCorrect('uncertainty');
ylabel('Uncertainty');

% now the repetition behaviour...
nbins = 4;
subplot(442); psychFuncShift_Bias('uncertainty', nbins, 1);
title('Correct'); ylim([0.45 0.6]);
subplot(443); psychFuncShift_Bias('uncertainty', nbins, 0);
title('Error'); ylim([0.45 0.6]);
subplot(444); b =  psychFuncShift_Bias('uncertainty', nbins, []);
title('All trls'); ylim([0.45 0.6]);

% correlate to mediation weights
load(sprintf('%s/Data/GrandAverage/mediationModel.mat', mypath));

subplot(445);
plot(b(:, end) - b(:, 1), grandavg.c_corr, '.'); axis tight; lsline;
xlabel('p(repeat) high - low'); ylabel('mediation c'); axis square;

print(gcf, '-dpdf', sprintf('%s/Figures/uncertaintyRepetition.pdf', mypath));

%% this looks super weird... redo for every individual!
close all; nbins = 4;

for sj = 1:27,
    clf;
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % first, establish the correlation between pupil, RT and uncertainty
    subplot(441); hold on;
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
    xlabel('Uncertainty'); ylabel('Pupil');
    
    subplot(442); hold on;
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
    xlabel('Uncertainty'); ylabel('RT');
    
    subplot(443); hold on;
    plot(data.decision_pupil(data.correct == 1), ...
        data.rt(data.correct == 1), 'g.');
    plot(data.decision_pupil(data.correct == 0), ...
        data.rt(data.correct == 0), 'r.');
    axis tight; axis square; box off;
    l = lsline;
    l(1).Color = l(1).Color * 0.8;
    l(2).Color = l(2).Color * 0.8;
    
    rho1 = corr(data.decision_pupil(data.correct == 1), ...
        data.rt(data.correct == 1), 'type', 'spearman');
    rho2 = corr(data.decision_pupil(data.correct == 0), ...
        data.rt(data.correct == 0), 'type', 'spearman');
    title(sprintf('%.2f / %.2f', rho1, rho2));
    xlabel('Pupil'); ylabel('RT');
    
    % then, show for all those 3
    subplot(445); hold on;
    psychFuncShift_Bias('uncertainty', nbins, 1, sj);
    psychFuncShift_Bias('uncertainty', nbins, 0, sj);
    set(gca, 'ytick', 0:0.05:1);
    ylims(1,:) = get(gca, 'ylim');
    
    subplot(446); hold on;
    psychFuncShift_Bias('rt', nbins, 1, sj);
    psychFuncShift_Bias('rt', nbins, 0, sj);
    set(gca, 'ytick', 0:0.05:1);
    ylims(2,:) = get(gca, 'ylim');
    
    subplot(447); hold on;
    psychFuncShift_Bias('pupil', nbins, 1, sj);
    psychFuncShift_Bias('pupil', nbins, 0, sj);
    set(gca, 'ytick', 0:0.05:1);
    ylims(3,:) = get(gca, 'ylim');
    
    % then, show for all those 3
    subplot(4,4,9); hold on;
    psychFuncShift_Bias('uncertainty', nbins, [], sj);
    set(gca, 'ytick', 0:0.05:1);
    ylims(4,:) = get(gca, 'ylim');
    
    subplot(4,4,10); hold on;
    psychFuncShift_Bias('rt', nbins, [], sj);
    set(gca, 'ytick', 0:0.05:1);
    ylims(5,:) = get(gca, 'ylim');
    
    subplot(4,4,11); hold on;
    psychFuncShift_Bias('pupil', nbins, [], sj);
    set(gca, 'ytick', 0:0.05:1);
    ylims(6,:) = get(gca, 'ylim');
    
    % make sure the axes are the same
    newylims(1) = min(ylims(:));
    newylims(2) = max(ylims(:));
    for i = [5:7 9:11],
        subplot(4,4,i); ylim(newylims);
    end
    
    print(gcf, '-dpdf', sprintf('%s/Figures/P%02d_uncertaintyRepetition.pdf', mypath, sj));
end
