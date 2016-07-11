%% 

% group split, also correct and error
subplot(4,6, 15); MedianSplit('pupil', 'correct');
ylim([-0.35 0.2]);
ylabel('Pupil * correct');
subplot(4,6,16); MedianSplit('pupil', 'incorrect');
set(gca, 'yaxislocation', 'right');
ylabel('Pupil * error');ylim([-0.35 0.2]);

% group split, also correct and error
subplot(4,6, 17); MedianSplit('rt', 'correct');
ylim([-0.35 0.2]);
ylabel('RT * correct');
subplot(4,6,18); MedianSplit('rt', 'incorrect');
set(gca, 'yaxislocation', 'right');
ylabel('RT * error');ylim([-0.35 0.2]);


close;
subplot(5,5,1); SjCorrelation('pupil', 'correct');  title('Correct');
ylim([-0.5 0.5]); set(gca, 'ytick', [-0.4 0 0.4]); 
subplot(5,5,2); SjCorrelation('pupil', 'incorrect'); title('Incorrect');
ylim([-0.5 0.5]); set(gca, 'ytick', [-0.4 0 0.4]); 

% test between the two corr
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));
y1 = dat.('correct_pupil')(:, 1);
y2 = dat.('incorrect_pupil')(:, 1);
x = dat.response(:, 1);
% one way, parametric Steiger's test
[rddiff,cilohi,p] = rddiffci(corr(x, y1), corr(x, y2), corr(y1, y2), 27, 0.05);
fprintf('correct pupil vs incorrect pupil, delta rho = %.3f, p = %.3f \n', rddiff, p);

subplot(5,5,6);  SjCorrelation('rt', 'correct'); ylim([-0.5 0.5]); set(gca, 'ytick', [-0.4 0 0.4]); 
subplot(5,5,7);  SjCorrelation('rt', 'incorrect'); ylim([-0.5 0.5]); set(gca, 'ytick', [-0.4 0 0.4]); 

y1 = dat.('correct_rt')(:, 1);
y2 = dat.('incorrect_rt')(:, 1);
x = dat.response(:, 1);
% one way, parametric Steiger's test
[rddiff,cilohi,p] = rddiffci(corr(x, y1), corr(x, y2), corr(y1, y2), 27, 0.05);
fprintf('correct rt vs incorrect rt, delta rho = %.3f, p = %.3f \n', rddiff, p);

print(gcf, '-dpdf', sprintf('%s/Figures/figureS_scatterIndividual.pdf', mypath));
