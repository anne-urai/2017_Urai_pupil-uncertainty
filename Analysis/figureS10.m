% This code reproduces the analyses in the paper
% Urai AE, Braun A, Donner THD (2016) Pupil-linked arousal is driven 
% by decision uncertainty and alters serial choice bias. 
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% If you use the Software for your own research, cite the paper.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
% DEALINGS IN THE SOFTWARE.
%
% Anne Urai, 2016
% anne.urai@gmail.com

clearvars -except mypath; clc; close all;
global mypath;

subplot(441); sjCorrelation('pupil', 'correct');  %title('Correct');
ylim([-0.7 0.5]); set(gca, 'ytick', [-0.4 0 0.4]);
xlim([-0.5 0.5]); set(gca, 'xtick', [-0.4 0 0.4]);
ylabel('Pupil x correct weight');
xlabel('Choice weight');

subplot(442); sjCorrelation('pupil', 'incorrect'); %title('Incorrect');
ylim([-0.7 0.5]); set(gca, 'ytick', [-0.4 0 0.4]);
xlim([-0.5 0.5]); set(gca, 'xtick', [-0.4 0 0.4]);
ylabel('Pupil x error weight');
set(gca, 'yaxislocation', 'right');

% test between the two corr
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));
y1 = dat.('correct_pupil')(:, 1);
y2 = dat.('incorrect_pupil')(:, 1);
x = dat.response(:, 1);

% one way, parametric Steiger's test
[rddiff,cilohi,p] = rddiffci(corr(x, y1), corr(x, y2), corr(y1, y2), 27, 0.05);
[p, rddiff] = permtest_correlation(x, y1, y2);
xlabel(sprintf('dr = %.3f, p = %.3f \n', rddiff, p));

% now for RT
subplot(4,4,3);  sjCorrelation('rt', 'correct');
ylim([-0.7 0.5]); set(gca, 'ytick', [-0.4 0 0.4]);
xlim([-0.5 0.5]); set(gca, 'xtick', [-0.4 0 0.4]);
ylabel('RT x correct weight');
set(gca, 'yaxislocation', 'left');
xlabel('Choice weight');

subplot(4,4,4);  sjCorrelation('rt', 'incorrect');
ylim([-0.7 0.5]); set(gca, 'ytick', [-0.4 0 0.4]);
xlim([-0.5 0.5]); set(gca, 'xtick', [-0.4 0 0.4]);
ylabel('RT x error weight');
set(gca, 'yaxislocation', 'right');

y1 = dat.('correct_rt')(:, 1);
y2 = dat.('incorrect_rt')(:, 1);
x = dat.response(:, 1);
% one way, parametric Steiger's test
[rddiff,cilohi,p] = rddiffci(corr(x, y1), corr(x, y2), corr(y1, y2), 27, 0.05);
[p, rddiff] = permtest_correlation(x, y1, y2);
fprintf('correct rt vs incorrect rt, delta rho = %.3f, p = %.3f \n', rddiff, p);
xlabel(sprintf('dr = %.3f, p = %.3f \n', rddiff, p));

print(gcf, '-dpdf', sprintf('%s/Figures/FigureS10.pdf', mypath));
