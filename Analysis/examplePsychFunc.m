function examplePsychFunc(sj)
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

global mypath;

data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
data = data(find(data.sessionnr > 1), :);
% outcome vector need to be 0 1 for logistic regression
data.resp(data.resp == -1) = 0;

[bias, slope, lapse] = fitLogistic(data.motionstrength, data.resp);
xvals = -5:0.1:5;
Logistic = @(p, x) p(3)+(1-p(3)-p(3)) * (1./(1+exp(-p(2).*(x+p(1)))));
psychCurve = Logistic([bias slope lapse], xvals);
[binnedx, binnedy] = divideintobins(data.motionstrength, data.resp, 15);

%% plot these two in the colors we will also use later on
hold on;
plot([0 0], [0 1], 'k', 'linewidth', 0.1);
plot([min(xvals) max(xvals)], [0.5 0.5], 'k', 'linewidth', 0.1);
plot(xvals, psychCurve, '-', 'color', 'k', 'linewidth', 1);
plot(binnedx, binnedy, '.k', 'markersize', 12);

xlim([min(xvals) max(xvals)]); set(gca, 'xtick', [min(xvals) 0 max(xvals)]);
ylim([-0.05 1]); set(gca, 'ytick', [0 0.5 1]);
xlabel('Sensory evidence (a.u.)');
ylabel('P(choice = 1)');
axis square;
box off;

end