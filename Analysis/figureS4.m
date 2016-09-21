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
close; figure;

% first, illustrate shift
sigmoid = @(x, mu) 1 ./ (1 + exp(-x + mu));
x = -5:0.01:5;
subplot(441); % tiny!
hold on;
plot(x, sigmoid(x, 0), 'k');
plot([-5 5], [0.5 0.5], 'k');
plot([0 0], [0 1], 'k');
axis tight; axis square;
ylabel('P(choice = 1)');
xlabel('Sensory evidence (a.u.)');
ylim([-0.05 1]); 

subplot(442); examplePsychFunc(20); title(sprintf('Participant %d', 20));
subplot(443); examplePsychFunc(13); title(sprintf('Participant %d', 13));
subplot(444); examplePsychFuncShift(10); title('Participant 10');

print(gcf, '-dpdf', sprintf('%s/Figures/FigureS4.pdf', mypath));
