
global mypath;

% grey shades
gr1 = [0.2 0.2 0.2];
gr2 = [0.6 0.6 0.6];

% sigma from datas
data = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
b = glmfit(data.motionstrength, (data.resp > 0), 'binomial', 'link', 'probit');

% standard deviation at these values is the inverse!
sigma = 1/b(2);

% our stimuli range from -6 to 6
stim = -7:0.01:20;

subplot(441);
% simulate two motionstrength items
s1 = normpdf(stim, 7, sigma) ./ max(normpdf(stim, 7, sigma));
p = plot(stim, s1, zeros(1, 2), [0 1], 'k');
axis off;
saveas(gcf, sprintf('~/Dropbox/Meetings/CertainNormal.eps'), 'epsc');

subplot(441);
% our stimuli range from -6 to 6

% simulate two motionstrength items
s1 = normpdf(stim, 7, sigma) ./ max(normpdf(stim, 7, sigma));
s2 = normpdf(stim, 7, sigma*2) ./ max(normpdf(stim, 7, sigma*2));
p = plot(stim, s1, stim, s2, zeros(1, 2), [0 max(s2)], 'k');
axis off;
saveas(gcf, sprintf('~/Dropbox/Meetings/UncertainNormal.eps'), 'epsc');

