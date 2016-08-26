function ExamplePsychFuncShift(sj)

global mypath;

data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
data = data(find(data.sessionnr > 1), :);
% outcome vector need to be 0 1 for logistic regression
data.resp(data.resp == -1) = 0;

% condition on two previous choices
resps = [0 1];
for r = 1:2,
    trls = find(data.resp == resps(r));
    lag = 1;
    
    laggedtrls = trls+lag;
    % exclude trials at the end of the block
    if any(laggedtrls > size(data, 1)),
        trls(laggedtrls > size(data, 1)) = [];
        laggedtrls(laggedtrls > size(data, 1)) = [];
    end
    
    % remove trials that dont match in block nr
    removeTrls = data.blocknr(laggedtrls) ~= data.blocknr(trls);
    laggedtrls(removeTrls) = [];
    trls(removeTrls) = [];
    
    % fit logistic regression
    thisdat = data(laggedtrls, :);
    
    % get an overall logistic fit
    b = glmfit(thisdat.motionstrength, thisdat.resp, ...
        'binomial','link','logit');
    xvals = -4:0.1:4;
    psychCurve(r, :) = glmval(b, xvals, 'logit');
end

%% plot these two in the colors we will also use later on

colors = cbrewer('qual', 'Set2', 8);
colors = colors([1 3], :);
hold on;

plot([0 0], [0 1], 'k', 'linewidth', 0.1);
plot([min(xvals) max(xvals)], [0.5 0.5], 'k', 'linewidth', 0.1);

for r = 1:2, 
p(r) = plot(xvals, psychCurve(r,:), '-', 'color', colors(r,:), 'linewidth', 1);
end

xlim([min(xvals) max(xvals)]); set(gca, 'xtick', [min(xvals) 0 max(xvals)]);
ylim([-0.05 1]); set(gca, 'ytick', [0 0.5 1]);
xlabel('Sensory evidence (a.u.)');
ylabel('P(choice = 1)');
axis square;
box off;
%legend(p, {'previous choice -1'; 'previous choice 1'});
%legend boxoff;

end