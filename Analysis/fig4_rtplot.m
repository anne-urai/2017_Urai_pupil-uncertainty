
close all; clear;
mods = {'pupil', 'rt'};
for m = 1:2,
    subplot(4,4,m);
    fig4_psychFuncShift_RT(mods{m}, 'all', 1);
    fig4_psychFuncShift_RT(mods{m}, 'all', 0);
end

mods = {'pupil', 'rt'}; m = 2;
%for m = 1:2,
subplot(4,4,3);

fig4_psychFuncShift_RT(mods{m}, 'switch', []);
subplot(4,4,4);

fig4_psychFuncShift_RT(mods{m}, 'repeat', []);
%end


print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty_paper/fig4_RTadvantage.pdf'));
