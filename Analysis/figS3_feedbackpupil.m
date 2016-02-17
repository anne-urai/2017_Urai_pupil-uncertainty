% figure 4 overview
global mypath;

lagGroups = 1;
mods = {'fbpupil', 'fb+decpupil'}; %, 'pupil-rt', 'feedbackpupil', 'fb-decpupil'};
nbins = 3;
close;
figure;
spcnt = 1;
for m = 1:length(mods),
    whichmodulator = mods{m};
    
    % middle row
    subplot(4,4,(m-1)*8+1); fig3c_psychFuncShift_Bias_byResp(whichmodulator, nbins); % todo
    set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
    
    subplot(4,4,(m-1)*8+2); fig3d_psychFuncShift_Bias_Slope(whichmodulator, nbins, []);
    set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
    
    %subplot(4,4,(m-1)*8+3); fig3d_psychFuncShift_Slope(whichmodulator, [], nbins);
   % set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
    
    subplot(4,8,(m-1)*16+7); fig3hi_HistoryPupil_Bar(whichmodulator);
    % set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
    
    switch m
        case 1
            suplabel('Post-feedback pupil response', 'x');
        case 2
            suplabel('Decision pupil included as covariate', 'x');
    end
    
end

print(gcf, '-dpdf', sprintf('%s/Figures/figureS2.pdf', mypath));
