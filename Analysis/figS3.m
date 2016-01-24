% figure 4 overview

lagGroups = 1;
mods = {'feedbackpupil', 'fb-decpupil'}; %, 'pupil-rt', 'feedbackpupil', 'fb-decpupil'};
nbins = 3;
close;
figure;
spcnt = 1;
for m = 1:length(mods),
    whichmodulator = mods{m};
    
    % middle row
    subplot(4,4,(m-1)*8+1); fig4c_psychFuncShift_Bias_byResp(lagGroups, whichmodulator, 'all', nbins); % todo
    set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
    
    subplot(4,4,(m-1)*8+2); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'all', nbins, []);
    % subplot(446); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'all', 0);
    set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
    
    subplot(4,4,(m-1)*8+3); fig4d_psychFuncShift_Slope(lagGroups, whichmodulator, 'all', [], nbins);
    %  subplot(447); fig4e_psychFuncShift_Slope(lagGroups, whichmodulator, 'all', 0);
    set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
    
    subplot(4,8,(m-1)*16+7); fig4f_HistoryPupil_Bar(lagGroups, whichmodulator, 'all');
    % set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
    
    switch m
        case 1
            suplabel('Post-feedback pupil response', 'x');
        case 2
            suplabel('Residual post-feedback - pre-feedback pupil responses', 'x');
    end
            
    
end

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty_paper/figureS3.pdf'));
