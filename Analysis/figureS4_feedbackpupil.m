% figure 4 overview
global mypath;

lagGroups = 1;
mods = {'fbpupil', 'fb+decpupil'}; 
nbins = 3;
close;
figure;

for m = 1:length(mods),
    whichmodulator = mods{m};
    subplot(4,4,(m-1)*8+1); fig3c_psychFuncShift_Bias_byResp(whichmodulator, nbins); 
    subplot(4,4,(m-1)*8+2); fig3d_psychFuncShift_Bias_Slope(whichmodulator, nbins, []);
    subplot(4,8,(m-1)*16+6); fig3hi_HistoryPupil_Bar(whichmodulator);
end

print(gcf, '-dpdf', sprintf('%s/Figures/figureS3.pdf', mypath));
