

subplot(441);
FruendKernels('pupil+rt', 'response_pupil')
print(gcf, '-dpdf', sprintf('%s/Figures/pupilLags.pdf', mypath));
