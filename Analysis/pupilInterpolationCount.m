% count

global mypath;

for sj = 1:27,
    
    load(sprintf('%s/Data/P%02d_alleyeNaN.mat', mypath, sj));
    
    percentageNaN = @(x) mean(isnan(x), 2);
    % for each trial, get the percentage of interpolated samples
    
    interpolatedSamples = cellfun(percentageNaN, data.trial, 'uniformoutput', 0);
    interpolatedSamples = cat(2, interpolatedSamples{:});
    interpolatedSamples = interpolatedSamples(find(strcmp(data.label, 'EyePupil')==1), :); % take the pupil chan
    
    % for each subject, compute the percentage of trials with > 50% rejected samples
    grandavg.rejectTrials(sj) = 100 * mean(interpolatedSamples > 0.5);
end
    
fprintf('range of trials with > 50% interpolation: %.3f to %.3f, mean %.3f \n', ...
    min(grandavg.rejectTrials), max(grandavg.rejectTrials), mean(grandavg.rejectTrials));