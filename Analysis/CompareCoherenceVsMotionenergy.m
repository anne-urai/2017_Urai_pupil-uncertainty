%% load in the same data but once with motion energy and one with coherence - why are they so different in filesize?

sj = 1;
coh = load(sprintf('~/Data/UvA_pupil/sim_backup/2ifc_sj%02d.txtresults.mat', sj));
men = load(sprintf('~/Data/UvA_pupil/sim_backup/2ifc_me_sj%02d.txtresults.mat', sj));

%% compare the weights that have been fit
clf; hold on;
plot(coh.model_w_hist.w(8:end));
plot(men.model_w_hist.w(8:end));
% these seem to look the same, have been corrected for slope


%% how about the bootstrap? is it okay if we just correct this with the
% slope?
men.bootstrap_corr(:,1:14)  = bsxfun(@times, men.bootstrap(:, 1:14), men.bootstrap(:, end-1));
coh.bootstrap_corr(:,1:14)  = bsxfun(@times, coh.bootstrap(:, 1:14), coh.bootstrap(:, end-1));
clf; hold on;
plot(men.bootstrap_corr(:, 1)); hold on; plot( coh.bootstrap_corr(:, 1))

%% histograms dont quite overlap... the motionenergy one has a higher value!
clf;
histogram(coh.bootstrap_corr(:, 2));
hold on;
histogram(men.bootstrap_corr(:, 2));
% this seems OK too, if we correct with the slope.

%% now, the files are much much bigger in motionenergy because of the permutation matrix. what happened here?

