modes = {'projected', 'hsvd', 'hsvd_fast', 'full', 'mat-vec'};

N = 50

for j = 1:length(modes)
	mode = modes{j};
	load(sprintf('fig_timing_nmr_%s.mat', mode))

	min_time = min(times(1:N,:), [], 1);
	mean_time = mean(times(1:N,:), 1);
	median_time = median(times(1:N,:), 1);
	max_time = max(times(1:N,:), [], 1);
	pgf_dump(sprintf('fig_timing_nmr_%s.dat',mode),{'n','min','mean','median', 'max'},[n_vec',min_time', mean_time', median_time', max_time']);
end

