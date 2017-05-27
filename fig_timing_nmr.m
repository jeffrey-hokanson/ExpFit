function fig_timing_nmr(mode, Niter, rerun)

f = [-86;-70;-54;152;168;292;308;360;440;490;530];
d = [50;50;50;50;50;50;50;25;285.7;25;200];
amp = [75;150;75;150;150;150;150;150;1400;60;500];
th = 135*ones(11,1);
n0 = 256;
dt0 = 1e-3*1/3; % DIANA's correction

p = 11;

if strcmp(mode,'projected')
	n_vec = 2.^([8:24]);
elseif strcmp(mode, 'hsvd')
	n_vec = 2.^([8:14]);
elseif strcmp(mode, 'hsvd_fast')
	n_vec = 2.^([8:23]);
elseif strcmp(mode, 'full')
	n_vec = 2.^([8:20]);
elseif strcmp(mode, 'mat-vec')
	n_vec = 2.^([8:24]);
end

if nargin < 3
	rerun = 0;
end

try 
	if rerun == 0
		load(sprintf('fig_timing_nmr_%s.mat', mode))
	else
		throw()
	end
catch
	times = zeros(0,length(n_vec));
	omega_vec = zeros(p,0,length(n_vec));
	a_vec = zeros(p,0, length(n_vec));
end


% How many iterations conducted previously
it = size(times,1)
mode
opts = struct();
opts.tolX = 1e-10;
opts.tolFun = 1e-9;
%opts.optimization_iter_per_efficiency_check = 100;



while it < Niter
	it = it+1;
	for jsize = 1:length(n_vec)
		n = n_vec(jsize);
		% Set random number generators
		try
			s = RandStream('mcg16807','Seed',it*n);
			RandStream.setDefaultStream(s);
		catch
			s = RandStream('mcg16807','Seed',it*n);
			RandStream.setGlobalStream(s);
		end

		dt = dt0*(n0/n);
		
		omega_hat = (2i*pi*f-d)*dt;
		omega_hat = real(omega_hat)+1i*mod(imag(omega_hat),2*pi);
		a_hat = amp.*exp(1i*th*pi/180);
		y_hat = make_V(n, omega_hat)*a_hat;

		g = (randn(n,1)+randn(n,1))/sqrt(2);
		
		% The noise in the original exmaple of VHB97
		y = y_hat + 15*g;

		if strcmp(mode, 'projected')
			tic,
			[omega, a] = projected_expfit(y, p, opts);
		elseif strcmp(mode, 'full')
			tic,
			[omega, a] = expfit(y, p, opts);
		elseif strcmp(mode, 'hsvd')
			tic,
			[omega,a] = hsvd(y,ones(p,1));
		elseif strcmp(mode, 'hsvd_fast')
			tic,
			[omega,a] = hsvd_fast(y,ones(p,1), 100);
		elseif strcmp(mode, 'mat-vec')
			omega = omega_hat;
			a = a_hat;
			tic,
			make_V(n, (1i-1/n)*ones(44,1))'*y;
		end
		
		times(it, jsize) = toc;

		omega_vec(:,it,jsize) = omega;
		a_vec(:, it, jsize) = a;
		[err,I] = marriage_norm(omega,omega_hat);
		%[omega(I), omega_hat]

		fprintf('Finished %s iteration %5d for n=%8d in %5f seconds with error=%5g\n',mode,it, n, times(it,jsize),err);
	end
	
	save(sprintf('fig_timing_nmr_%s.mat', mode),'times','omega_vec','a_vec', 'n_vec')
end
fprintf('\n')

