function fig_nmr_pert(algorithm, Niter, rerun)

if nargin < 3
	rerun = 0;
end

f = [-86;-70;-54;152;168;292;308;360;440;490;530];
d = [50;50;50;50;50;50;50;25;285.7;25;200];
amp = [75;150;75;150;150;150;150;150;1400;60;500];
th = 135*ones(11,1);
n0 = 256;
dt0 = 1e-3*1/3; % DIANA's correction

p = 11;

%n = 2^10;
%n = 2^8;
n = 2^12

fname = sprintf('fig_nmr_pert_%d_%s.mat', n, algorithm);

try 
	if rerun == 0
		load(fname)
	else
		throw()
	end
catch
	times = zeros(0,length(n));
	omega_vec = zeros(p,0,length(n));
	a_vec = zeros(p,0, length(n));
end


% How many iterations conducted previously
it = size(times,1)
algorithm

opts = struct();
opts.tolX = 1e-10;
opts.tolFun = 1e-9;

while it < Niter
	it = it+1;
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

	if strcmp(algorithm, 'projected')
		tic,
		[omega, a] = projected_expfit(y, p, opts);
	elseif strcmp(algorithm, 'full')
		tic,
		[omega, a] = expfit(y, p, opts);
	elseif strcmp(algorithm, 'hsvd')
		tic,
		[omega,a] = hsvd(y,ones(p,1));
	elseif strcmp(algorithm, 'hsvd_fast')
		tic,
		[omega,a] = hsvd_fast(y,ones(p,1), 100);
	elseif strcmp(algorithm, 'naive')
		tic,
		[omega, a] = expfit_naive(y, omega_hat, a_hat);
	end
	
	times(it, 1) = toc;

	omega_vec(:,it,1) = omega;
	a_vec(:, it, 1) = a;
	[err,I] = marriage_norm(omega,omega_hat);
	%[omega(I), omega_hat]

	fprintf('Finished %s iteration %5d for n=%8d in %5f seconds with error=%5g\n',algorithm,it, n, times(it,1),err);
	save(fname,'times','omega_vec','a_vec', 'n')
end
fprintf('\n')

