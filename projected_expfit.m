function [omega, a] = projected_expfit(y, p, opts)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default options 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default_opts = struct();
% Number of steps to take before updating subspace
default_opts.opt_steps = 10;
% Convergence tolerance for movement of parameter values
default_opts.tolX = 1e-16;
% Convergence tolerance for change in objective function
default_opts.tolFun = 1e-16;
% If 'on', display per iteration results
default_opts.display = 'off';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options for the optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
	opts = struct();
end
fields = fieldnames(default_opts);
for i = 1:length(fields)
	if ~isfield(opts, fields{i})
		opts = setfield(opts, fields{i}, getfield(default_opts, fields{i}));
	end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configure optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt_lsqnonlin = optimoptions(@lsqnonlin, 'Jacobian','on',...
		 'Display', opts.display, ...
		 'MaxIterations', opts.opt_steps, ...
		 'TolX', opts.tolX, ...
		 'TolFun', opts.tolFun);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega = [];
a = [];
n = length(y);
freqs = (0:n-1)'/n*2*pi;
mu = [];
Vmu_y = [];
R = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup for subspace 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = [-inf; -1.424175; -6.718913e-1; -3.537494e-1; -1.818821e-1; -9.206280e-2];
n_rings = floor(log2(n))+1;
if length(alpha) >= n_rings
	alpha = alpha(1:n_rings);
else
	alpha = [alpha; -2.9720*(0.5.^(length(alpha):n_rings-1)')];
end
alpha(end) = 0;

for p_current = 1:p 

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Estimate new frequency
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	append_omega();

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Main optimization loop for fixed set of frequencies
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	exitflag = -1;
	subspace_unchanged = -1;
	it = 0;
	while ~ ( (exitflag > 0) & (subspace_unchanged) )
		update_subspace();
		
		[omega, ~, ~, exitflag] = lsqnonlin(@(x) objfun(x), omega, [], [], opt_lsqnonlin);
		% Update iteration counter and terminate if we've taken too many steps
		it = it + 1;
		if it > 20; break; end
	end
end

	function update_subspace()
		mu_new = [];
		for k = 1:length(omega)
			ell = max(find(real(omega(k)) < alpha));
			if isempty(ell); ell = length(alpha); end;
			ell = min(length(alpha) -1, ell);	% can't define outside of unit circle
			
			if ell == 0
				mu_new = [mu_new; alpha(1); alpha(2); alpha(2) + 1i*pi];
			else
				% Inner ring always has points added at roots of 2^ell
				kk = floor(imag(omega(k))*2^ell/(2*pi));
				mu_new = [mu_new; alpha(ell) + 2i*pi*kk/2^ell; alpha(ell) + 2i*pi*(kk+1)/2^ell];

				if ell == length(alpha) - 1
					% Outer ring is origin ring
					kk = floor(imag(omega(k))*n/(2*pi));
					mu_new = [mu_new; alpha(ell+1) + 2i*pi*kk/n; alpha(ell+1) + 2i*pi*(kk+1)/n];
				else	
					kk = floor(imag(omega(k))*2^(ell+1)/(2*pi));
					mu_new = [mu_new; alpha(ell+1) + 2i*pi*kk/2^(ell+1); alpha(ell+1) + 2i*pi*(kk+1)/2^(ell+1)];
				end	
			end
		end
		% Remove repeated points
		mu_new = unique(mu_new);
		% Remove points that we already have
		mu_new = setdiff(mu_new, mu);

		if ~ isempty(mu_new)
			% Update basis
			mu = [mu; mu_new];
			% Compute new inner-products
			Vmu_y = [Vmu_y; (make_V(n,mu_new)'*y)];	

			% Update norm-correction matrix
			VmuVmu = make_VV(n, mu, mu);
		
			[U, ew] = eig(VmuVmu);
			ew = diag(ew);
			I = find(ew>1e-7);
			R = U(:,I)*diag(1./sqrt(ew(I)));

			%norm(R'*VmuVmu*R - eye(length(mu)), 'inf')		% should be zero
			%svd(make_V(n, mu)*R) 		% should be a vector of ones
			subspace_unchanged = 0;
		else
			subspace_unchanged = 1;
		end

	end % update_subspace

	function append_omega()
		% Compute residual
		if length(omega) == 0
			r = -y;
		else
			r = make_V(n, omega)*a - y;
		end	
		Fr = fft(r);
		[trash,idx] = max(abs(Fr));
		idx = idx(1);
		% add new exponential
		omega = [omega;1i*freqs(idx)];
	end % append_omega
	
	% This nested function has access to variables inside the parent scope
	function [r, J] = objfun(omega)
		% Force the real part not to be too far into the right half plane
		omega = min(real(omega), 10/n) + 1i * imag(omega);
		if nargout == 1
			VV = make_VV(n, mu, omega);
		else	
			[VV, VVp] = make_VV(n, mu, omega);
		end

		WV = R'*VV;
		Wy = R'*Vmu_y;

		[Q, T] = qr(WV, 0);
		a_ = T\(Q'*Wy); %a_ = V\y;
		r = WV*a_ - Wy;

		a = a_; % Pass amplitude out to global scope
		if nargout < 2; return; end;
		
		% Compute the Jacobian
		WVp = R'*VVp;
		WVpa = WVp * diag(a_);
		L = WVpa - Q*(Q'*WVpa);
		K = -Q*((T')\diag(WVp'*r));
		J = L + K;
		return
	end % objfun
end % expfit




