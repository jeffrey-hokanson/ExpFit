function [omega, a] = expfit(y, p, opts)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default options 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default_opts = struct();
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
		 'TolX', opts.tolX, ...
		 'TolFun', opts.tolFun ...
		 ); %'CheckGradients', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega = [];
a = [];
n = length(y);
freqs = (0:n-1)'/n*2*pi;

for p_current = 1:p 

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Estimate new frequency
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Main optimization loop for fixed set of frequencies
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[omega, trash, trash, exitflag] = lsqnonlin(@(x) objfun(x), omega,[],[],opt_lsqnonlin);
end


	
	% This nested function has access to variables inside the parent scope
	function [r, J] = objfun(omega)
		% Force the real part not to be too far into the right half plane
		omega = min(real(omega), 10/n) + 1i * imag(omega);
		
		V = make_V(n, omega);
		[Q, T] = qr(V, 0);
		%a_ = V\y;
		a_ = T\(Q'*y);
		r = V*a_ - y;

		a = a_; % Pass amplitude out to global scope
		if nargout < 2; return; end;

		Vp = make_Vp(n, omega);
		
		% Compute the Jacobian
		b = Q'*y;
		Py = Q*b;
		Vpa = Vp * diag(a_);
		L = Vpa - Q*(Q'*Vpa);
		K = -Q*((T')\diag(Vp'*r));
		J = L + K;
		
		return
	end % objfun


end % expfit




