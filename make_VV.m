function [VV, VVp] = make_VV(n, mu, omega)
%
% Computes the entires of V(mu)'*V(omega) where V(mu) in C^(n\times |mu|)
% is a Vandermonde matrix with entries
%
%  [ V(mu) ]_{j,k} = e^{ j * mu_k}
%
% If a second output parameter is requested, the product V(mu)'* Vp(omega) 
% is computed where
%
%  [ Vp(omega) ]_{j,k} = j * e^{ j * omega_k}
%
% If only mu is provided, it sets omega = mu

% Initial setup
if nargin < 3; omega = mu; end

p = length(omega);
m = length(mu);

mu = reshape(mu, m,1);
omega = reshape(omega, p, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the product V(mu)'* V(omega)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta = bsxfun(@plus, conj(mu), omega.');
ndelta = n*delta;
expm1_delta = expm1(delta);
expm1_ndelta = expm1(ndelta);
VV = expm1_ndelta./expm1_delta;

I = (abs(expm1_delta) < 1e-30);
VV(I) = n*(1 + (n-1)/2*delta(I));

% Handle infinity properly
I = (isinf(mu));
VV(I,:) = 1;
I = (isinf(omega));
VV(:,I) = 1;


% stop computation if we don't need the second term
if nargout <= 1; return; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the product V(mu)'* Vp(omega)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default case, exp(delta) > exp(0.25/n^2)
VVp = n*exp(n*delta)./expm1_delta - exp(delta).*expm1_ndelta./expm1_delta.^2;

% Second case, 0 < exp(delta) < exp(0.25/n^2)
%I = find(abs(exp(delta(:))) < exp(0.25/n^2));
I = (sqrt( real(delta(:)).^2 + (mod(imag(delta(:)), 2*pi)).^2) < 0.25/n^2);

coeff = (abs(delta(I))> 1e-30) .* expm1_ndelta(I)./expm1_delta(I)...
		 + (abs(delta(I))< 1e-30) .* n .* (1 + (n-1).*delta(I)./2.0);

ndelta = ndelta(I);
ndelta2 = ndelta.*ndelta;
ndelta3 = ndelta2.*ndelta;
ndelta5 = ndelta3.*ndelta2;
ndelta7 = ndelta5.*ndelta2;
ndelta9 = ndelta7.*ndelta2;
ndelta11 = ndelta9.*ndelta2;

delta2 = delta(I).*delta(I);
delta3 = delta2.*delta(I);
delta5 = delta3.*delta2;
delta7 = delta5.*delta2;
delta9 = delta7.*delta2;
delta11 = delta9.*delta2;


VVp_alt = coeff.*(...
			n * (0.5 + ndelta/12.0 - ndelta3/720.0 + ndelta5/30240.0 - ndelta7/1209600.0 + ndelta9/47900160.0 - ndelta11*691.0/1307674368000.0)...
			-(0.5 + delta(I)/12.0 - delta3/720.0 + delta5/30240.0 - delta7/1209600.0 + delta9/47900160.0 - delta11*691.0/1307674368000.0)...
			);

VVp(I) = VVp_alt;

% Third special cases exp(delta) = 1
I = (abs(exp(delta(:)) - 1) < 1e-13);
VVp(I) = n * (n - 1)/2.0;

% Handle infinity properly
I = (isinf(mu));
VVp(I,:) = 0;
I = (isinf(omega));
VVp(:,I) = 0;

