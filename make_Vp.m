% function Vp = mkVp(omega,n)
% 
% Constructs the matrix Vp(omega) in C^{n \times |omega|} where
%
%   [ Vp(omega) ]_{j,k} = j * exp( j * omega_k )
% 
% n      - number of rows to create
% omega  - m parameters, each corresponding to one column
% Vp     - n x m output matrix

function V = make_Vp(n, omega)
p = length(omega);

sc = log( (1:n-1)');
V = zeros(n,p);
for j = 1:p
	% actual formula should read:
	% V(2:end,j) = (0:n-1)'.*exp(omega(j)*(1:n-1)'.*sc);
	% To avoid possible numerical instability, we pull the scaling 
	% into the exponential.
	V(2:end,j) = exp(omega(j)*(1:n-1)' + sc);
end
