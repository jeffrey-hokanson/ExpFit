% function V = make_V(n, omega)
%
% Returns the matrix V(omega) in C^(n \times |omega|) where
%
%   [ V(omega) ]_{j,k} = exp(j * omega_k)
% 
% n      - number of rows to create
% omega  - p parameters, each corresponding to one column
% V      - n x m output matrix
% 

function V = make_V(n, omega)
m = length(omega);

V = zeros(n,m);
for j = 1:m
	V(:,j) = exp(omega(j)*(0:n-1)');
end
