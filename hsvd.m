% HSVD as described by HCDH94

function [omega,a]=hsvd(y,omega0,L)

m = length(omega0);
n = length(y);

if nargin <3
	%L = n -m +1;
	L = floor(n/2);
else
	if(L<m)
		fprintf('L parameter too small\n');
	end
	if(L>n-m+1)
		fprintf('L parameter too large\n');
	end
end

H = hankel(y(1:L),y(L:n));

% we drop the k subscript
[U,S,V] = svd(H);

% form the matrix containing U_k lower and U_k upper
Ub =  U(1:L-1,1:m);
Ut = U(2:L,1:m);


ut = U(end,1:m);	% \tilde{u}_{N-M+1} in the paper
% this uses the Sherman--Forrison inverse formula to simplify construction

A = Ub'*Ut;
Z = A + ut'*(ut*A)/(1-ut*ut');

% the eigenvalues of Z correspond to e^omega, so we take the log and then 
% move all the eigenvalues into the upper half plane
omega = log(eig(Z));
omega = real(omega) + 1i*mod(imag(omega),2*pi);


if nargout>1
	a = make_V(n,omega)\y;
end
