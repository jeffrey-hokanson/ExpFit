% HSVD as described by HCDH94

function [omega,a]=hsvd_fast(y,omega0, L)

p = length(omega0);
n = length(y);

if nargin < 3
	L = ceil(n/2);
else
	L = max(ceil(n/2), L);
end
n = L*2-1;

Fc1 = hankel_prep(y(1:L), y(L:n));
Fc2 = hankel_prep(conj(y(1:L)), conj(y(L:n)));

if 0
	H = hankel(y(1:L),y(L:n));
	Z = zeros(size(H));
	B = [ Z, H'; H, Z];
	x = randn(2*L,1);
	y1 = B*x;
	y2 = mkB(x, Fc1, Fc2);
	norm(y1-y2)
end

% we drop the k subscript
%[U,S,V] = svds(H,L);
opts.isreal = 0;
opts.issym = 1;
%opts.v0 = y(1:2*);
opts.disp = 0;
opts.tol = 1e-20;
[V,E] = eigs(@(x) mkB(x, Fc1, Fc2), 2*L, p, 'lr', opts);
%ew = diag(E)
%non_zero = length(find(ew>0));
%keep = min(L, non_zero)
U = V(1:L,:);
[U,R] = qr(U,0);

% form the matrix containing U_k lower and U_k upper
Ub = U(1:L-1,1:p);
Ut = U(2:L,1:p);

%Z = Ub\Ut;
ut = U(end,1:p);	% \tilde{u}_{N-M+1} in the paper
% this uses the Sherman--Forrison inverse formula to simplify construction
A = Ub'*Ut;
Z = A + ut'*(ut*A)/(1-ut*ut');
%norm(Z1-Z,'fro')

% the eigenvalues of Z correspond to e^omega, so we take the log and then 
% move all the eigenvalues into the upper half plane
omega = log(conj(eig(Z)));
omega = real(omega) + 1i*mod(imag(omega),2*pi);

if nargout>1
	a = make_V(n,omega)\y(1:n);
end

function y = mkB(x, Fc1, Fc2)
% Return [0 H';H 0]*x where H is a Hankel matrix
n = length(x);
y1 = hankel_product(Fc2, x(n/2+1:end));
y2 = hankel_product(Fc1, x(1:n/2));
y = [y1;y2];



