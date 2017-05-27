% If x and y are two vectors, return the precomputed information for a fast-hankel vector product
% (give in same order as hankel)

% Currently only works with square Hankel Matrices
function Fc = hankel_prep(x,y)

x = reshape(x,length(x),1);
y = reshape(y,length(y),1);
% first column in circulant matrix
c = [ y; x(1:end-1)];
Fc = fft(c);
