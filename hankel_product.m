
function HX = hankel_product(Fc,X);
n = length(Fc);
m = size(X,1);
X = [ flipud(X);zeros(n-m,size(X,2))];


%HX = (ifft( diag(Fc)*fft(X)));
HX = ifft(bsxfun(@times,Fc,fft(X)));
HX = (HX(1:m,:));
