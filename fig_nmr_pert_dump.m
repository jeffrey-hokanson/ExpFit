%modes = {'full','reduced','hsvd', 'hsvd_fast'}
modes = {'projected', 'hsvd_fast', 'full'}


f = [-86;-70;-54;152;168;292;308;360;440;490;530];
d = [50;50;50;50;50;50;50;25;285.7;25;200];
amp = [75;150;75;150;150;150;150;150;1400;60;500];
th = 135*ones(11,1);
n0 = 256;
dt0 = 1e-3*1/3; % DIANA's correction

%n = 1024;
n = 4096;
%n = 256;

dt = dt0*(n0/n);

omega_hat = (2i*pi*f-d)*dt;
omega_hat = real(omega_hat)+1i*mod(imag(omega_hat),2*pi);
a_hat = amp.*exp(1i*th*pi/180);

J = [make_Vp(n, omega_hat)*diag(a_hat), make_V(n,omega_hat)];
RIJ = [real(J),-imag(J); imag(J), real(J)];

[Q,Gamma] = qr(RIJ,0);

epsilon = 15;


figure(1),clf
x = linspace(0,100,500)';
y = chi2pdf(x,44);
plot(x,y,'k')
hold on

pgf_dump(sprintf('fig_nmr_pert_chi2.dat'),{'x','y'},[x,y]);

colors = {'b','r','g','k'}

for j = 1:length(modes)
	mode = modes{j};
	load(sprintf('fig_nmr_pert_%d_%s.mat',n, mode))
	err = zeros(length(times),1);
	theta_hat = [omega_hat;a_hat];
	for k = 1:length(times)
		[p,II] = marriage_norm(omega_hat,omega_vec(:,k));
		theta = [omega_vec(II,k);a_vec(II,k)];
		err(k) = norm(Gamma*([real(theta-theta_hat);imag(theta-theta_hat)]))^2;
	end

	%[den, x] = density_est_lin_1D(err/epsilon^2*2, [0,100],40)
	[bw, den, x] = kde(err/epsilon^2*2, 2^8, 0, 250,0.06);
	x = reshape(x,length(x),1);
	den = reshape(den,length(x),1);
	
	pgf_dump(sprintf('fig_nmr_pert_%d_%s.dat',n, mode),{'x','y'},[x,den]);
	plot(x,den, colors{j})
end


