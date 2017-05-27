%Niter = 4000;
%fig_timing_nmr_pert('full', Niter);
%fig_timing_nmr_pert('projected', Niter);
%fig_timing_nmr_pert('hsvd', Niter);
%fig_timing_nmr_pert('hsvd_fast', Niter);


Niter = 50;
%fig_timing_nmr('projected', Niter);
%fig_timing_nmr('hsvd_fast', Niter);
fig_timing_nmr('mat-vec', Niter)
%fig_timing_nmr('full', Niter);
%fig_timing_nmr('hsvd', Niter);

fig_timing_nmr_dump
