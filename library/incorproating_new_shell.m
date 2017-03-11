% get irreducible and primitive shells; this function has goals: 
%   1) save [v] orbit reps
%   2) save [s_ck] the space symmetries which stabilize the bond
%   3) save [g_ck] the space symmetries which generate the orbit
%   4) save [Qi]   the space symmetry   which maps the prototypical
%                  irreducible pair onto the rep [v]




clear all;clc

import am_lib.*

% load restart.mat

tiny = am_lib.tiny;

cutoff=3;
fname='POSCAR';
% fname='POSCAR.BMg2';
% fname='outfile.supercell';
fname='infile.supercell';
flags='';

% get cells
[uc,pc,ic] = get_cells(fname,flags);

% get irreducible shells
[bvk,pp] = get_bvk(cutoff,pc,uc,'infile.force_position.4.00-300.short',flags);

% %% numerical
bzp = get_bz_path(pc,'fcc');
plot_bvk_dispersion(bvk,bzp)

%%

clear;clc; import am_lib.*

cutoff=2.5;
fname='POSCAR';
% fname='POSCAR.BMg2';
% fname='outfile.supercell';
fname='infile.supercell';
flags='';

% get cells
[uc,pc,ic] = get_cells(fname,flags);

% get tight-binding
spdf={'d','p'}; nskips=5; fname='EIGENVAL.ibz'; Ef=0;
[tb,pp] = get_tb(pc,uc,cutoff,spdf,nskips,Ef,fname,flags)

% plot electron band structure along path
bzp = get_tb_dispersion(tb,bzp);

% define figure properties
fig_ = @(h)       set(h,'color','white');
axs_ = @(h,qt,ql) set(h,'Box','on','XTick',qt,'Xticklabel',ql);

figure(1); fig_(gcf); plot(bzp.x,sort(real(bzp.E)),'-k');
axs_(gca,bzp.qt,bzp.ql); axis tight; ylabel('Energy E'); xlabel('Wavevector k');




