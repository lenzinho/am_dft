clear;clc

import am_lib.*

cutoff=5; dt=2; %ps 

% get cells
[uc,pc,ic] = get_cells('POSCAR'); 

% load md
[md] = load_md(uc,'infile.force_position',dt);

% get shells
[bvk,pp] = get_bvk(pc,uc,md,cutoff);
%%
% % plot dispersion aling high symmetry path
path={'hex','fcc-short','fcc'}; path=path{1};
plot_bvk_dispersion(bvk,get_bz_path(pc,31,path));

%%

% dt = 0.1; nsteps = 1000; Q = 10; T = 300; [md] = run_bvk_md(bvk,pp,uc,dt,nsteps,Q,T)
% 
% %%
% F = plot_md_cell(md);

%%


% n=[4,4,4]; kpt=[0;0;0]; amp=4; mode=4; nsteps=31;
% [dc,idc] = get_displaced_cell(pc,bvk,n,kpt,amp,mode,nsteps); 
% clf; F = plot_md_cell(idc); movie(F,3)


%%



nbands = 1878;
en_equi = load_band_energies(nbands,'infile.electron_energies_ref');
en_aimd = load_band_energies(nbands,'infile.electron_energies');

nEs = 101; E = linspace(-5,15,nEs);
degauss=1;
A_equi = permute(sum(gauss_((en_equi-permute(E(:),[3,2,1]))./degauss)./degauss,1),[3,2,1]);
A_aimd = permute(sum(gauss_((en_aimd-permute(E(:),[3,2,1]))./degauss)./degauss,1),[3,2,1]);
%%
n=[5;5;5];
fbz = get_fbz(pc,n);
fbz = get_bvk_dispersion(bvk,fbz);
q2u = get_bvk_normal_atransform(bvk,uc,fbz);
q_sk = q2u\reshape(md.tau,[],md.nsteps);


%%
surf(A_equi-A_aimd,'edgecolor','none')

% en_equi

%%

n=[4,4,4]; kpt=[1/2;1/2;1/2]; amp=4; mode=3; nsteps=51;
[dc,idc] = get_displaced_cell(pc,bvk,n,kpt,amp,mode,nsteps); 
% clf; F = plot_md_cell(idc); movie(F,3)

















