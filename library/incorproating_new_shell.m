% get irreducible and primitive shells; this function has goals: 
%   1) save [v] orbit reps
%   2) save [s_ck] the space symmetries which stabilize the bond
%   3) save [g_ck] the space symmetries which generate the orbit
%   4) save [Qi]   the space symmetry   which maps the prototypical
%                  irreducible pair onto the rep [v]
%
% NOTE: requires p2u to point to atoms who's uc coordinates share a common
% closest primitive lattice point. set this in get_primitive_cell
% 
% NOTE: application of space symmerties requires the permutation matrices
% P1s and P2s to not have any zero values. all atoms must be mapped onto
% another atom.
%
% NOTE: The primitive shells as defined here involve only two primitive
% cell atoms for each given orbit. Although some space symmetries may map
% atoms onto other primitive cell atoms, those orbits are decomposed into
% two unique primitive shells which are related to the prototypical
% irreducible shell via Qi.



clear;clc

import am_lib.*

% load restart.mat

tiny = am_lib.tiny;

cutoff=5;
fname='POSCAR';
% fname='POSCAR.BMg2';
% fname='outfile.supercell';
% fname='infile.supercell';
flags='';

% get cells
[uc,pc,ic] = get_cells(fname,flags);

% get irreducible shells
[pp]         = get_pairs(pc,uc,cutoff)

% force constant model
% bvk = get_bvk_model(pc,pp);

% [bvk,pp] = get_bvk(cutoff,pc,uc,'infile.force_position.4.00-300',flags);

% bzp = get_bz_path(pc,'fcc');

% plot_bvk_dispersion(bvk,bzp);
%%
plot_bvk_dispersion(bvk,bzp,pp,uc) % MATCHES SYMBOLIC VERSION






