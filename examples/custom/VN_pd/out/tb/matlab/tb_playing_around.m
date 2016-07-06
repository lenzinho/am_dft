%% step 1: fresh start
clear;clc;
tiny=1.0E-6;
% step 2: oad symmetry rep
dumpdir = '/Volumes/Lenzinho/MyLinux/calc.vasp.5.3.3/development/am_lib/examples/custom/VN_pd';
load_dumps
% dumpdir = '/Volumes/Lenzinho/MyLinux/calc.vasp.5.3.3/development/am_lib/examples/custom/VN_pd/';
% load_dumps
%% step 3: set Fermi energy, load dispersion
EF = 9.50006808;
[dr,bz]=load_vasp_eigenval('../../../EIGENVAL');
dr.E=dr.E;%-EF;
%% step 4: load primitive cell, convert k-points to cart 
[pc]  = load_poscar('../../uc/outfile.primitive');
bz.kpt = pc.recbas*bz.kpt;
%% step 5: plot
X=[...
% +1.28     -0.49     -3.94     +1.01     +1.03     -0.12     -0.22     -0.08     -0.55     +0.05     +0.14     -0.13     +0.44     +0.30
 +10.78     +9.01     +5.56     +1.01     +1.03     -0.12     -0.22     -0.08     -0.55     +0.05     +0.14     -0.13     +0.44     +0.30
  +3.02     +3.47     +2.37     +0.13     +0.62     -0.31     -0.12     +0.06     +0.25     -0.20     +0.05     -0.05     +0.48     +0.10
  +3.67     +2.49     +2.91     +0.06     +0.53     -0.28     -0.06     -0.00     -0.04     -0.10     +0.26     -0.14     +0.35     +0.03
  +3.58     +2.66     +2.81     +0.18     +0.50     -0.26     -0.12     -0.14     +0.02     -0.13     +0.25     -0.16     +0.39     +0.10
  +3.58     +2.66     +2.81     +0.18     +0.50     -0.26     -0.12     -0.14     +0.02     -0.13     +0.25     -0.16     +0.39     +0.10
  +3.58     +2.66     +2.81     +0.18     +0.50     -0.26     -0.12     -0.14     +0.02     -0.13     +0.25     -0.16     +0.39     +0.10
  +3.61     +2.62     +2.84     +0.16     +0.51     -0.26     -0.12     -0.10     +0.02     -0.13     +0.25     -0.16     +0.38     +0.09
  +3.62     +2.60     +2.83     +0.12     +0.49     -0.25     -0.07     -0.03     -0.02     -0.18     +0.24     -0.17     +0.39     +0.09
  +3.59     +2.67     +2.80     +0.08     +0.49     -0.27     -0.13     -0.04     +0.02     -0.21     +0.24     -0.16     +0.40     +0.10
  +3.62     +2.61     +2.83     +0.12     +0.49     -0.25     -0.10     -0.03     -0.00     -0.19     +0.25     -0.17     +0.39     +0.09
  +3.62     +2.60     +2.83     +0.12     +0.49     -0.24     -0.08     -0.03     -0.02     -0.18     +0.25     -0.17     +0.39     +0.09];

%%

for j = 1:10

for i = 1:bz.nkpts
    H(:,:,i) = get_H_numeric_cart(tbpg, X(j,:), bz.kpt(:,i),[1:100]);
end
[V,D]=eigenshuffle(H);

clf;
plot(real(D'))
hold on;
plot(dr.E',':')
drawnow;
pause(1)
end
