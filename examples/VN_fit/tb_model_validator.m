
clear;clc;clf;
load('EeV.mat')
tiny=1.0E-6;

% kpoint path [frac]
kstart=[
    0.000000   0.000000   1.000000
    0.500000   0.500000   1.000000
    0.000000   0.000000   0.000000
    ]';
kend=[
    0.000000   0.500000   0.500000
    0.000000   0.000000   0.000000
    0.000000   0.500000   0.000000
    ]';

v=[ 1 1 10 ];

[pg] = load_tb_point_group();
[pc] = load_poscar('outfile.supercell');
[bz] = get_kpoint_path(pc.bas,kstart,kend,40);

for i = 1:bz.npaths
for j = 1:bz.ndivs
    H = getH(pg, v, bz.path(i).kpt_cart(:,j));
    bz.path(i).D(:,j) = sort(real(eig(H)));
end
end

for i = 1:bz.npaths
    plot([1:length(bz.path(i).x)]+40*(i-1),bz.path(i).D,'r')
    hold on;
end
hold off;
axis tight;
box on;
%

hold on;
plot(mod(1:size(EeV,1),121),EeV,'.','color',[1,1,1]*0.65)