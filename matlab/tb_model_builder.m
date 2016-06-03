
clear;clc;
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

v=zeros(13,1);
v(1)=10;
v(2)=60;
v(3)=-90;
v(4)=5;
% v(6)=-8;
% v(7)= 10;
% v(8)= -5
v(9) = 10;
% v(10) = -10
% v(11) = 10

[pg] = load_tb_point_group();
[pc] = load_poscar('outfile.POSCAR.primitive');
[bz] = get_kpoint_path(pc.bas,kstart,kend,40);

for i = 1:bz.npaths
for j = 1:bz.ndivs
    H = getH(pg, v, bz.path(i).kpt_cart(:,j));
    bz.path(i).D(:,j) = sort(real(eig(H)));
end
end

for i = 1:bz.npaths
    plot(bz.path(i).x,bz.path(i).D,'-')
    hold on;
end
hold off;
axis tight;
box on;
