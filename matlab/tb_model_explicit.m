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

[pg] = load_tb_point_group();
[pp] = load_tb_matrix_elements();
[pc] = load_poscar('outfile.primitive');
[bz] = get_kpoint_path(pc.bas,kstart,kend,40);

for i = 1:bz.npaths
for j = 1:bz.ndivs
    H = get_H_explicit(pp, bz.path(i).kpt_cart(:,j));l
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



