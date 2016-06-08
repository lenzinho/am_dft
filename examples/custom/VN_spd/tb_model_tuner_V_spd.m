
clear;clc;
tiny=1.0E-6;

% kpoint path [frac]
kstart=[
    0.500000   0.500000   0.500000
    0.000000   0.000000   0.000000
    1.000000   0.500000   0.500000
    ]';
kend=[
    0.000000   0.000000   0.000000
    0.000000   0.500000   0.500000
    0.000000   0.000000   0.000000
    ]';
klabel={'L','G','X','G'};

v(1) = 1;
v(2) = 2;
v(3) = 0.1;
v(4) = 0.5;
v(5) = 0.5;
v(6) = 3; % shell 3
v(7) = 4; % shell 3

[pg] = load_tb_point_group();
[pc] = load_poscar('outfile.primitive');
[bz] = get_kpoint_path(pc.bas,kstart,kend,40);

for i = 1:bz.npaths
for j = 1:bz.ndivs
%     H = get_H_numeric_cart(pg, v, bz.path(i).kpt_cart(:,j));
    H = get_H_model_frac(v, bz.path(i).kpt_frac(:,j)*2);
    if (abs(norm(H-H'))>tiny) 
        fprintf('H is not hermitian\n')
        return
    end
    bz.path(i).D(:,j) = sort(real(eig(H)));
end
end

for i = 1:bz.npaths
    plot(bz.path(i).x,bz.path(i).D,'-')
    hold on;
end
hold off; axis tight; box on;

% get ticks
ticks(1) = 0;
for i = 1:bz.npaths
    ticks(i+1) = bz.path(i).x(end);
end
set(gca,'XTick',ticks);
set(gca,'Xticklabel',klabel);