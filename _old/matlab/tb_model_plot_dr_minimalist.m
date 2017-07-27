clear;clc


odir='/Volumes/Lenzinho/MyLinux/calc.vasp.5.3.3/development/tmn/LaN';

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

[v]  = load_tb_matrix_elements([odir,'/out/tb/outfile.tb_matrix_elements_irreducible']);
[pc] = load_poscar([odir,'/vasp/POSCAR']);
[pg.sym] = load_dump([odir,'/out/tb/outfile.tb_pg.sym']);
[bz] = get_kpoint_path(pc.bas,kstart,kend,40);

for i = 1:bz.npaths
    for j = 1:bz.ndivs
        H = get_H_numeric_frac(pg, v, bz.path(i).kpt_frac(:,j),[1:100]);
        bz.path(i).D(:,j) = sort(real(eig(H)));
    end
end

%plot
for i = 1:bz.npaths
    plot(bz.path(i).x,bz.path(i).D,'-')
    hold on;
end
hold off; axis tight; box on;
%ticks
ticks(1) = 0;
for i = 1:bz.npaths
    ticks(i+1) = bz.path(i).x(end);
end
set(gca,'XTick',ticks);
set(gca,'Xticklabel',klabel);