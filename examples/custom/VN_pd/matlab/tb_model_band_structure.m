
clear;clc;
tiny=1.0E-6;

% % kpoint path [frac]
% kstart=[
%     0.000000   0.000000   1.000000
%     0.500000   0.500000   1.000000
%     0.000000   0.000000   0.000000
%     ]';
% kend=[
%     0.000000   0.500000   0.500000
%     0.000000   0.000000   0.000000
%     0.000000   0.500000   0.000000
%     ]';
% klabel={'G','X','G','L'};

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
v(1:11) = 0;
v=[...
-1.96E+01
-1.98E+01
-1.99E+01
-2.00E+01
-2.01E+01
-2.01E+01
-2.00E+01
-2.00E+01
-1.99E+01
-1.98E+01];

% tb=xml_read('save.tb');
% reshape(str2num(tb.group.sym.value),tb.group.sym.shape);

[pc] = load_poscar('../POSCAR');
[bz] = get_kpoint_path(pc.bas,kstart,kend,40);

%%

if (~tb_check_hamiltonian_at_gamma(pg)) 
    fprintf('WARNING: HAMILTONIAN DOES NOT COMMUTE WITH POINT SYMMETRIES AT GAMMA!\n')
end

for i = 1:bz.npaths
for j = 1:bz.ndivs
    H = getH(pg, v, bz.path(i).kpt_cart(:,j));
    if (abs(norm(H-H'))>tiny) 
        fprintf('H is not hermitian\n')
        bz.path(i).kpt_cart(:,j)
        H
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