clear;clc;
tiny=1.0E-6;

[pg] = load_tb_point_group();
[pp] = load_tb_matrix_elements();
[H]  = get_tb_hamiltonian(pp,[0,0,0]);

% check H = H'
if (abs(norm(H-H'))>tiny)
    fprintf('ERROR: H is not hermitian!\n')
end

% check [H,R] = 0
j=0;
for i = 1:pg.nsyms
    R = pg.sym(:,:,i);
    if (abs(norm(H*R-R*H))>tiny)
        j=j+1;
        fprintf('ERROR: [H,R] ~= 0!\n')
    end
end
if (j~=0)
    fprintf('ERROR: %i symmetry operations found to not commute with H.\n',j)
end


%%

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
    H = get_tb_hamiltonian(pp, bz.path(i).kpt_cart(:,j));
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
hold off;
axis tight;
box on;



