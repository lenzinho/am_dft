clear;clc;
tiny=1.0E-6;
%
[pp] = load_tb_matrix_elements();
[pg] = load_tb_point_group();

kpt = [0,0,0];
[H] = get_tb_hamiltonian(pp,kpt);

% check [H,R] = 0
for i = 1:pg.nsyms
    R = pg.sym(:,:,i);
    if (abs(norm(H*R-R*H))>tiny)
       asdfs
    end
end

