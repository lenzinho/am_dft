clear;clc;
tiny=1.0E-6;
%
[pg] = load_tb_point_group();
[pp] = load_tb_matrix_elements();

kpt = [0,0,0];
[H] = get_tb_hamiltonian(pp,kpt);

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
fprintf('ERROR: %i symmetry operations found to not commute with H.\n',j)
