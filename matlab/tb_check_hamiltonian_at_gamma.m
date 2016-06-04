function [pass] = tb_check_hamiltonian_at_gamma(pg)
    tiny=1.0E-6;
    % initialize random variables
    v=rand(100,1);
    % get hamiltonian
    H = getH(pg,v,[0,0,0]);
    % initialize pass output
    pass = true;
    % check commutator
    j=0;
    for i = 1:pg.nsyms
        R = pg.sym(:,:,i);
        if (abs(norm(H*R-R*H))>tiny)
            j=j+1;
            fprintf('ERROR: [H,R] ~= 0!\n')
            pass = false
        end
    end
    if (j~=0)
        fprintf('ERROR: %i symmetry operations found to not commute with H.\n',j)
    end
end