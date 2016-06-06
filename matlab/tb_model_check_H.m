function [bool] = tb_model_check_H(pg,H)
    % [pg] = load_tb_point_group();
    % [pp] = load_tb_matrix_elements();
    % [H]  = get_H_explicit(pp,[0,0,0]);
    % defaults
    tiny=1.0E-6;
    bool = true;
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
        bool = false;
        fprintf('ERROR: %i symmetry operations found to not commute with H.\n',j)
    end
    % check H = H'
    if (abs(norm(H-H'))>tiny) 
        fprintf('H is not hermitian\n');
        bool = false;
    end
end