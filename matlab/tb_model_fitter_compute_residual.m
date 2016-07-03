function [R] = tb_model_fitter_compute_residual(bz,dr,pg,x,selector_kpoint,selector_shell,band_skip)
    for j = 1:bz.nkpts
        if (any(selector_kpoint==j))
            % get hamiltonian
            H = get_H_numeric_cart(pg, x, bz.kpt(:,j), selector_shell);
            % get eigenvalues
            D = real(eig(H));
            % get cost matrix for munkres assigment
            % c = abs(dr.E(1:pg.nbases,j)*ones(1,pg.nbases) - (D*ones(1,pg.nbases))');
            % evaluate assignemtn
            % [assign,~] = munkres(c);
            [~,assign] = sort(D);
            % get residual vector
            inds =  [1:pg.nbases] + (j-1)*pg.nbases;
            R(inds,1) = sort(dr.E([1:pg.nbases]+band_skip,j)) - D(assign);
        end
    end
end