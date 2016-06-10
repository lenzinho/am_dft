function [R] = tb_model_fitter_compute_residual(bz,dr,pg,x,mask)
    for j = 1:bz.nkpts
        if (mask(j))
            % get hamiltonian
            H = get_H_numeric_cart(pg, x, bz.kpt(:,j));
            % get eigenvalues
            D = real(eig(H));
            % get cost matrix for munkres assigment
            % c = abs(dr.E(1:pg.nbases,j)*ones(1,pg.nbases) - (D*ones(1,pg.nbases))');
            % evaluate assignemtn
            % [assign,~] = munkres(c);
            [~,assign] = sort(D);
            % get residual vector
            inds =  [1:pg.nbases] + (j-1)*pg.nbases;
            R(inds,1) = sort(dr.E(1:pg.nbases,j)) - D(assign);
        end
    end
end