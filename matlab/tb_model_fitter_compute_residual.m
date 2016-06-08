function [R] = tb_model_fitter_compute_residual(bz,dr,pg,x,mask)
for j = 1:bz.nkpts
    if mask(j)
        H = get_H_numeric_cart(pg, x, bz.kpt(:,j));
        inds = (j*pg.nbases):((j+1)*pg.nbases-1);
        R(inds) = dr.E(1:pg.nbases,j) - sort(real(eig(H)));
    end
end
end

