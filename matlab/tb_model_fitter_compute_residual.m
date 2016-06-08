function [R] = tb_model_fitter_compute_residual(bz,dr,pg,x0,x,xinds,mask,isplot)
    
    x_internal        = x0;
    x_internal(xinds) = x;
    
    for j = 1:bz.nkpts
        if mask(j)
            H = get_H_numeric_cart(pg, x_internal, bz.kpt(:,j));
            inds = (j*pg.nbases):((j+1)*pg.nbases-1);
            R(inds) = dr.E(1:pg.nbases,j) - sort(real(eig(H)));
        end
    end
    if (isplot)
        for j = 1:bz.nkpts
            H = get_H_numeric_cart(pg, x_internal, bz.kpt(:,j));
            tbdr.E(:,j) = sort(real(eig(H)));
        end
        plot(1:bz.nkpts,dr.E);
        hold on;
        plot(1:bz.nkpts,tbdr.E,'--');
        hold off;
        drawnow;
    end
end