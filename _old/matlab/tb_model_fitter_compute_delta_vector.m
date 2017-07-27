function [d,F1,A,B] = tb_model_fitter_compute_delta_vector(bz,dr,pg,x,mask)
    [J,F1,~] = tb_model_fitter_compute_jacobian(bz,dr,pg,x,mask);
    A = transpose(J)*J;
    B = transpose(J)*F1;
    d(:,1) = A\B;
end