function [J,F1,F2] = tb_model_fitter_compute_jacobian(bz,dr,pg,x,mask)
    
    delta = 1E-2;
    delta_mat = eye(length(x))*delta;
    
    for i = 1:length(x)
        x1 = x;
        x2 = x(:) + delta_mat(:,i);
        F1(:,1) = tb_model_fitter_compute_residual(bz,dr,pg,x1,mask);
        F2(:,1) = tb_model_fitter_compute_residual(bz,dr,pg,x2,mask);
         J(:,i) = (F2-F1)/delta;
    end
end