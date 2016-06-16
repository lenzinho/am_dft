function [U] = orthonormalize(V)
    m = size(V,1);
    n = size(V,2);
    U(1:m,1:n) = 0;
    for k = 1:n
        U(:,k) = V(:,k);
        for j = 1:(k-1)
            U(:,k) = U(:,k) - dot(U(:,j),V(:,k)) .* U(:,j);
        end
        U(:,k) = U(:,k)/norm(U(:,k));
    end
end