function V_perm = permute_nonzero_entries(V,degree)
    tiny = 1E-6;
    nbases = length(V);
	inds = [1:nbases];
    inds_zero    = inds( abs(V)<tiny );
    inds_nonzero = inds( abs(V)>tiny );
    % apply shift
    inds_nonzero_perm = circshift(inds_nonzero,degree,2);
    inds(inds_nonzero) = inds(inds_nonzero_perm);
    %
    V_perm = V(inds);
    %
end