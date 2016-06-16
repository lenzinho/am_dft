function V = cycle_structure2eigenvector(cycle_structure)
    i2pi = sqrt(-1)*2*pi;
    tiny = 1E-6;
    nbases = length(cycle_structure);
    nonzeros = sum(cycle_structure~=0);
    V(1:nbases,1) = 0;
    j = 0;
    for i = 1:nbases
        if (abs(cycle_structure(i)>tiny))
            j = j + 1;
            V(i) = exp(i2pi*(j-1)/nonzeros);
        end
    end
end