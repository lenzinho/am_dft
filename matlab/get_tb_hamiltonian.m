function [H] = get_tb_hamiltonian(pp,kpt)
    minS = min(pp.S);
    maxE = max(pp.E);
    H(minS:maxE,minS:maxE) = 0;
    i2pi = 2*pi*sqrt(-1);
    for k = 1:pp.nshells
    SEm = pp.S(pp.shell(k).m):pp.E(pp.shell(k).m);
    SEn = pp.S(pp.shell(k).n):pp.E(pp.shell(k).n);
    for p = 1:pp.shell(k).npairs
        H(SEm, SEn) = H(SEm, SEn) + pp.shell(k).pair(p).V * exp(i2pi*dot(pp.shell(k).pair(p).tau,kpt));
    end
    end
end



