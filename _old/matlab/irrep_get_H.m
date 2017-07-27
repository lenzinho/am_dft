function H = irrep_get_H(rr)
    tiny = 1E-6;
    nbases=size(rr,1);
    nsyms =size(rr,3);

    for r = 1:nbases
    for s = 1:nbases
        Hrs(1:nbases,1:nbases) = 0;
        if     r==s
            Hrs(r,s) = 1;
        elseif r>s
            Hrs(r,s) = 1;
            Hrs(s,r) = 1;
        elseif r<s
            Hrs(r,s) = sqrt(-1);
            Hrs(s,r) =-sqrt(-1);
        end
        
        H(1:nbases,1:nbases) = 0;
        for i = 1:nsyms
        H = H + transpose(conj(rr(:,:,i))) * Hrs * rr(:,:,i);
        end
        H = H / nsyms;
        
        if (any(abs((H(1,1)*eye(nbases) - H))>tiny))
            return
        end
    end
    end
    
    H = eye(nbases);
end