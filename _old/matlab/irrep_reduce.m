function [V,C] = irrep_reduce(rep,plot)
    % set tiny
    tiny=1E-6;

    % get dimensions
    nbases= size(rep,1);
    nsyms = size(rep,3);

    % initialize decomposition loop
    M = rep;
    V = eye(nbases);
    inds = ones(nbases,1);
    ninds = 1;

    % initialize plot
    if (plot)
    figure(1); set(gcf,'color','white'); clf; box on; clr = {'b','k','r','y','c','g'}; hold on;
    end
    
    % six loops is enough to decompose the irrep
    for k = 1:6
        % loop over cycle structures
        for j = 1:ninds
            mask = inds==j;

            H = irrep_get_H( M(mask,mask,:) );
    %         [Vp,~] = schur(H)
    
            if any(abs(H-H')>tiny) 
                fprintf('ERROR: H is not hermitian\n')
                return
            end
            [Vp,D] = eig(H);
            D = diag(D);

            Dunq = unique(D);
            for i = 1:length(Dunq)
                mask2 = abs(D-Dunq(i))<tiny;
                Vp(:,mask2) = orth( Vp(:,mask2) );
            end

            for i = 1:nsyms
                M(mask,mask,i) = Vp'*M(mask,mask,i)*Vp;
            end
            V(:,mask) = (V(:,mask)*Vp);
        end
            C=sum(abs(M),3);
            C(C<tiny)=0;
            C(C>tiny)=1;
            
            if (plot)
                spy(C,clr{k});
            end
%             k
%             C
            
            C = rref(C);
            ninds = sum(any(C~=0,2));
            for i = 1:ninds
                inds(C(i,:)~=0) = i;
            end
    end
    
    
end