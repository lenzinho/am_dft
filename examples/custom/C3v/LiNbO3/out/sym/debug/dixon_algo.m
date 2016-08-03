clear;clc;
load_dumps
%% select irrep
j = 9
S=cumsum([0, pg.ct.irrep_dim(1:(pg.ct.nirreps-1)).^2])+1;
E=cumsum(pg.ct.irrep_dim.^2);

clear rr_block
for i = 1:pg.nsyms
    rr_block(:,:,i) = pg.ct.irrep_proj_V'*pg.mt.rr(:,:,i)*pg.ct.irrep_proj_V;
end

for j = 1:pg.ct.nirreps
    H = irrep_get_H(rr_block(S(j):E(j),S(j):E(j),:));
    [V,~] = eig(H);
    for i = 1:pg.nsyms
        T = (pg.ct.irrep_proj_V(:,S(j):E(j)) * V);
        rr_block(S(j):E(j),S(j):E(j),i) = T' * pg.mt.rr(:,:,i) * T;
    end
end
%%
rr_block(abs(rr_block)<0.00001)=0;
%%
spy_tensor(rr_block)



