clear;clc;
%%
load_dumps
%%
B = tbpg.ct.block_proj;
% B = orth(rand(9));
% B = eye(size(tbpg.sym,1));


for i = 1:size(tbpg.sym,3)
A(:,:,i) = B' * tbpg.sym(:,:,i) * B;
end

[V,C] = irrep_reduce(A,true)

%%
V = sg.ct.block_proj;
for i = 1:size(sg.sym,3)
A(:,:,i) = V' * sg.sym(:,:,i) * V;
end

spy(sum(abs(A),3)>1E-14)


%% recreate tbpg
for i = 1:tbpg.nsyms
    rr_block(:,:,i) = tbpg.ct.irrep_proj_V'*tbpg.sym(:,:,i)*tbpg.ct.irrep_proj_V;
end
phi(1:tbpg.nbases,1:tbpg.nbases) = 0;
for j = 1:tbpg.ct.nirreps
    mask = tbpg.ct.irrep_id==j;
    H = irrep_get_H( rr_block(mask,mask,:) );
    [V,D] = eig(H);
    phi(:,mask) = (tbpg.ct.irrep_proj_V(:,mask) * V);
end
for i = 1:tbpg.nsyms
    rr_block(:,:,i) = phi(:,:)' * tbpg.sym(:,:,i) * phi(:,:);
end
rr_block(abs(rr_block)<0.00001)=0;
spy_tensor(rr_block)
%% check tbpg
clear rr_block
for i = 1:tbpg.nsyms
    rr_block(:,:,i) = tbpg.ct.block_proj' * tbpg.sym(:,:,i) * tbpg.ct.block_proj;
end
rr_block(abs(rr_block)<0.00001)=0;
spy_tensor(rr_block)





%% recreate pg
for i = 1:pg.nsyms
    rr_block(:,:,i) = pg.ct.irrep_proj_V'*pg.sym(:,:,i)*pg.ct.irrep_proj_V;
end
phi(1:pg.nbases,1:pg.nbases) = 0;

for j = 1:pg.ct.nirreps
    mask = pg.ct.irrep_id==j;
    H = irrep_get_H( rr_block(mask,mask,:) );
    [V,D] = eig(H);
    phi(:,mask) = (pg.ct.irrep_proj_V(:,mask) * V);
end
for i = 1:pg.nsyms
    rr_block(:,:,i) = phi(:,:)' * pg.sym(:,:,i) * phi(:,:);
end
rr_block(abs(rr_block)<0.00001)=0;
spy_tensor(rr_block)
%% check pg
clear rr_block
for i = 1:pg.nsyms
    rr_block(:,:,i) = pg.ct.block_proj' * pg.sym(:,:,i) * pg.ct.block_proj;
end
rr_block(abs(rr_block)<0.00001)=0;
spy_tensor(rr_block)




%% recreate sg
for i = 1:sg.nsyms
    rr_block(:,:,i) = sg.ct.irrep_proj_V'*sg.sym(:,:,i)*sg.ct.irrep_proj_V;
end
phi(1:sg.nbases,1:sg.nbases) = 0;
for j = 1:sg.ct.nirreps
    mask = sg.ct.irrep_id==j;
    H = irrep_get_H( rr_block(mask,mask,:) );
    [V,D] = eig(H);
    phi(:,mask) = (sg.ct.irrep_proj_V(:,mask) * V);
    
    
    if j==3
        H
        V
        D
        
        asdf
    end
end
for i = 1:sg.nsyms
    rr_block(:,:,i) = phi(:,:)' * sg.sym(:,:,i) * phi(:,:);
end
rr_block(abs(rr_block)<0.00001)=0;
spy_tensor(rr_block)
%% check sg
clear rr_block
for i = 1:sg.nsyms
    rr_block(:,:,i) = sg.ct.block_proj' * sg.sym(:,:,i) * sg.ct.block_proj;
end
rr_block(abs(rr_block)<0.00001)=0;
spy_tensor(rr_block)




