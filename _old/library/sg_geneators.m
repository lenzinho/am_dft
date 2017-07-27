clc

% i=sg_possibilties(1);
k=225;

g = reshape(uniquecol_(reshape(S{k}(1:3,1:3,:),9,[])),3,3,[]);
G = get_multiplication_table(g); ngs = size(G,1);
h = reshape(uniquecol_(reshape(S_p(1:3,1:3,:),9,[])),3,3,[]);
H = get_multiplication_table(h); nhs = size(H,1);

for i = 1%:3
Gp = identify_point_symmetries(g);
Hp = identify_point_symmetries(h);

Gc = identify_classes(G); fwd=rankc_([Gp;Gc(:).']); rev(fwd)=[1:ngs]; G=rev(G(fwd,fwd)); Gp=Gp(fwd); Gc=Gc(fwd); 
Hc = identify_classes(H); fwd=rankc_([Hp;Hc(:).']); rev(fwd)=[1:nhs]; H=rev(H(fwd,fwd)); Hp=Hp(fwd); Hc=Hc(fwd);
end
[Gp;Hp]

[Gc,Hc].'

[sum([Gc].'==[1:10].',2),sum([Hc].'==[1:10].',2)]

% sasdf
% 
% [Gc,fwd(Gc(rev)).']
% asdf
% Hc = identify_classes(H); %Hc = sum(Hc(:).'==unique(Hc),2).'*(Hc(:).'==unique(Hc));
% asd
% fwd=rankc_([Gc(:).';Gp]); rev(fwd)=[1:ngs]; G(fwd,fwd)=fwd(G); Gc = Gc(fwd); Gp = Gp(fwd); g = g(:,:,fwd);
% fwd=rankc_([Hc(:).';Hp]); rev(fwd)=[1:nhs]; H=fwd(H(rev,rev)); Hc = Hc(fwd); Hp = Hp(fwd); h = h(:,:,fwd);

subplot(2,1,1); imagesc(G); daspect([1 1 1]);
subplot(2,1,2); imagesc(H); daspect([1 1 1]);

%%


V = sum(Gp==unique(Gp(:)),2); nVs=numel(V);

J = int8(perms(1:nVs)); E = cumsum(V); S = cumsum(V)-V+1;
F = factorial(V); EF = cumsum(F); SF = cumsum(F)-F+1;



%%



% try trivial: phi: G -> H
phi = zeros(ngs,1);
phi = [1:nhs];
% phi(1) = 1; 

self_consistent = true;
for i = 1:nhs; if self_consistent==true
for j = 1:nhs; if self_consistent==true
    if phi(g(i,j))==0
        phi(g(i,j)) = h(i,j);
    else
        % isomorphism is defined as phi: G -> H such that 
        % phi(ij) = phi(i)phi(j) for every i,j in G
        if phi(g(i,j)) ~= g(phi(i),phi(j))
            self_consistent=false;
            [i,j,h(i,j),g(i,j)]
        end
    end
end; end
end; end
sum(phi)
fprintf('%3i %3i %3i %3i %3i %3i %3i %3i %3i %3i %3i %3i %3i \n',phi); fprintf('\n');








aasdf

%%


clear;clc
%%
% clear;clc
% % generate databases using:
% import am_lib.*
for i = 1:237
%     % generate space symmetries and multiplication table
%     [S{i},MT{i}] = generate_sg(i); 
%     % determine pg corresponding to each sgi %4i , ...\n',cell2mat(pg))
    pg{i} = identify_pointgroup( reshape(uniquecol_( reshape(S{i}(1:3,1:3,:),9,[]) ),3,3,[]) );
    % count number of symmetries
    nsyms{i} = size(MT{i},1);
    % indentify number of classes
    nclasses{i} = numel(unique(identify_classes(MT{i})));
end
fprintf('pg_database = [ ... \n'); fprintf('%4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i , ...\n',cell2mat(pg)); fprintf(']; \n');
fprintf('nsyms_database = [ ... \n');fprintf('%4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i , ...\n',cell2mat(nsyms)); fprintf(']; \n');
fprintf('nclasses_database = [ ... \n');fprintf('%4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i , ...\n',cell2mat(nclasses)); fprintf(']; \n');

%%
[~,pc] = get_cells('../examples/custom/VN_pd/POSCAR');
% generate space symmetries and multiplication table
[T_p,H_p,S_p,R_p] = get_symmetries(pc); identify_spacegroup(identify_pointgroup(R_p));
%%

% get list of space groups with matching point group
possibility = identify_spacegroup(identify_pointgroup(R_p)); npossibilities=numel(possibility);
% get possible centering vectors (append identity, remove all zeros)
m=1; 
for i = 1:numel(possibility)
    A=uniquecol_(reshape(S{possibility(i)}(1:3,4,:),3,[])); n=size(A,2);
    T(1:3,m:(m+n-1))=A; m=m+n;
end
T=uniquecol_(T); T=T(:,sum(T,1)~=0);T=[T,eye(3)];nTs=size(T,2);
% choose three vectors at a time
ijk=nchoosek_(nTs,3);nijks=size(ijk,2);
% loop over each basis
for i = 1:npossibilities
for j = 1:nijks
    A = zeros(3,3,nsyms{possibility(i)});
for k = 1:nsyms{possibility(i)}
    if det(T(:,ijk(:,k)))
    A(:,:,j) = T(:,ijk(:,k)) * S{possibility(i)}(1:3,1:3,j) * inv(T(:,ijk(:,k)));
end
asdf
end
end


    


