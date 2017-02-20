% get irreducible shells (shells for each irreducible atom)
% goals:
%   1) save [v] orbit reps (flip the bond, ordering primitive atom indices, if necessary)
%   2) save [s_ck] the space symmetries which stabilize the bond 
%   3) save [g_ck] the space symmetries which generate the orbit
%
% NOTE: requires p2u to point to atoms who's uc coordinates share a common
% closest primitive lattice point. set this in get_primitive_cell
% NOTE: application of space symmerties requires the permutation matrices
% P1s and P2s to not have any zero values. all atoms must be mapped onto
% another atom.

clear;clc

import am_lib.*


tiny = am_lib.tiny;

cutoff=5;
fname='POSCAR';
fname='POSCAR.BMg2';
% fname='outfile.supercell';
flags='';

[uc,~,ic] = get_cells(fname,flags);




% readjust cutoff based on unitcell
cutoff = min([normc_(uc.bas)/2,cutoff]);

% get conversion from primitive to supercell basis
uc2pc = ic.bas/uc.bas; pc2uc=inv(uc2pc);

% define function to apply symmetries to position vectors
sym_apply_ = @(S,tau) reshape(matmul_(S(1:3,1:3,:),tau),3,[],size(S,3)) + S(1:3,4,:);

% define function to change basis of seitz symmetry
sym_rebase_ = @(B,S) [[ matmul_(matmul_(B,S(1:3,1:3,:)),inv(B)), reshape(matmul_(B,S(1:3,4,:)),3,[],size(S,3))]; S(4,1:4,:)];

% define function to get the clostest primitive lattice vector in supercell basis
G_ = @(tau) tau - reshape(matmul_(uc2pc,mod_(matmul_(pc2uc,tau))),size(tau));

% define function to get bond vector from unit cell atom pair indicies
diff_ = @(V) uc.tau(:,V(2,:))-uc.tau(:,V(1,:));

% get point symmetries in [supercell basis]
S = ic.S; nsyms = size(S,3); S = sym_rebase_(uc2pc,S);

% determine primitive atoms shells which are connected
P1s=[];P2s=[]; npcs=numel(uc.p2u); npairs=zeros(1,npcs);
for i = 1:npcs
    % identify primitive atom in unit cell
    Delta = uc.tau(:,uc.p2u(i)); 

    % compute cutoff distances, exclude atoms above cutoff, count pairs involving the i-th primitive atom
    d=normc_(uc2ws(uc.bas*(uc.tau-Delta),uc.bas)); ex_=[d<cutoff]; npairs(i)=sum(ex_);

    % first atom in pair (must be in primitive cell); record G which shifts tau back to itself 
    tau = sym_apply_(S,Delta); G=G_(tau); P1 = member_(mod_(tau-G),uc.tau);
    
    % second atom in pair (SEE NOTEs)
    tau = sym_apply_(S,uc.tau(:,ex_));    P2 = member_(mod_(tau-G),uc.tau);

    % save results
    P1s = [P1s;repelem(P1,npairs(i),1)]; P2s = [P2s;P2]; %#ok<AGROW>
end

% define function to simultaneously apply operation (helps to not forget about one)
bundle_ = @(ex_,A,B) deal(A(:,ex_),B(ex_));

% pair vector and indicies
V = [ P1s(:), P2s(:) ].'; Vinds = [1:(sum(npairs)*nsyms)];

% check that all atoms are mapped
if any(V(:)==0,1); error('something is wrong; space symmetry must map atoms onto atoms without exceptions!'); end

% create a unique pair label
[~,i2p,p2i]=unique(sort(V).','rows','stable');

% get permutation representation
PM = zeros(sum(npairs),nsyms); PM(Vinds) = p2i;
% CHECK: PM(all(diff(PM,1,2)==0,2),:), V(:,PM(all(diff(PM,1,2)==0,2),1))

% get map connecting pairs (sort A based on distances)
A = get_binary_rep(PM); A=A( rankcol_( normc_(uc2ws(uc.bas*diff_( V(:,i2p(findrow_(A))) ),uc.bas)) ), :);
% CHECK: spy(A(:,~all(A==0,1)))
% CHECK: diff_ = @(V) uc.tau(:,V(2,:))-uc.tau(:,V(1,:));
% CHECK: normc_(uc2ws(uc.bas*diff_( V(:,i2p(A(30,:))) ),uc.bas))

Z=[];
for i = 1:ic.natoms
    Y=[];
    for p2i = 1:size(A,1)
        [U,Uinds] = bundle_(i2p(A(p2i,:)),V,Vinds); ic_id = ic.i2u(i); ex_ = any(U==ic_id,1);
        if any(ex_)
           rep = find(ex_,1);  % an arbitrary pair in the orbit
            xy = U(:,rep);     % uc indicies
            mn = uc.u2p(xy).'; % pc indicies
            ij = uc.u2i(xy).'; % ic indicies
            r  =  diff_(xy); r = uc2ws(uc.bas*r,uc.bas); % bond vector
            d  =   normc_(r); % bond length
            w  =    sum(ex_); % number of pairs in orbit


    %             s_ck = PM(find(A(2,:),1),:)== ;
    %  PM(find(A(2,:),1),:) == PM(find(A(p2i,:),1),24)
            if (ij(1)~=ic.u2i(ic_id))
                Y = [Y,[d(:);-r(:);w(:);ij([2,1]);mn([2,1]);-p2i]];
            else
                Y = [Y,[d(:);+r(:);w(:);ij([1,2]);mn([1,2]);+p2i]];
            end
        end
    end

    bar_ = @(x) repmat('-',[1,x]);
    fprintf('irreducible atom %i at [%.5f,%.5f,%.5f] has %i shells\n', i,ic.bas*ic.tau(:,i),size(Y,2));
    fprintf('%-10s   %-30s   %-4s   %-7s   %-7s   %-4s\n', 'd [cart]','bond [cart]','#','ic(i,j)','pc(m,n)','flp?'); 
    fprintf('%-10s   %-30s   %-4s   %-7s   %-7s   %-4s\n', bar_(10),bar_(30),bar_(4),bar_(7),bar_(7),bar_(4));
    fprintf('%10.5f  %10.5f %10.5f %10.5f   %4i   %-3i-%3i   %-3i-%3i   %3i\n', Y(1:10,:) );
    fprintf('\n');

    Z=[Z,Y];
end

% get irreducible shells (a negative p2i means the pair was flipped)
p2i = Z(end,:); [~,i2p,~]=unique(abs(p2i),'stable');

bar_ = @(x) repmat('-',[1,x]);
fprintf('%i irreducible shells\n',numel(i2p));
fprintf('%-10s   %-30s   %-4s   %-7s   %-7s   %-4s\n', 'd [cart]','bond [cart]','#','ic(i,j)','pc(m,n)','flp?'); 
fprintf('%-10s   %-30s   %-4s   %-7s   %-7s   %-4s\n', bar_(10),bar_(30),bar_(4),bar_(7),bar_(7),bar_(4));
fprintf('%10.5f  %10.5f %10.5f %10.5f   %4i   %-3i-%3i   %-3i-%3i   %3i\n', Z(1:10,i2p(rankcol_( Z([1],i2p)))) );
fprintf('\n');





% pairs = V(:,i2p(findrow_(A)));
% r = uc2ws(uc.bas*diff_(pairs),uc.bas);
% d = normc_(r);
% Z = [d;r;sum(A,2).';uc.u2i(pairs);uc.u2p(pairs);1:numel(d)];
% Z = Z(:, rankcol_([d;r]) );
% Z = Z(:, i2p(rankcol_( Z([8,1],i2p))) ) );


% V(:,i2p(A(5,:)))
% fprintf(' %10f %10f %10f %10f %3i %3i %3i %3i %3i %3i \n', Z)



% pairs = uc.u2p(V(:,i2p(A(6,:))))







