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
% fname='POSCAR.BMg2';
% fname='outfile.supercell';
flags='';

[uc,pc,~] = get_cells(fname,flags);




% define function to apply symmetries to position vectors
sym_apply_ = @(S,tau) reshape(matmul_(S(1:3,1:3,:),tau),3,[],size(S,3)) + S(1:3,4,:);

% define function to change basis of seitz symmetry
sym_rebase_ = @(B,S) [[ matmul_(matmul_(B,S(1:3,1:3,:)),inv(B)), reshape(matmul_(B,S(1:3,4,:)),3,[],size(S,3))]; S(4,1:4,:)];

% define function to get the clostest primitive lattice vector in supercell basis
G_ = @(uc2pc,tau) tau - reshape(matmul_(uc2pc,mod_(matmul_(inv(uc2pc),tau))),size(tau));



% readjust cutoff based on unitcell
cutoff = min([normc_(uc.bas)/2,cutoff]);

% get conversion from primitive to supercell basis
uc2pc = pc.bas/uc.bas;

% define function to get bond vector from unit cell atom pair indicies
diff_ = @(V) uc.tau(:,V(1,:))-uc.tau(:,V(2,:));

% get point symmetries in [supercell basis]
[S,~] = get_symmetries(pc); nSs = size(S,3); S = sym_rebase_(uc2pc,S);

% determine primitive atoms shells which are connected
P1s=[];P2s=[]; npcs=numel(uc.p2u); npairs=zeros(1,npcs);
for i = 1:npcs
    % identify primitive atom in unit cell
    Delta = uc.tau(:,uc.p2u(i)); 

    % compute cutoff distances, exclude atoms above cutoff, count pairs involving the i-th primitive atom
    d=normc_(uc2ws(uc.bas*(uc.tau-Delta),uc.bas)); ex_=[d<cutoff]; npairs(i)=sum(ex_);

    % first atom in pair (must be in primitive cell); record G which shifts tau back to itself 
    tau = sym_apply_(S,Delta); G=G_(uc2pc,tau); P1 = member_(mod_(tau-G),uc.tau);

    % second atom in pair (SEE NOTEs)
    tau = sym_apply_(S,uc.tau(:,ex_));          P2 = member_(mod_(tau-G),uc.tau);

    % save results
    P1s = [P1s;repelem(P1,npairs(i),1)]; P2s = [P2s;P2]; %#ok<AGROW>
end

% create a unique pair label
V=[P1s(:),P2s(:)].'; [~,i2p,p2i]=unique(sort(V).','rows','stable'); V=V(:,i2p);

% get permutation representation (entries are unique pair indicies)
PM = reshape(p2i,[sum(npairs),nSs]);

% get map connecting pairs (sort A based on distances)
A = get_binary_rep(PM); 

% sort rows of A based on pair distances
d_ = @(r) normc_(uc2ws(uc.bas*r,uc.bas)); 
A = A( rankcol_( d_(diff_(V(:,findrow_(A)))) ),:);

% find the symmetry corresponding to the identity
identity = member_(reshape(eye(4),16,1),reshape(S,16,[]));

% bundle
bundle_ = @(ex_,A,B) deal(A(:,ex_),B(ex_));

Z=[];
for i = 1:pc.natoms
    pc_id = pc.p2u(i);
    
    Y=[];
    for irres = 1:size(A,1)
        % get all connected shells which involve irreducible atom i
        [U,Uinds]=bundle_(A(irres,:),V,[1:size(V,2)]); ex_=any(U==pc_id,1);  
        
        % if such a shell exists, analyze it
        if any(ex_)
            % trim
            U=U(:,ex_);
            
            % choose the rep as the bond involving the greatest number of primitive cell atoms 
            [~,rep] = max(sum(ismember(U,pc.p2u),1)); 
            
            % see if it is flipped
            if U(1,rep)==pc_id; isflipped=+1; else; isflipped=-1; end
            
            % sort indicies to put pc_id at "center"
            U(U==pc_id)=0; U=sort(U); U(U==0)=pc_id;

            % record basic things
            xy  = U(:,rep);     % uc indicies
            mn  = uc.u2p(xy).'; % pc indicies
            ij  = uc.u2i(xy).'; % ic indicies
            v   = diff_(xy); v = uc2ws(uc.bas*v,uc.bas); % bond vector
            d   = normc_(v);    % bond length
            w   = sum(ex_);     % number of pairs in orbit

            % record stabilizer and generators
            row = find(PM(:,identity)==Uinds(rep),1);
            [~,b,c] = unique( PM(row,:) ) ;
            s_ck = [c(:).'==c(identity)];
            g_ck = false(nSs,1); g_ck(b(:)) = true;

            % rep = Uinds(find(ex_,1));
            % row = find(PM(:,identity)==rep,1);
            % x = mod_(diff_(V(:,PM(row,:))));
            % y = mod_(reshape(sym_apply_(S,diff_(V(:,rep))),3,[]));
            % max(abs(x(:)-y(:))) % should be equal to zero (and is)
            
            % [ d(1), r(2,3,4), w(5), ij(6,7), mn(8,9), irres(10), s_ck, g_ck ]
            Y = [Y,[d(:);+v(:);w(:);ij([1,2]);mn([1,2]);+irres*isflipped;s_ck(:);g_ck(:)]];
        end
    end

    bar_ = @(x) repmat('-',[1,x]);
    fprintf('primitive atom %i at [%.5f,%.5f,%.5f] has %i shells\n', i,pc.bas*pc.tau(:,i),size(Y,2));
    fprintf('%-10s    %-30s   %-4s   %-7s   %-7s   %-4s\n', 'd [cart]','bond [cart]','#','ic(i,j)','pc(m,n)','flp?'); 
    fprintf('%-10s    %-30s   %-4s   %-7s   %-7s   %-4s\n', bar_(10),bar_(30),bar_(4),bar_(7),bar_(7),bar_(4));
    fprintf('%10.5f  %10.5f %10.5f %10.5f   %4i   %-3i-%3i   %-3i-%3i   %3i\n', Y(1:10,:) );
    fprintf('\n');

    Z=[Z,Y];
end

% get irreducible shells (a negative p2i means the pair was flipped)
p2i = Z(10,:); [~,i2p,~]=unique(abs(p2i),'stable');

% save structures
v_ = @(Z) Z(2:4,:); w_ = @(Z) Z(5,:); i_ = @(Z) Z(6,:); j_ = @(Z) Z(7,:); m_ = @(Z) Z(8,:); n_ = @(Z) Z(9,:); 
s_ck_ = @(Z) Z(10+[1:nSs],:); g_ck_ = @(Z) Z(10+nSs+[1:nSs],:);

pair_ = @(uc,S,Z) struct('units','frac','bas',uc.bas,'recbas',uc.recbas, ...
    'symb',{{uc.symb{:}}},'species',uc.species,'mass',uc.mass,...
    'nshells',size(Z,2),'norbits',w_(Z),'v',v_(Z),'i',i_(Z),'j',j_(Z),'m',m_(Z),'n',n_(Z), ...
    's_ck',s_ck_(Z),'g_ck',g_ck_(Z),'nRs',size(S,3),'S',S);

pp = pair_(uc,S,Z);
ip = pair_(uc,S,Z(:,i2p(rankcol_( Z([1],i2p)))));

% print irreducible pairs
bar_ = @(x) repmat('-',[1,x]);
fprintf('%i irreducible shells\n',numel(i2p));
fprintf('%-10s    %-30s   %-4s   %-7s   %-7s   %-4s\n', 'd [cart]','bond [cart]','#','ic(i,j)','pc(m,n)','flp?'); 
fprintf('%-10s    %-30s   %-4s   %-7s   %-7s   %-4s\n', bar_(10),bar_(30),bar_(4),bar_(7),bar_(7),bar_(4));
fprintf('%10.5f  %10.5f %10.5f %10.5f   %4i   %-3i-%3i   %-3i-%3i   %3i\n', Z(1:10,i2p(rankcol_( Z([1],i2p)))) );
fprintf('\n');








