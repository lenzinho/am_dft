% get irreducible and primitive shells; this function has goals: 
%   1) save [v] orbit reps
%   2) save [s_ck] the space symmetries which stabilize the bond
%   3) save [g_ck] the space symmetries which generate the orbit
%   4) save [Qi]   the space symmetry   which maps the prototypical
%                  irreducible pair onto the rep [v]
%
% NOTE: requires p2u to point to atoms who's uc coordinates share a common
% closest primitive lattice point. set this in get_primitive_cell
% 
% NOTE: application of space symmerties requires the permutation matrices
% P1s and P2s to not have any zero values. all atoms must be mapped onto
% another atom.
%
% NOTE: The primitive shells as defined here involve only two primitive
% cell atoms for each given orbit. Although some space symmetries may map
% atoms onto other primitive cell atoms, those orbits are decomposed into
% two unique primitive shells which are related to the prototypical
% irreducible shell via Qi.



clear;clc


import am_lib.*


tiny = am_lib.tiny;

cutoff=5;
fname='POSCAR';
% fname='POSCAR.BMg2';
% fname='outfile.supercell';
flags='';

[uc,pc,ic] = get_cells(fname,flags);

[ip,pp,p2i,i2p] = get_shells(pc,uc,cutoff);

%%

% define function to apply symmetries to position vectors
sym_apply_ = @(S,tau) reshape(matmul_(S(1:3,1:3,:),tau),3,[],size(S,3)) + S(1:3,4,:);

% define function to change basis of seitz symmetry
sym_rebase_ = @(B,S) [[ matmul_(matmul_(B,S(1:3,1:3,:)),inv(B)), reshape(matmul_(B,S(1:3,4,:)),3,[],size(S,3))]; S(4,1:4,:)];

% define function to get the closest primitive lattice vector in supercell basis
G_ = @(uc2pc,tau) tau - reshape(matmul_(uc2pc,mod_(matmul_(inv(uc2pc),tau))),size(tau));

% define function to get distance [cart] from positions [frac]
d_ = @(r) normc_(uc2ws(uc.bas*r,uc.bas));

% define function to get bond vector from unit cell atom pair indicies
diff_ = @(V) uc.tau(:,V(1,:))-uc.tau(:,V(2,:));

% define pair saving function
pair_ = @(uc,Q,Z) struct('units','frac','bas',uc.bas,'recbas',uc.recbas, ...
    'symb',{{uc.symb{:}}},'species',uc.species,'mass',uc.mass,...
    'nshells',size(Z,2),'norbits',Z(5,:),'v',Z(2:4,:),'i',Z(6,:),'j',Z(7,:),'m',Z(8,:),'n',Z(9,:), ...
    's_ck',Z([10+[1:size(Q,3)]],:),'g_ck',Z([10+size(Q,3)+[1:size(Q,3)]],:),'Qi',Z(end,:),'nQs',size(Q,3),'Q',Q);

% ----------------------------------------------------------------------- %


% readjust cutoff based on unitcell
cutoff = min([normc_(uc.bas)/2,cutoff]);

% get conversion from primitive to supercell basis
uc2pc = pc.bas/uc.bas;

% get space symmetries in [supercell basis]
[S,~] = get_symmetries(pc); nSs = size(S,3); S = sym_rebase_(uc2pc,S);

% designate forth row of S to denote permutations, save permutation+space symmetry as Q
M = [[1,2,0,0];[2,1,0,0]].';
Q = repmat(S,1,1,size(M,2));
Q(4,:,:) = permute(repelem(M,1,nSs),[3,1,2]);
nQs=size(S,3)*size(M,2);

% determine primitive atoms shells which are connected
P1s=[];P2s=[]; npcs=numel(uc.p2u); npairs=zeros(1,npcs);
for i = 1:npcs
    % identify primitive atom in unit cell
    Delta = uc.tau(:,uc.p2u(i)); 

    % compute cutoff distances, exclude atoms above cutoff, count pairs involving the i-th primitive atom
    d = normc_(uc2ws(uc.bas*(uc.tau-Delta),uc.bas)); ex_ = [d<cutoff]; npairs(i) = sum(ex_);

    % compute action of space symmetries on pair positions
    tau{1} = repmat(sym_apply_(S,Delta),[1,npairs(i),1]); tau{2} = sym_apply_(S,uc.tau(:,ex_));
    
    % at least one atom must be in the primitive cell
    G=G_(uc2pc,tau{1}); P11 = member_(mod_(tau{1}-G),uc.tau); P12 = member_(mod_(tau{2}-G),uc.tau);
    G=G_(uc2pc,tau{2}); P21 = member_(mod_(tau{1}-G),uc.tau); P22 = member_(mod_(tau{2}-G),uc.tau); % these bonds are flipped below

    % save results (P1s are always in the primitive cell) (MECHANISM A)
    P1s = [P1s;[P11,P22]]; P2s = [P2s;[P12,P21]];
end

% create a unique pair label
[V,~,p2i]=unique([P1s(:),P2s(:)],'rows'); V=V.';

% get permutation representation (entries are unique pair indicies)
PM = reshape(p2i,size(P1s));

% get map connecting pairs (sort A based on distances)
A = get_connectivity(PM); 

% sort rows of A based on pair distances [cart]
A = A( rankcol_( d_(diff_(V(:,findrow_(A)))) ),:);

% save irrep pair information
nirreps = size(A,1); PMi = findrow_(A); PMs = PM(PMi,:); 
w    = sum(diff(sort(PMs.').',1,2)~=0,2).'+1; % weights
xy   = V(:,PMi);   % uc indicies
mn   = uc.u2p(xy); % pc indicies
ij   = uc.u2i(xy); % ic indicies
v    = uc2pc\diff_(xy);  % bond vector [primitive frac]
d    = normc_(v);  % bond length [frac]
s_ck = [PMs==PMi].';
g_ck = zeros(size(s_ck));
Qi   = repelem(find(prod(s_ck,2)),1,size(s_ck,2));

% save pair
Z    = [d;v;w;ij;mn;[1:nirreps];s_ck;g_ck;Qi];
ip   = pair_(uc,Q,Z); 

% print irreducible pairs (NOTE: bonds and distances are converted to cart for printing)
Z(2:4,:) = uc2ws(uc.bas*uc2pc*Z(2:4,:),uc.bas); Z(1,:) = normc_(Z(2:4,:));
bar_ = @(x) repmat('-',[1,x]);
fprintf('primitive atom %i at [%.5f,%.5f,%.5f] has %i shells\n', i,pc.bas*pc.tau(:,i),sum(ex_));
fprintf('%-10s    %-30s   %-4s   %-7s   %-7s   %-4s   %-4s\n', 'd [cart]','bond [cart]','#','ic(i,j)','pc(m,n)','irr.','Qi'); 
fprintf('%-10s    %-30s   %-4s   %-7s   %-7s   %-4s   %-4s\n', bar_(10),bar_(30),bar_(4),bar_(7),bar_(7),bar_(4),bar_(4));
fprintf('%10.5f  %10.5f %10.5f %10.5f   %4i   %-3i-%3i   %-3i-%3i   %4i   %4i\n', Z([1:10,end],:) );
fprintf('\n');

Z=[];
for i = 1:pc.natoms
    for j = 1:nirreps
        % get pairs centered on irreducible atom i that belong to irreducible shell j
        ex_ = [V(1,PMs(j,:))==pc.p2u(i)];
        % record the first occurance of each value in PMs(j,:)
        ey_ = uniquemask_(PMs(j,:));
        % if such a shell exists, analyze it
        if any(ex_)
            % get pairs which involve only two different primitive atoms
            for k = 1:pc.natoms
                ez_ = [uc.u2p(V(2,PMs(j,:)))==k];
                g_ck = and( and(ex_ ,ez_), ey_ );
                if any(g_ck)

                    w   = sum(g_ck);      % number of pairs in orbit
                    Qi  = find(g_ck,1);   % record operation which takes PMi(j) -> xy
                    PMr = PMs(j,Qi);      % rep index
                    xy  = V(:,PMr);       % uc indicies
                    mn  = uc.u2p(xy).';   % pc indicies
                    ij  = uc.u2i(xy).';   % ic indicies
                    v   = uc2pc\diff_(xy);% bond vector [primitive frac]
                    d   = normc_(v);      % bond length [primitive frac]
                    s_ck=[PMs(j,:)==PMr]; % record stabilizers

%                     AA = mod_(sym_apply_( Q(:,:,Qi) , diff_(  V(Q(4,1:2,Qi),PMi(j)) ) ));
%                     BB = mod_(diff_(xy));
%                     max(abs(AA(:)-BB(:)))

                    % [ d(1), r(2,3,4), w(5), ij(6,7), mn(8,9), irres(10), s_ck, g_ck ]
                    Z = [Z,[d(:);v(:);w(:);ij(:);mn(:);j;s_ck(:);g_ck(:);Qi]];
                end
            end
        end
    end
end

% save primitive pairs
pp = pair_(uc,Q,Z);

% print primitive pairs
Z(2:4,:) = uc2ws(uc.bas*uc2pc*Z(2:4,:),uc.bas); Z(1,:) = normc_(Z(2:4,:));
nprims = size(Z,2); fprintf('%i primitive shells\n',nprims);
for i = 1:pc.natoms
    bar_ = @(x) repmat('-',[1,x]); ex_ = Z(8,:)==i;
    fprintf('atom %i at [%.5f,%.5f,%.5f] has %i shells\n', i,pc.bas*pc.tau(:,i),sum(ex_));
    fprintf('%-10s    %-30s   %-4s   %-7s   %-7s   %-4s   %-4s\n', 'd [cart]','bond [cart]','#','ic(i,j)','pc(m,n)','irr.','Qi'); 
    fprintf('%-10s    %-30s   %-4s   %-7s   %-7s   %-4s   %-4s\n', bar_(10),bar_(30),bar_(4),bar_(7),bar_(7),bar_(4),bar_(4));
    fprintf('%10.5f  %10.5f %10.5f %10.5f   %4i   %-3i-%3i   %-3i-%3i   %4i   %4i\n', Z([1:10,end],ex_) );
    fprintf('\n');
end










