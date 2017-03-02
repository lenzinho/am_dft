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

% get irreducible shells
[ip,pp,p2i,i2p] = get_shells(pc,uc,cutoff); 

% save shells maps
pp.i2p = i2p; pp.p2i = p2i;
%%



import am_lib.*

% set sym digits
 digits(5);
 
% define function to change basis of rotational symmetry
sym_rebase_ = @(B,S) [[ matmul_(matmul_(B,S(1:3,1:3,:)),inv(B)), reshape(matmul_(B,S(1:3,4,:)),3,[],size(S,3))]; S(4,1:4,:)];

% get cart symmetries
Qc{1} = sym_rebase_(pc.bas,ip.Q{1}); Qc{2} = ip.Q{2};

% get form of force constants for irreducible prototypical bonds
for i = 1:ip.nshells
    % use stabilzer group to determine crystallographic symmetry relations; A*B*C' equals kron(C,A)*B(:)
    W = sum(kron_( Qc{1}(1:3,1:3,ip.s_ck(:,i)) , Qc{1}(1:3,1:3,ip.s_ck(:,i)))-eye(9),3);

    % enforce intrinsic symmetry (immaterial order of differentiation: c == c.')
    F = zeros(9,9); F(sub2ind([9,9],[1:9],[1,4,7,2,5,8,3,6,9])) = 1; W = W + F-eye(9);
    
    % get linearly-independent nullspace and normalize to first nonzero element
    W = null(W); W = frref_(W.').'; for j = 1:size(W,2); W(:,j) = W(:,j)./W(find(W(:,j),1),j); end; 

    % define parameters
    c = sym(sprintf('c%02i_%%d%%d',i),[3,3],'real'); c = c(findrow_(double(W).'));

    % get symmetry adapted force constants
    phi = reshape( sym(W,'d')*c(:),[3,3]); 

    % clean up zeros and ones in phi, sym-adapted FC 
    for j = 1:numel(phi); [C,T]=coeffs(phi(j)); C(C<1E-8)=0; C(abs(C-1)<1E-8)=1; phi(j) = sum(C.*T); end
    
    % print sym-adapted force constants
    fprintf('v = [%10f,%10f,%10f] \n', ip.v(:,i)); disp( phi )

    % save important stuff (sort W to be in line with c, matlabFunction sorts D variables)
    [sav.c{i},n] = sort(c(:).'); sav.W{i} = W(:,n); sav.shell{i} = ones(1,numel(c))*i; sav.phi(:,:,i) = phi;
end

% construct the dynamical matrix
nbands=3*max([pp.m,pp.n]); D=sym(zeros(nbands)); kvec=sym('k%d',[3,1],'real'); mass=sym('m%d',[1,max([pp.i,pp.j])]);
for p = 1:pp.nshells
    % take proto phi to rep phi
    phi = permute( Qc{1}(1:3,1:3,pp.Qi(p)) * sav.phi(:,:,p2i(p)) * Qc{1}(1:3,1:3,pp.Qi(p)).' , Qc{2}(:,pp.Qi(p)) );

    % get block of dynamical matrix 
    B = zeros(3,3); 
    for i = find(pp.g_ck(:,p).')
        if pp.Q{2}(1,i) > pp.Q{2}(2,i); sign_=+1; else; sign_=-1; end
        B = B + Qc{1}(1:3,1:3,i) * phi * Qc{1}(1:3,1:3,i)' .* exp( sym(2i*pi) * dot( sym(pp.Q{1}(1:3,1:3,i) * sign_ * pp.v(:,p),'d'), kvec) ); 
    end

    % update block matrix with masses
    B = B ./ sqrt(prod(mass([pp.i(p),pp.j(p)])));

    % augment dynamical matrix 
    mp = 3*(pp.m(p)-1)+[1:3]; 
    np = 3*(pp.n(p)-1)+[1:3];
    D( mp,np ) = D( mp,np ) + B; 
end

% simplify the dynamical matrix
D = simplify(D,'Steps',500);

% create bvk structure
bvk_ = @(ip,sav,D) struct('units','cart,frac-recp','bas',ip.bas,'recbas',ip.recbas,'natoms',numel(unique(ip.m)),'mass',ip.mass, ...
    'nshells',size(sav.W,2),'W',{sav.W},'shell',{sav.shell},'nbands',nbands,'D',matlabFunction(D));
bvk = bvk_(ip,sav,D);











