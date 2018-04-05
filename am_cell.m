classdef am_cell
    
    properties (Constant)
        cifdir    = '/Users/lenzinho/Developments/_materials/cif/';
    end
    
    properties
        % basic properties
        formula  = []; % 
        units    = []; % frac/cart
        bas      = []; % basis
        type     = []; % type of atom
        ntypes   = []; % number of types of atoms
        natoms   = []; % number of atoms
        nspecies = []; % number of atoms of each type
        species  = []; % species identifier
        tau      = []; % atomic coordinates
        tol      = []; % numerical tolernce
        % structural info
        mass_density     = []; % [g/cm3] mass density
        numb_density     = []; % [atoms/nm3] number density
        form_density     = []; % [f.u./nm3] formula unit density
        mole_weight      = []; % [amu] molecular weight
        % symmetry info
        bravais          = []; % bravais type
        holohodry        = []; % holohodry
        point_group      = []; % point group
        space_group      = []; % space group  
        % charge density info
        nchg = []; % [n1,n2,n3]
        chg  = [];
        % symmetry functions
    end
    
    methods (Static)

        function [uc] = define(varargin) % (bas,type,species,tau) or (type,nspecies,mass_density)
            

            switch nargin
                case {4}
                    % create dft cell
                    [bas,type,species,tau]= deal(varargin{1:4});
                    % bas      = []; % basis
                    % atom     = []; % atoms
                    % species  = []; % species identifier
                    % tau      = []; % atomic coordinates
                    % tol      = []; % numerical tolernce
                    uc          = am_cell;
                    uc.units    = 'frac';
                    uc.bas      = bas;
                    uc.ntypes   = numel(type);
                    uc.type     = type;
                    uc.natoms   = size(tau,2);
                    uc.nspecies = sum(unique(species)==species(:),1);
                    uc.species  = species;
                    uc.tau      = tau;
                case {3}
                    % create xrr cell
                    [type,nspecies,mass_density]= deal(varargin{1:3});
                    amu2gram  = 1/6.02214E23;
                    %
                    if ischar([type{:}]); type = am_atom.define(type); end
                    uc          = am_cell;
                    uc.units = 'frac';
                    uc.bas      = eye(3) * ( sum( [type(:).mass] .* nspecies ) / mass_density * amu2gram ).^(1/3) * 1E7;
                    uc.ntypes   = numel(type);
                    uc.type     = type;
                    uc.natoms   = sum(nspecies);
                    uc.nspecies = nspecies;
                    uc.species  = repelem(1:numel(type),uc.nspecies);
                    uc.tau      = [];
                otherwise
                    error('unknown input');
            end
            uc.tol          = 1E-6;
            uc.formula      = uc.get_formula();
            uc.mass_density = uc.get_mass_density();
            uc.mole_weight  = uc.get_molecular_weight();
            uc.numb_density = uc.get_atomic_density();
            uc.form_density = uc.get_formula_density();
        end
       
    end
    
    methods  % convert between cells

        function [pc,p2u,u2p]    = get_primitive_cell(uc)
            % [pc,p2u,u2p] = get_primitive_cell(uc)
            % NOTE: saves p2u entries which share a common closest
            % primitive lattice vector, not just the first primitive atoms
            % produced by the matrix A. When building shells, this property
            % is exploited.

            import am_lib.*

            % build permutation matrix for atoms related by translations
            T = get_symmetries(uc); nTs=size(T,2); PM=zeros(uc.natoms,nTs);
            for i = [1:nTs]; PM(:,i)=rankc_( [mod_(uc.tau(:,:,1)+T(1:3,i));uc.species] ); end

            % construct a sparse binary representation
            A=zeros(uc.natoms); A(sub2ind([1,1]*uc.natoms,repmat([1:uc.natoms].',nTs,1),PM(:)))=1; A=frref_(A); A=A(~all(A==0,2),:);

            % find the smallest primitive cell volume (the three smallest vectors which preserve periodic boundary conditions)
            inds=[0,0,0]; T_cart = uc.bas*T; fwd = rankc_(normc_(T_cart)); T_cart = T_cart(:,fwd); T = T(:,fwd);
            for j =           1:nTs; if sum(abs(T_cart(:,j)))                         >uc.tol; inds(1)=j; break; end; end
            for j = (inds(1)+1):nTs; if sum(abs(cross(T_cart(:,inds(1)),T_cart(:,j))))>uc.tol; inds(2)=j; break; end; end
            for j = (inds(2)+1):nTs; if abs(det(T_cart(:,[inds(1:2),j])+eye(3)*eps))  >uc.tol; inds(3)=j; break; end; end
            if any(inds==0); error('basis not found'); end

            % save basis
            B=T(:,inds); vol=abs(det(uc.bas*B));

            % refined the primitive cell vectors to get (if possible)
            % 1) angles between vectors to equal each other or 30, 60, 90, or 120 deg
            % 2) lengths of primitive vectors to be the same
            if true
                % augment lattice vectors
                grid_ =  [1, 0,1,1, 1,0,0, -1, 1, 1, 0, -1, 1,-1, 1, 0, 0, -1,-1, 1, -1, 0, 0, -1,-1, 0, -1; ...
                          1, 1,0,1, 0,1,0,  1,-1, 1, 0,  0, 0, 1,-1, 1,-1, -1, 1,-1,  0,-1, 0, -1, 0,-1, -1; ...
                          1, 1,1,0, 0,0,1,  1, 1,-1, 0,  1,-1, 0, 0,-1, 1,  1,-1,-1,  0, 0,-1,  0,-1,-1, -1];
                T = osum_(T,grid_,2);  T = uniquec_(T); T_cart = uc.bas*T;

                % (to speed things up) keep vectors whose norm are smaller than 2 times the largest norm in B the tentative basis
                ex_ = normc_(T_cart)<2*max(normc_(uc.bas*B)); nTs = sum(ex_); T = T(:,ex_); T_cart = T_cart(:,ex_);

                % find all combination of unit vectors that the smallest volume possible
                ijk = nchoosek_(nTs,3); m = size(ijk,2);
                v = zeros(1,m); for i = 1:m; v(i) = det(T_cart(:,ijk(:,i)));  if v(i)<0; ijk(:,i)=flipud(ijk(:,i)); v(i)=abs(v(i)); end; end
                ex_ = ~eq_(v,0); ex_(ex_) = eq_(v(ex_),vol); ijk = ijk(:,ex_); m = size(ijk,2);

                % find basis which produces the highest symmetry
                n = zeros(1,m);
                for i = 1:m
                    % get lattice with most symmetric metric tensor
                    M=T_cart(:,ijk(:,i));M=M.'*M.*[1 2 2; 2 1 2; 2 2 1]; n(1,i) = numel(uniquetol(M(:), uc.tol));
                    % for metric tensors with equal symmetries, get the one which has the
                    % most number of angles at 30,60,90,120 deg.
                    n(2,i) = numel(uniquetol([30,60,90,120,bas2abc(T_cart(:,ijk(:,i)))], uc.tol));
                end
                inds_=rankc_(n); B = T(:,ijk(:,inds_(1)));
            end

            % set identifiers (see NOTE: cannot simply use p2u = findrow_(A)!)
            p2u = member_(mod_(B*mod_(B\uc.tau(:,findrow_(A),1))),mod_(uc.tau(:,:,1))).'; u2p = ([1:size(A,1)]*A);

            % create structure
            pc = am_cell.define(uc.bas*B, uc.type, uc.species(p2u), mod_(B\uc.tau(:,p2u)) );
        end

        function [ic,i2p,p2i]    = get_irreducible_cell(pc)
            % idenitifes irreducible atoms

            import am_lib.*

            % get seitz matrices
            [~,~,S] = pc.get_symmetries();
            % define function to apply symmetries to position vectors
            seitz_apply_ = @(S,tau) mod_(reshape(matmul_(S(1:3,1:3,:),tau),3,[],size(S,3)) + S(1:3,4,:), pc.tol);
            % get permutation matrix and construct a sparse binary representation
            PM = member_(seitz_apply_(S,pc.tau),pc.tau, pc.tol*1.01); [~,i2p,p2i] = get_connectivity(PM);
            % create structure
            ic = am_cell.define(pc.bas, pc.type, pc.species(i2p), pc.tau(1:3,i2p) );
            % unassign these parameters which are meaningless
            ic.mass_density = [];
            ic.numb_density = [];
            ic.form_density = [];
            ic.mole_weight  = [];
        end

        function [xc,x2i,i2x]    = get_expanded_cell(ic, S)
            
            import am_lib.*
            
            % space groups
            nSs = size(S,3);
            % define function to apply symmetries to position vectors
            seitz_apply_ = @(S,tau) mod_(reshape(matmul_(S(1:3,1:3,:),tau),3,[],size(S,3)) + S(1:3,4,:));
            % expand atomic positions over symmetry basis
            X = [repmat(1:ic.natoms,1,nSs);reshape(seitz_apply_(S,ic.tau),3,[])]; [~,~,jc] = uniquec_(X); 
            % symmeterize atomic positions
            natoms=max(jc); tau = zeros(3,natoms); x2i = zeros(1,natoms);
            for i = 1:natoms
                tau(:,i) = mean(X(2:4,jc==i),2); 
                x2i(i) = X(1,find(jc==i,1));
            end
            % reorder to put irreducible atoms close together
            fwd = rankc_([x2i;tau]); x2i = x2i(fwd); tau = tau(:,fwd);
            % get mapping 
            [~,i2x,~]=unique(x2i,'stable');
            % create structure
            xc = am_cell.define(ic.bas, ic.type, ic.species(x2i), tau );
        end

        function [uc,u2p,p2u]    = get_supercell(pc, B)
            % [uc,u2p,p2u] = get_supercell(pc, B)
            % Example: to make a cell pc have the same shape as a reference cell ref, do this:
            % [~,ref] = load_cells('185_P63cm.poscar');
            % [~,pc]  = load_cells('194_P63mmc.poscar');
            % pc.bas = rotzd_(-30)*pc.bas; T = round(pc.bas\ref.bas);
            % sc = get_supercell(pc,T);
            %
            import am_lib.* am_dft.*
            
            % simple lattice transformations
            if ischar(B)
               switch B
                   case 'bcc'; B = (ones(3) -   eye(3))/2;
                   case 'fcc'; B =  ones(3) - 2*eye(3);
                   otherwise; error('Invalid lattice transformation.');
               end
            end
            
            % check if only diagonal entries are supplied
            if numel(B) == 3; B = diag(B); end
            
            % check size
            if size(B,1) ~= 3 || size(B,2) ~= 3; error('B must be a 3x3 matrix.'); end

            % check that B has integers
            if abs(mod_(det(B)))>am_lib.eps; error('determinant of B must be an integer'); end

            % generate primitive lattice vectors
            n=ceil(normc_(B.')); [Y{1:3}]=ndgrid(0:n(1),0:n(2),0:n(3)); nLs=prod(n+1); L=reshape(cat(3+1,Y{:})-1,[],3).';

            % expand atoms, coordinates supercell fractional, and reduce to primitive supercell
            % Here, uniquec_ should not use 'stable', since atoms in the primitive cell should go first
            % uniquec_ without stable has an inherent sort; alternatively, use uniquec_ (with 
            % 'stable' built-in and then sort after); this is important later when getting pairs and triplets
            X = uniquec_( [reshape(repmat([1:pc.natoms],nLs,1),1,[]); mod_(inv(B)*osum_(L,pc.tau,2))] );
            X = X(:,rankc_(X));

            % create mapping
            u2p = X(1,:); [~,p2u]=unique(u2p); p2u=p2u(:).';

            % create structure
            [uc] = am_cell.define(pc.bas*B, pc.type, pc.species(u2p), X(2:4,:));
        end
        
    end
    
    methods % cell properties
        
        function formula         = get_formula(uc)
            formula = ''; x = uc.nspecies./am_lib.gcd_(uc.nspecies);
            for j = 1:numel(x)
                % get formula
                formula = [formula,strtrim(uc.type(j).symb)];
                % remove 1's if present
                if x(j)~=1; formula=[formula,strtrim(num2str(x(j)))]; end
            end
        end

        function mass_density    = get_mass_density(uc)
            % mass_density = get_mass_density(uc) 
            % [g/cm3]
            amu2gram  = 1/6.02214E23; 
            mass_density = sum( [uc.type(uc.species).mass] ) * amu2gram / det(uc.bas * 1E-7);
        end

        function mole_weight     = get_molecular_weight(uc)
            mole_weight = sum([uc.type(uc.species).mass]); 
        end
        
        function num_density     = get_atomic_density(uc)
            % num_density = get_atomic_density(uc)
            % [atoms/nm3]
            num_density = uc.natoms / det(uc.bas);
        end
        
        function formula_density = get_formula_density(uc)
            % formula_density = get_formula_density(uc)
            % [f.u./nm3]
            formula_density = uc.get_atomic_density()/(uc.natoms/am_lib.gcd_(uc.nspecies));
        end

        function metric          = get_metric(uc)
            % metric = get_metric(bas)
            % Compute metric tensors from unit cell column basis vector.
            switch 2
                case 1
                    % explicit
                    metric = zeros(3,3);
                    for i = 1:3
                    for j = 1:i
                        metric(i,j)=dot(uc.bas(:,i),uc.bas(:,j)); metric(j,i)=metric(i,j);
                    end
                    end
                case 2
                    % matrix multiplication
                    metric = uc.bas.'*uc.bas;
            end
        end
        
        function [bas,T]         = get_niggli(uc)
            % [niggli_bas,T]  = get_niggli_(bas); 
            % niggli_bas == bas * T;
            %
            % Refs: I. Krivy and B. Gruber, Acta Crystallographica Section A 32, 297 (1976).
            %   	R. W. Grosse-Kunstleve, N. K. Sauter, and P. D. Adams, Acta
            %   	Crystallographica Section a Foundations of Crystallography 60, 1 (2004). 
            % 
            % Note: The exmaple in the first reference uses rounded numbers.
            %       bas=abc2bas([3 5.196 2 103.55 109.28 134.53]);
            %       M = get_metric(bas); 
            %       S = [M(1,1),M(2,2),M(3,3),M(2,3),M(1,3),M(1,2)];
            %       S = round(S); 
            %
            
            T = eye(3); [x,y,z,a,b,c,l,m,n,bas] = update_(uc.bas, eye(3), uc.tol);
           
            % begin reduction procedure
            for counter = 1:100
                if (a > b + uc.tol || (~ abs(a - b) > uc.tol && abs(x) >  abs(y) + uc.tol))
                    % Procedure A1
                    A = [0, -1, 0; -1, 0, 0; 0, 0, -1];
                    [x,y,z,a,b,c,l,m,n,bas] = update_(bas, A, uc.tol); T = T * A; %#ok<ASGLU>
                end
                if (b > c + uc.tol || (~ abs(b - c) > uc.tol && abs(y) >  abs(z) + uc.tol))
                    % Procedure A2
                    A = [-1, 0, 0; 0, 0, -1; 0, -1, 0];
                    [x,y,z,a,b,c,l,m,n,bas] = update_(bas, A, uc.tol); T = T * A;
                    continue
                end
                if l * m * n == 1
                    % Procedure A3
                    if l == -1; i = -1; else; i = 1; end
                    if m == -1; j = -1; else; j = 1; end
                    if n == -1; k = -1; else; k = 1; end
                    A = [i, 0, 0; 0, j, 0; 0, 0, k];
                    [x,y,z,a,b,c,l,m,n,bas] = update_(bas, A, uc.tol); T = T * A; %#ok<ASGLU>
                else
                    % Procedure A4
                    if l == 1; i = -1; else; i = 1; end
                    if m == 1; j = -1; else; j = 1; end
                    if n == 1; k = -1; else; k = 1; end
                    if i * j * k == -1
                        if l == 0; i = -1; end
                        if m == 0; j = -1; end
                        if n == 0; k = -1; end
                    end
                    A = [i, 0, 0; 0, j, 0; 0, 0, k];
                    [x,y,z,a,b,c,l,m,n,bas] = update_(bas, A, uc.tol); T = T * A; %#ok<ASGLU>
                end
                if ( abs(x) > b + uc.tol || (~ abs(b - x) > uc.tol && 2 * y < z - uc.tol) || (~ abs(b + x) > uc.tol && z < -uc.tol))
                    % Procedure A5
                    A = [1, 0, 0; 0, 1, - sign(x); 0, 0, 1];
                    [x,y,z,a,b,c,l,m,n,bas] = update_(bas, A, uc.tol); T = T * A;
                elseif ( abs(y) > a + uc.tol || (~ abs(a - y) > uc.tol && 2 * x < z - uc.tol) || (~ abs(a + y) > uc.tol && z < -uc.tol))
                    % Procedure A6
                    A = [1, 0, - sign(y); 0, 1, 0; 0, 0, 1];
                    [x,y,z,a,b,c,l,m,n,bas] = update_(bas, A, uc.tol); T = T * A;
                elseif ( abs(z) > a + uc.tol || (~ abs(a - z) > uc.tol && 2 * x < y - uc.tol) || (~ abs(a + z) > uc.tol && y < -uc.tol))
                    % Procedure A7
                    A = [1, - sign(z), 0; 0, 1, 0; 0, 0, 1];
                    [x,y,z,a,b,c,l,m,n,bas] = update_(bas, A, uc.tol); T = T * A;
                elseif (x + y + z + a + b < -uc.tol || (~ abs(x + y + z + a + b) > uc.tol && 2 * (a + y) + z > uc.tol))
                    % Procedure A8
                    A = [1, 0, 1; 0, 1, 1; 0, 0, 1];
                    [x,y,z,a,b,c,l,m,n,bas] = update_(bas, A, uc.tol); T = T * A;
                else
                    break;
                end
            end
            
            % check for completion
            if counter==100; error('Failed to reduce to Niggli cell. \n'); end
            
            function [x,y,z,a,b,c,l,m,n,bas] = update_(bas, A, tol)
                bas = bas*A; metric = bas.'*bas;
                a = metric(1,1); x = 2 * metric(2,3); l = 0;
                b = metric(2,2); y = 2 * metric(1,3); m = 0;
                c = metric(3,3); z = 2 * metric(1,2); n = 0;
                if x < -tol ; l = -1; elseif x > tol; l = 1; end
                if y < -tol ; m = -1; elseif y > tol; m = 1; end
                if z < -tol ; n = -1; elseif z > tol; n = 1; end
                % fprintf('%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n',a,b,c,x,y,z);
            end
        end

        function [T,H,S,R,W]     = get_symmetries(pc)
            % [T,H,S,R,W] = get_symmetries(pc, tol=am_dft.tiny)
            % T = all possible translations which restore the crystal to iteself
            % H = holohogries (all possible rotations which restore the bravais lattice onto iteself)
            % S = space group symmetries
            % R = point group symmetries

            import am_lib.*

            tol = pc.tol;
            
            % define function to check first two dimensions
            check_ = @(A) all(all(abs(A)<tol,1),2);

            % define function to sort atoms and species into a unique order (reference)
            X_ = @(tau,species) sortc_([mod_(tau);species]); X = X_(pc.tau(:,:,1),pc.species);

            % get all vectors connecting atom N to all other atoms
            % maybe expanding over supercells isn't necessary?
            % sc = pc;
            sc = pc.get_supercell([2,2,2]);
                N = 1; V =    mod_( sc.tau(:,sc.species==sc.species(N))-sc.tau(:,N) + 1/2, tol) - 1/2; 
                N = 1; V = [V,mod_( sc.tau(:,sc.species==sc.species(N))+sc.tau(:,N) + 1/2, tol) - 1/2]; 
                nVs=size(V,2);
            V = V*2;
            
            % find out which subset of vectors V preserve periodic boundary conditions
            ex_=false(1,nVs); 
            for j = 1:nVs; ex_(j) = check_( X_(pc.tau(1:3,:,1)-V(:,j),pc.species)-X ); end
            T=[V(:,ex_),eye(3)]; T=T(:,rankc_(normc_(T)));
            
            % condense and get unique 
                % maybe this isn't necessary?
                T = am_lib.mod_(T); T = am_lib.uniquec_(T); T = wdv_(T); T = [T,eye(3)];

            if nargout == 1; return; end

                % get arithmetic holodries (symmetries for which R'*g*R = g; g = bas'*bas)
                N=9; Q=[-1:1]; nQs=numel(Q);[Y{N:-1:1}]=ndgrid(1:nQs); L=reshape(Q(reshape(cat(N+1,Y{:}),[],N)).',3,3,[]);
                get_holodries_frac_ = @(M) L(:,:,check_(matmul_(matmul_(permute(L,[2,1,3]),M.'*M),L)-M.'*M));
                H = get_holodries_frac_(pc.bas); nHs = size(H,3);
                id = member_(flatten_(eye(3)),reshape(H,3^2,[])); H(:,:,[1,id])=H(:,:,[id,1]);

            if nargout == 2; return; end

                % get seitz operators which leave the atomic basis invariant
                S = zeros(4,4,nHs*nVs); S(4,4,:)=1; nSs=0;
                for i = 1:nHs; for j = 1:nVs
                    if check_( X_(H(:,:,i)*pc.tau+V(:,j),pc.species) - X ); nSs=nSs+1; S(1:3,1:4,nSs)=[ H(:,:,i), V(:,j) ]; end
                end; end; S = S(:,:,1:nSs);

                % condense and get unique 
                    % maybe this isn't necessary?
                    S(1:3,4,:) = am_lib.mod_(S(1:3,4,:));
                    S = reshape(am_lib.uniquec_(reshape(S,16,[])),4,4,[]); nSs = size(S,3);

                % set well defined values
                S=am_lib.wdv_(S);
                 
                % set identity first
                id = member_(flatten_(eye(4)),reshape(S,4^2,[])); S(:,:,[1,id])=S(:,:,[id,1]);

            if nargout == 3; return; end

                % get point symmetries
                R  = reshape(uniquec_( reshape(S(1:3,1:3,:),[9,nSs]) ),3,3,[]);

                % set identity first
                id = member_(flatten_(eye(3)),reshape(R,3^2,[])); R(:,:,[1,id])=R(:,:,[id,1]);
            
            if nargout == 4; return; end
            
                % get double group
                j=1/2; W = get_wigner(j,R,'spherical'); 

                % Add plus and minus in accordance with 
                % V. Heine, Group Theory in Quantum Mechanics (Elsevier, 2014), p 62, eq 8.24.
                W = cat(3,W,-W); 

                % remove numerical noise
                W = wdv_(W);

        end
        
    end
    
    methods % plotting

        function [h]             = plot(pc)

            import am_lib.* am_dft.*

            % initialize figure
            set(gcf,'color','w'); hold on;

            % plot atoms
            clist=am_lib.colormap_('spectral',max(pc.species));
            for i = 1:max(pc.species)
                ex_ = pc.species==i; radius = ones(1,sum(ex_)) * pc.type(i).r_ionic * 5000;
                h(i) = scatter3_(pc.bas*pc.tau(:,ex_), radius,'MarkerEdgeColor','k','MarkerFaceColor',clist(i,:));
            end
            
            % plot pc boundaries
            plothull_(pc.bas*[0,1,0,1,0,1,0,1;0,0,1,1,0,0,1,1;0,0,0,0,1,1,1,1]);
            
            hold off; daspect([1 1 1]); box on;

            % legend
            lh_ = legend(h,{pc.type(:).symb}); lh_.Box='off'; axis off;

        end

        function [F]             = plot_md(md, varargin)
            % n=[4;4;4]; kpt=[0;0;1/4]; amp=10; mode=6; nsteps=51;
            % [~,md] = get_displaced_cell(pc,bvk,n,kpt,amp,mode,nsteps);
            % clf; [F]=plot_md_cell(md,'view',[0;1;0]); movie(F,3); % write_poscar(md,'POSCAR_test')

            import am_lib.*

            if varargin{1}=='view'; v_xyz = varargin{2}; else; v_xyz = 3; end

            % initialize figure
            set(gcf,'color','w'); hold on;

            % plot cell boundaries
            plothull_(md.bas*[0,1,0,1,0,1,0,1;0,0,1,1,0,0,1,1;0,0,0,0,1,1,1,1]);
            daspect([1 1 1]); box on; axis tight; view(v_xyz); fixaxis=axis;

            % get coordinates in [cart]
            tau = matmul_(md.bas,md.tau);

            % plot paths for each atom
            hold on; plot3_(reshape(tau(:,:,:),3,[]),'.','markersize',5);

            % plot first point
            hold on; h = scatter3_(tau(:,:,1),50*sqrt(md.mass(md.species)),md.species(:),'filled','MarkerEdgeColor','k'); hold off;
            axis(fixaxis); drawnow; F(md.nsteps) = struct('cdata',[],'colormap',[]); F(1) = getframe;

            if md.nsteps>1
            for i = 2:md.nsteps
                [h.XData,h.YData,h.ZData] = deal(tau(1,:,i),tau(2,:,i),tau(3,:,i)); F(i) = getframe;
            end
            end
        end

        
    end
    
    % aux library

    methods (Static)

        % aux unit cells
        
        
        function uc           = load_material(material)
            switch material
            % toys (NOTE: CONTINUOUS SYMETRIES ARE NOT CHECKED FOR.)
            case '1D-chain';    bas = abc2bas([10,1],'tetra');                    tau = [[0;0;0]];                                                                                                                    symb = {'H'};                            sg_code = 1;
            case '1D-dimer';    bas = abc2bas([10,1],'tetra');                    tau = [[0;0;0],[0;0;0.5]];                                                                                                          symb = {'H','He'};                       sg_code = 1;
            case '2D-square';   bas = abc2bas([1,10],'tetra');                    tau = [[0;0;0],[1;1;0]/2,[1;0;0]/2,[0;1;0]/2,[0;0;1]/4];                                                                            symb = {'H','H','H','H','H'};            sg_code = 1;    % the last atom at [0,0,1/4] breaks z-mirror symmety so that only in-plane symmetries are considered
            case '2D-BN';       bas = abc2bas([1,1,10],'hex');                    tau = [[0;0;0],[2/3;1/3;0]];                                                                                                        symb = {'B','N'};                        sg_code = 187;
            case '2D-graphene'; bas = abc2bas([1,1,10],'hex');                    tau = [[2/3;1/3;0]];                                                                                                                symb = {'C'};                            sg_code = 191;
            case '3D-cube';     bas = abc2bas(1,'cubic');                         tau = [[0;0;0]];                                                                                                                    symb = {'H'};                            sg_code = 1;
            case '3D-NaCl';     bas = abc2bas(1,'cubic');                         tau = [[0;0;0], [1;1;1]/2];                                                                                                         symb = {'Na','Cl'};                      sg_code = 225;  
            % metals
            case 'fcc-Co';      bas = abc2bas(0.35441,'cubic');                   tau = [[0;0;0]];                                                                                                                    symb = {'Co'};                           sg_code = 225;  % ICSD  44989
            case 'hcp-Co';      bas = abc2bas([0.25054,0.40893],'hex');           tau = [[1/3;2/3;1/4]];                                                                                                              symb = {'Co'};                           sg_code = 194;  % ICSD  44990
            case 'CsCl-CoFe';   bas = abc2bas(0.28570,'cubic');                   tau = [[0;0;0],[1;1;1]/2];                                                                                                          symb = {'Co','Fe'};                      sg_code = 221;  % ICSD  56273
            case 'Cu';          bas = abc2bas(0.36151,'cubic');                   tau = [[0;0;0]];                                                                                                                    symb = {'Cu'};                           sg_code = 225;  % ICSD  43493
            case 'fcc-Fe';      bas = abc2bas(0.36468,'cubic');                   tau = [[0;0;0]];                                                                                                                    symb = {'Fe'};                           sg_code = 225;  % ICSD  44862
            case 'bcc-Fe';      bas = abc2bas(0.29315,'cubic');                   tau = [[0;0;0]];                                                                                                                    symb = {'Fe'};                           sg_code = 229;  % ICSD  44863
            % salts
            case 'NaCl';        bas = abc2bas(0.54533,'cubic');                   tau = [[0;0;0], [1;1;1]/2];                                                                                                         symb = {'Na','Cl'};                      sg_code = 225;  
            % oxides
            case 'Al2O3';       bas = abc2bas([0.47617,1.29990],'hex');           tau = [[0;0;0.3522],[0.6936;0;0.2500]];                                                                                             symb = {'Al','O'};                       sg_code = 167;  % ICSD  10425
            case 'gamma-Al2O3'; bas = abc2bas(0.79110,'cubic');                   tau = [[1;1;1]*0.37970,[5;5;5]/8,[0;0;0],[1;1;1]*0.15220];                                                                          symb = {'O','Al','Al','Al'};             sg_code = 227;  % ICSD  66558 - CIF has different origin choice
            case 'eta-Al2O3';   bas = abc2bas(0.79140,'cubic');                   tau = [[1;1;1]*0.37990,[5;5;5]/8,[1;0;0]*0.77390,[1;1;1]*0.19490];                                                                  symb = {'O','Al','Al','Al'};             sg_code = 227;  % ICSD  66559 - CIF has different origin choice
            case 'BiAlO3';      bas = abc2bas([0.537546,1.33933],'hex');          tau = [[0;0;0],[0;0;0.2222],[0.5326;0.0099;0.9581]];                                                                                symb = {'Bi','Al','O'};                  sg_code = 161;  % ICSD 171708
            case 'Bi2Al4O9';    bas = abc2bas([0.77134,0.81139,0.56914],'orth');  tau = [[0.1711;0.1677;0],[0.5;0;0.2645],[0.3545;0.3399;0.5],[0;0;0.5],[0.3718;0.2056;0.2503],[0.1364;0.412;0.5],[0.1421;0.4312;0]]; symb = {'Bi','Al','Al','O','O','O','O'}; sg_code = 55;   % ICSD  88775
            case 'SrTiO3';      bas = abc2bas(0.39010,'cubic');                   tau = [[0;0;0], [1;1;1]/2, [1;1;0]/2];                                                                                              symb = {'Sr','Ti','O'};                  sg_code = 221;  % ICSD   8087
            case 'PbTiO3';      bas = abc2bas([0.3902,0.4156],'tetra');           tau = [[0;0;0],[0.5;0.5;0.5377],[0.5;0.5;0.1118],[0;0.5;0.6174]];                                                                   symb = {'Pb','Ti','O','O'};              sg_code = 99;   % ICSD  61168
            % scandates
            case 'TbScO3';      bas = abc2bas([5.72920,7.91700,5.46540,90,90,90]);tau = [[0.0595;0.2500;0.0164],[0;0;.5],[0.4449;0.25;0.8798],[0.2993;0.4436;0.3082]];                                                symb = {'Tb','Sc','O','O'};              sg_code = 62;   
            % semiconductors
            case 'Si';          bas = abc2bas(0.54305,'cubic');                   tau = [[0;0;0]];                                                                                                                    symb = {'Si'};                           sg_code = 227;  % ICSD  51688
            case 'GaAs';        bas = abc2bas(0.5652,'cubic');                    tau = [[0;0;0], [1;1;1]/4];                                                                                                         symb = {'Ga','As'};                      sg_code = 216;  % ICSD 107946 
            % nitrides            
            case 'VN';  		bas = abc2bas(0.4134,'cubic'); 					  tau = [[0;0;0],[1;1;1]/2]; 																									      symb = {'V','N'};  					   sg_code = 225;
            case 'ScN'; 		bas = abc2bas(0.4501,'cubic'); 					  tau = [[0;0;0],[1;1;1]/2]; 																									      symb = {'Sc','N'}; 					   sg_code = 225;
            case 'TiN'; 		bas = abc2bas(0.4240,'cubic'); 					  tau = [[0;0;0],[1;1;1]/2]; 																									      symb = {'Ti','N'}; 					   sg_code = 225;
            case 'ZrN'; 		bas = abc2bas(0.4573,'cubic'); 					  tau = [[0;0;0],[1;1;1]/2]; 																									      symb = {'Zr','N'}; 					   sg_code = 225;
            case 'HfN'; 		bas = abc2bas(0.4524,'cubic'); 					  tau = [[0;0;0],[1;1;1]/2]; 																									      symb = {'Hf','N'}; 					   sg_code = 225;
            case 'CeN'; 		bas = abc2bas(0.5043,'cubic'); 					  tau = [[0;0;0],[1;1;1]/2]; 																									      symb = {'Ce','N'}; 					   sg_code = 225;
            case 'CrN'; 		bas = abc2bas(0.4162,'cubic'); 					  tau = [[0;0;0],[1;1;1]/2]; 																									      symb = {'Cr','N'}; 					   sg_code = 225;
            % otherwise
            otherwise; error('ERROR 1508123: unknown material');
            end
            uc = am_dft.create_cell(bas,tau,symb,sg_code);
        end
        

        function [uc,pc,ic,cc]= load_cell(flag, arg, tol)
            % [uc,pc,ic,cc]   = load_cell(fposcar)
            % fposcar can be a poscar or cif file

            import am_lib.* am_dft.*
            
            % stupid matlab requires these symbolic variables to be initialized
            x=[];y=[];z=[]; %#ok<NASGU>
            
            % set default numerical tolerance
            if nargin < 3; tol = am_dft.tiny; end
            
            % time
            fprintf(' ... getting cell'); tic
            
            % validate input
            validatestring(flag,{'poscar','cif','material'});
            switch flag
                case {'poscar','cif'};  if exist(arg,'file')~=2; fprintf('\n'); error('File does not exist: %s',arg); end
            end
            
            % convert bas from [Ang] to [nm] when loading from cif/poscar
            switch flag
                case 'poscar';   [uc]     = load_poscar(arg); % already converted within the function
                case 'cif';      [uc,str] = load_cif(arg);      uc.bas = uc.bas*0.1;
                case 'material'; [uc]     = load_material(arg);
                case 'create';   [uc]     = am_dft.create_cell(arg{:});
            end

            % get primitive cell
            [pc,p2u,u2p] = get_primitive_cell(uc, tol);

            % get irreducible cell
            [ic,i2p,p2i] = get_irreducible_cell(pc, tol);

            % get conventional cell
            [cc,c2p,p2c] = get_conventional_cell(pc,tol);

            % complete mapping
            u2i = p2i(u2p); i2u = p2u(i2p);
            u2c = p2c(u2p); c2u = p2u(c2p);
            c2i = p2i(c2p); i2c = p2c(i2p);

            % sace mapping to cells
            pc.p2i = p2i; pc.p2u = p2u; pc.p2c = p2c;
            pc.i2p = i2p; pc.u2p = u2p; pc.c2p = c2p; 

            ic.i2p = i2p; ic.i2u = i2u; ic.i2c = i2c;
            ic.p2i = p2i; ic.u2i = u2i; ic.u2i = c2i;

            uc.u2p = u2p; uc.u2i = u2i; uc.u2c = u2c;
            uc.p2u = p2u; uc.i2u = i2u; uc.c2u = c2u;

            cc.c2i = c2i; cc.c2u = c2u; cc.c2p = c2p;
            cc.i2c = i2c; cc.u2c = u2c; cc.p2c = p2c;

            % save bas2pc and tau2pc to convert [uc/cc-frac] to [pc-frac]
            uc.bas2pc = pc.bas/uc.bas; uc.tau2pc = pc.bas\uc.bas;
            cc.bas2pc = pc.bas/cc.bas; cc.tau2pc = pc.bas\cc.bas;

            % print basic symmetry info
            [~,H,~,R] = get_symmetries(pc, tol);
            bv_code = identify_bravais_lattice(pc.bas, tol);

            % holohodry should give same info as bravais lattice
            hg_code = identify_pointgroup(H); 
            pg_code = identify_pointgroup(R); 
            sg_code = identify_spacegroup(pg_code); % BETA

            % print relevant information
            verbose = true;
            if verbose
                fprintf(' (%.3f s) \n',toc);
                fprintf('     %-16s = %s\n','formula',get_cell_formula(uc));
                fprintf('     %-16s = %s\n','primitive',decode_bravais(bv_code));
                fprintf('     %-16s = %s\n','holohodry',decode_holohodry(hg_code));
                fprintf('     %-16s = %s\n','point group',decode_pg(pg_code));
                fprintf('     %-16s = %s\n','laue group',decode_laue(identify_laue(pg_code)));
                fprintf('     %-16s = %s\n','space group',strrep(cell2mat(join(decode_sg(sg_code),',')),' ',''));
                fprintf('     %-16s = %-8.3f [g/cm3] \n','mass density',get_cell_mass_density(uc));
                fprintf('     %-16s = %-8.3f [atoms/nm3]\n','number density',get_cell_atomic_density(uc));
                fprintf('     %-16s = %-8.3f [f.u./nm3]\n','formula density',get_cell_formula_density(uc));
                fprintf('     %-16s = %-8.3f [amu/f.u.]\n','molecular weight',get_cell_molecular_weight(uc));
                if contains(flag,'cif')
                fprintf('     %-16s = %s\n','create command',str);
                end
            end

            % sub functions

            function [uc, str]    = load_cif(fcif)

                % load file into memory
                str = am_lib.load_file_(fcif);

                % exclude blank lines
                ex_=strcmp(strtrim(str),''); str=strtrim(str(~ex_)); nlines=numel(str);

                % get name
                name_aliases = {...
                    '_chemical_name_mineral',...
                    '_chemical_name_systematic',...
                    '_chemical_formula_structural',...
                    '_chemical_formula_sum'};
                for alias = name_aliases
                    mineral = am_lib.extract_token_(str,alias{:});
                    if ~isempty(mineral), break; end
                end

                % get lattice basis
                abc = [...
                       am_lib.extract_token_(str,'_cell_length_a',true) ...
                       am_lib.extract_token_(str,'_cell_length_b',true) ...
                       am_lib.extract_token_(str,'_cell_length_c',true)];
                angles = [...
                       am_lib.extract_token_(str,'_cell_angle_alpha',true) ...
                       am_lib.extract_token_(str,'_cell_angle_beta' ,true) ...
                       am_lib.extract_token_(str,'_cell_angle_gamma',true)];
                if length(abc)<3;    abc    = [1 1 1]; end
                if length(angles)<3; angles = [90 90 90]; end
                bas = am_dft.abc2bas([abc,angles]);

                % get atomic positions and species
                symb_dataset={...
                    'He' ,'Li' ,'Be' ,'Ne' ,'Na' ,'Mg' ,'Al' ,'Si' ,'Cl' ,'Ar' ,'Ca' ,'Sc' ,'Ti' ,'Cr' ,'Mn' ,'Fe' ,'Co' ,'Ni' , ...
                    'Cu' ,'Zn' ,'Ga' ,'Ge' ,'As' ,'Se' ,'Br' ,'Kr' ,'Rb' ,'Sr' ,'Zr' ,'Nb' ,'Mo' ,'Tc' ,'Ru' ,'Rh' ,'Pd' ,'Ag' , ...
                    'Cd' ,'In' ,'Sn' ,'Sb' ,'Te' ,'Xe' ,'Cs' ,'Ba' ,'La' ,'Ce' ,'Pr' ,'Nd' ,'Pm' ,'Sm' ,'Eu' ,'Gd' ,'Tb' ,'Dy' , ...
                    'Ho' ,'Er' ,'Tm' ,'Yb' ,'Lu' ,'Hf' ,'Ta' ,'Re' ,'Os' ,'Ir' ,'Pt' ,'Au' ,'Hg' ,'Tl' ,'Pb' ,'Bi' ,'Po' ,'At' , ...
                    'Rn' ,'Fr' ,'Ra' ,'Ac' ,'Th' ,'Pa' ,'Np' ,'Pu' ,'Am' ,'Cm' ,'Bk' ,'Cf' ,'Es' ,'Fm' ,'Md' ,'No' ,'Lr' ,'Rf' , ...
                    'Db' ,'Sg' ,'Bh' ,'Hs' ,'Mt' ,'Ds' ,'Rg' ,'Uub','Uut','Uuq','Uup','Uuh','H'  ,'B'  ,'C'  ,'N'  ,'O'  ,'F'  , ...
                    'P'  ,'S'  ,'K'  ,'V'  ,'Y'  ,'I'  ,'W'  ,'U'  }; nelements = numel(symb_dataset); % NOTE: single letters last!
                position_aliases = {...
                    '_atom_site_attached_hydrogens', ... % icsd
                    '_atom_site_U_iso_or_equiv', ... % crystal maker
                    '_atom_site_fract_z'};
                for alias = position_aliases; ex_ = contains(str,alias{:}); if any(ex_)
                    j=find(ex_); i=1; tau=[];
                    switch alias{:}
                        case {'_atom_site_fract_z','_atom_site_U_iso_or_equiv'}
                            while and(~any(am_lib.strmatchi_(str{j+i},{'_','#','loop_'})),lt(i+j,nlines))
                                buffer = strsplit(str{j+i},' ');
                                tau(1,i) = sscanf(buffer{4},'%f');
                                tau(2,i) = sscanf(buffer{5},'%f');
                                tau(3,i) = sscanf(buffer{6},'%f');
                                for k = 1:nelements
                                    if contains(buffer{2},symb_dataset{k})
                                        species_symb{i}=symb_dataset{k}; break;
                                    end
                                end
                                i=i+1;
                            end
                        case '_atom_site_attached_hydrogens'
                            while lt(i+j-1,nlines) && ~any(am_lib.strmatchi_(str{j+i},{'_','#','loop_'}))
                                buffer = strsplit(str{j+i},' ');
                                tau(1,i) = sscanf(buffer{5},'%f');
                                tau(2,i) = sscanf(buffer{6},'%f');
                                tau(3,i) = sscanf(buffer{7},'%f');
                                for k = 1:nelements
                                    if contains(buffer{2},symb_dataset{k})
                                        species_symb{i}=symb_dataset{k}; break;
                                    end
                                end
                                i=i+1;
                            end
                    end
                    if ~isempty(tau); break; end
                end; end
                % replace symbolic species with a numerical value
                [~,u2i_species,species]=unique(species_symb,'stable'); species=species(:).'; symb = species_symb(u2i_species);
                % make sure there is only one atom at each site
                [~,i]=am_lib.uniquec_(tau); tau=tau(:,i); species=species(i); 
                % [~,i,j]=unique(species); symb=symb(i); species=j(:).';

                % get space group symmetries
                position_aliases = {...
                    '_symmetry_equiv_pos_as_xyz', ... % crystal maker
                    '_symmetry_Int_Tables_number', ... % icsd
                    };
                for alias = position_aliases; ex_ = contains(str,alias{:}); if any(ex_)
                    switch alias{:}
                        case '_symmetry_equiv_pos_as_xyz'
                            j=find(ex_); i=1; syms x y z;
                            while ~isempty(strtrim(str{j+i})) && ~any(am_lib.strmatchi_(str{j+i},{'_','#','loop_'})) && lt(i+j,nlines)
                                buffer = strsplit(str{j+i},'''');
                                buffer = strsplit(buffer{2},','); 
                                tmp=coeffs(eval(buffer{1}),'all'); T(1) = double(tmp(end));
                                tmp=coeffs(eval(buffer{2}),'all'); T(2) = double(tmp(end));
                                tmp=coeffs(eval(buffer{3}),'all'); T(3) = double(tmp(end));
                                S(1:3,1:3,i) = double(equationsToMatrix(...
                                    [eval(buffer{1}),eval(buffer{2}),eval(buffer{3})] ));
                                S(1:3,4:4,i) = am_lib.mod_(T(:));
                                S(4:4,1:3,i) = 0;
                                S(4:4,4:4,i) = 1;
                                i=i+1;
                            end
                            sg_id = [];
                        case '_symmetry_Int_Tables_number'
                            % true == generate from memory
                            sg_id = am_lib.extract_token_(str,'_symmetry_Int_Tables_number',true); S_generated = am_dft.generate_sg(sg_id);
                            if any(any( ~am_lib.eq_(am_lib.sortc_(reshape(S,16,[])),am_lib.sortc_(reshape(S_generated,16,[]))) ))
                                fprintf('\n'); warning('Generate space symmetries do not match those in cif file!');
                            end
                    end
                end; end
                nSs = size(S,3);
                
                % save this cell?
                str_bas = sprintf('am_dft.abc2bas([%g,%g,%g,%g,%g,%g])',[abc/10,angles]);
                str_tau = sprintf('[%g;%g;%g],',tau); str_tau = ['[',str_tau(1:(end-1)),']'];
                str_spc = sprintf('''%s'',',symb{species}); str_spc = ['{',str_spc(1:(end-1)),'}'];
                str = sprintf('create_cell(%s,%s,%s,%i)', str_bas, str_tau, str_spc, sg_id);

                % generate all positions based on symmetry operations
                seitz_apply_ = @(S,tau) am_lib.mod_(reshape(am_lib.matmul_(S(1:3,1:3,:),tau),3,[],size(S,3)) + S(1:3,4,:));
                tau = reshape(am_lib.mod_(seitz_apply_(S,tau)),3,[]); species=repmat(species,1,nSs);
                [tau, j] = am_lib.uniquec_(tau); species=species(j);
                [species,i]=sort(species); tau=tau(:,i);

                % define primitive cell creation function and make structure
                uc_ = @(bas,symb,species,tau) struct('units','frac','bas',bas, ...
                    'symb',{symb},'mass',am_dft.get_atomic_mass(symb),'nspecies',sum(unique(species).'==species,2).', ...
                    'natoms',numel(species),'tau',tau,'species',species);
                uc = uc_(bas,symb,species,tau);
                
                uc_gen = [];
                eval(['uc_gen = am_dft.',str]);
                if uc_gen.natoms ~= uc.natoms || ...
                   am_lib.any_(~am_lib.eq_(am_lib.sortc_(uc_gen.tau),am_lib.sortc_(uc.tau)))
                    str = 'Mismatch';
                end
            end
            
        end
        
        function [uc]         = load_poscar(fposcar,flag) % can also be used to load the charge density 
            if nargin < 2; flag=''; end
            fid=fopen(fposcar,'r'); if fid==-1; error('Unable to open %s',fposcar); end
                header=fgetl(fid);                 % read header (often times contains atomic symbols)
                latpar=sscanf(fgetl(fid),'%f');    % read lattice parameter
                a1=sscanf(fgetl(fid),'%f %f %f');  % first basis vector
                a2=sscanf(fgetl(fid),'%f %f %f');  % second basis vector
                a3=sscanf(fgetl(fid),'%f %f %f');  % third basis vector
                bas=latpar*[a1,a2,a3];             % construct the basis and convert to nm (column vectors)
                bas=bas/10;                        % convert basis from angstroms to nm
                buffer=fgetl(fid);                 % check vasp format
                if ~isempty(sscanf(buffer,'%f'))
                    % vasp 4 format (symbols are missing, jumps straight to species)
                    symb=regexp(header, '([^ \s][^\s]*)', 'match');
                    nspecies=sscanf(buffer,repmat('%f' ,1,length(symb)))';
                else
                    % vasp 5 format (symbols are present)
                    symb=regexp(buffer, '([^ \s][^\s]*)', 'match');
                    nspecies=sscanf(fgetl(fid),repmat('%f' ,1,length(symb)))';
                end
                % create atom type
                type = am_atom(); 
                for i = 1:numel(symb)
                    type(i) = am_atom.define(symb{i},1); 
                end
                % continue parsing
                coordtype=lower(strtrim(fgetl(fid)));
                l=0;
                for i=1:length(nspecies)
                    for j=1:nspecies(i); l=l+1;
                        tau(:,l)=sscanf(fgetl(fid),'%f %f %f');
                        species(l)=i;
                    end
                end
                if ~strcmp(coordtype(1),'d'); tau=bas\tau*latpar; end
                % make cell structure
                uc = am_cell.define(bas,type,species,tau);
                % load charge density
                if contains(flag,'charge')
                    fgetl(fid);
                    n = sscanf(fgetl(fid),'%i');
                    t = textscan(fid,'%f');
                    uc.nchg = n(:).';
                    uc.chg  = reshape(t{:},uc.n);
                end
            fclose(fid);
        end

        function [cc,c2p,p2c] = get_conventional_cell(pc, tol, algo)
            % [vc,v2p,p2v] = get_conventional_cell(pc,tol)
            
            import am_dft.* am_lib.*

            if nargin < 2; tol = am_dft.tiny; end
            if nargin < 3; algo = 2; end 
            % 2 seems more robust, at least for VN, Fe2O3, and Al2O3
            
            % get point symmetries
            [~,~,~,R] = get_symmetries(pc);

            % this boils down on how to get the centering matrix
            switch algo
                case 0
%                     % test
%                     fposcar='~/Developments/_materials/poscar/VN_Fm-3m.poscar';
%                     % fposcar='./Al2O3_167_R-3c_corundum.poscar';
%                     % fposcar='Fe3O4_Fd-3m_magnetite_spinel.poscar';
%                     [uc,pc,ic,cc]=load_cells(fposcar);
%                     [~,~,~,g] = get_symmetries(pc);
                    
                case 1
                    % find the transformation matrix mapping the pc point group to the one in standard setting
                    % find a transformation matrix which takes g to the standard setting h
                    [~,M] = find_pointgroup_transformation(R,generate_pg(identify_pointgroup(R)),1);
                    % get integer inverse matrices
                    C = inv(M); C = C/det(C);
                    
                case 2
                    % Sets the non-collinear triplet of highest rotation axes as edges of the
                    % cell. The algorithm is based on:
                    %
                    % Refs: V. L. Himes and A. D. Mighell, Acta Crystallographica Section a
                    %         Foundations of Crystallography 43, 375 (1987). 
                    % 
                    
                    % exclude identity and pure inversion
                    nRs = size(R,3); inds = [1:nRs]; 
                    [~,tr,dt] = identify_point_symmetries(R);
                    filter_ = @(tr,dt,inds,ex_) deal(tr(ex_),dt(ex_),inds(ex_));
                    [tr,dt,inds] = filter_(tr,dt,inds, tr ~=-3 );
                    [tr,dt,inds] = filter_(tr,dt,inds, tr ~= 3 );

                    % get the rotation axis of each proper rotation
                    nTs = numel(inds); T = zeros(3,nTs);
                    for i = 1:numel(inds)
                        % get null space (axis of rotation)
                        T(:,i) = R_axis_(R(:,:,inds(i)));
                        % convert vector to all integers
                        % T(:,i) = round_(T(:,i));
                    end

                    % sort rotation axes:
                    %   1) proper rotation axes first
                    %   1) higher order rotation axes first
                    %   2) vectors in positive octant first
                    cc_bas = pc.bas*T; 
                    fwd = rankc_([-sign(dt);-tr;-max(sign(cc_bas));sign(cc_bas)]);
                    filter_ = @(tr,dt,T,vc_bas,inds,ex_) deal(tr(ex_),dt(ex_),T(:,ex_),vc_bas(:,ex_),inds(ex_));
                    [~,~,T,~,~] = filter_(tr,dt,T,cc_bas,inds, fwd );

                    % append identity
                    T = [T,eye(3)]; cc_bas = pc.bas*T; nTs=nTs+3;
                    
                    % find the three highest symmetry non-collinear vectors and centering
                    % transformation matrix from primitive to conventional cell basis  
                    % i.e. C = pc.bas -> vc_bas; absolute highest symmetry is z, then y, then x, hence [k,j,i] and not [i,j,k]
                    for i = (  1):nTs; X = cc_bas(:,i);                    if any(~eq_(X, 0, tol)); break; end; end
                    for j = (i+1):nTs; X = cross(cc_bas(:,i),cc_bas(:,j)); if any(~eq_(X, 0, tol)); break; end; end
                    for k = (j+1):nTs; X = det(cc_bas(:,[i,j,k]));         if any(~eq_(X, 0, tol)); break; end; end
                    if any([k,j,i]==0); error('Conventional basis not found!'); else; C = T(:,[k,j,i]); end
                    
            end
            
            % construct conventional cell
            [cc,c2p,p2c] = get_supercell(pc,C);
        end

        function [uc]         = create_cell(bas,tau,symb,sg_code)
            % create_cell(bas,tau,symb,sg_code) creates cell based on wyckoff positions and standard crystallographic setting
            %
            % Example input for BiFeO3 P63cm (hypothetical polymorph):
            %     % define basis, atomic positions, species, and elements
            %     a = 6.200; c = 12.076;
            %     bas = abc2bas([a,a,c,90,90,120]);
            %     tau=[0.00000 0.00000 0.48021; 0.33333 0.66667 0.01946; ...
            %       0.29931 0.00000 0.15855; 0.63457 0.00000 0.34439; ...
            %       0.33440 0.00000 0.00109; 0.00000 0.00000 0.27855; ...
            %       0.33333 0.66667 0.23206].';
            %     symb = {'O','O','O','O','Fe','Bi','Bi'}; sg_code=185;
            %     % create cell
            %     [uc] = create_cell(bas,tau,symb,sg_code); [bragg] = get_bragg(uc);
            %     plot_bragg(bragg); tabulate_bragg(bragg,10);
            %     % save poscar
            %     write_poscar(uc,'BiFeO3_P63cm_dft.poscar');
            %

            import am_dft.* am_lib.*

            % identify atomic types and assign a number label
            [symb,~,species]=unique(symb,'stable'); species=species(:).';

            % get space symmetries in conventional setting
            S = am_dft.generate_sg(sg_code);

            % define function to sort atoms and species into a unique order (reference)
            X_ = @(tau,species,c_i2u) sortc_([mod_(tau);species;c_i2u]); 
            
            % apply symmetry operations to all atoms            
            natoms = size(tau,2); X = X_(tau,species,[1:natoms]); X = apply_symmetry(S,X); X(1:3,:)=mod_(X(1:3,:));

            % get unique species
            [~,ind] = am_lib.uniquec_(X(1:4,:)); tau = X(1:3,ind); species = X(4,ind); c_i2u = X(5,ind);

            % sort by species
            [~,fwd] = sort(c_i2u); tau = tau(:,fwd); species = species(fwd); c_i2u = c_i2u(fwd);

            % define irreducible cell creation function and make structure
            uc_ = @(bas,tau,symb,species) struct('units','frac','bas',bas,...
                'symb',{symb},'mass',am_dft.get_atomic_mass(am_dft.get_atomic_number(symb)),...
                'nspecies',sum(unique(species).'==species,2).', ...
                'natoms',size(tau,2),'tau',tau,'species',species);
            uc = uc_(bas,tau,symb,species);
        end
        

        
        function [uc]         = stack_cell(uc)
            % put one cell on top of the other along the z direction
            if numel(uc)==1; return; end
            
            % stack along z = [0;0;1] direction
            bas = [uc{1}.bas(:,1:2), uc{1}.bas(:,3) + uc{2}.bas(:,3)];
            tau = [uc{1}.tau,[0;0;1]+uc{1}.bas\[uc{1}.bas(:,1:2),uc{2}.bas(:,3)]*uc{2}.tau];
            tau = bas\uc{1}.bas*tau;

            % get symbols and species
            symb = [uc{1}.symb(uc{1}.species),uc{2}.symb(uc{2}.species)];
            [~,u2i,species]=unique(symb,'stable'); symb=symb(u2i); species=species(:).';

            % define primitive cell creation function and make structure
            uc_ = @(bas,symb,species,tau) struct('units','frac','bas',bas, ...
                'symb',{symb},'mass',am_dft.get_atomic_mass(symb),'nspecies',sum(unique(species).'==species,2).', ...
                'natoms',numel(species),'tau',tau,'species',species);
            uc = uc_(bas,symb,species,tau);
        end
        
        function [uc,inds]    = match_cell(uc,uc_ref)

            import am_lib.* am_dft.*

            res = Inf;
            for i = 1:uc_ref.natoms
                % place first atom ontop of reference atom i
                shift = -uc.tau(:,1,1)+uc_ref.tau(:,i,1);
                tau = mod_(uc.tau+shift);
                % find closest matching atom
                [m,fwd] = min( reshape(normc_( ...
                    mod_(osum_(-tau(:,:,1),uc_ref.tau(:,:,1),2)+.5)-.5 ),...
                                                        uc.natoms,uc.natoms) );
                % compute residual
                m = norm(m);
                % compare
                if m < res
                    res = m; inds = fwd; s_shift = shift;
                end
            end
            uc.tau = mod_(uc.tau+s_shift);

            % shuffle uc atoms so that order matchs uc_ref
            for f = {'tau','species'}
                if isfield(uc,f); uc.(f{:}) = uc.(f{:})(:,inds,:); end
            end
            for f = {'u2p','u2i'}
                if isfield(uc,f); uc.(f{:}) = uc.(f{:})(inds); end
            end
            for f = {'p2u','i2u'}
                if isfield(uc,f); uc.(f{:}) = inds(uc.(f{:})); end
            end
        end
        
    end
    
    

    
end
