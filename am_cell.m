classdef am_cell < dynamicprops % required for xrr simulation
    
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
        % properties
        nchg             = []; % [n1,n2,n3]
        chg              = []; % charge density
        % elastic constants
        cijkl            = []; % elastic constants
        % xrr simulation
        % [dynamic] thickness 
        % [dynamic] filling
        % [dynamic] roughness
    end
    
    methods (Static)

        function [uc]            = define(varargin) % (bas,type,species,tau) or (type,nspecies,mass_density)
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

        function [pc,p2u,u2p]    = get_primitive(uc)
            % [pc,p2u,u2p] = get_primitive(uc)

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
                    n(2,i) = numel(uniquetol([30,60,90,120,am_cell.bas2abc(T_cart(:,ijk(:,i)))], uc.tol));
                end
                inds_=rankc_(n); B = T(:,ijk(:,inds_(1)));
            end

            % set identifiers (NOTE: cannot simply use p2u = findrow_(A)!)
            p2u = member_(mod_(B*mod_(B\uc.tau(:,findrow_(A),1))),mod_(uc.tau(:,:,1))).'; u2p = ([1:size(A,1)]*A);

            % create structure
            pc = am_cell.define(uc.bas*B, uc.type, uc.species(p2u), mod_(B\uc.tau(:,p2u)) );  pc.tol = uc.tol;
        end

        function [ic,i2p,p2i]    = get_irreducible(pc)
            % idenitifes irreducible atoms

            import am_lib.*

            % get seitz matrices
            [~,~,S] = pc.get_symmetries();
            % define function to apply symmetries to position vectors
            seitz_apply_ = @(S,tau) mod_(reshape(matmul_(S(1:3,1:3,:),tau),3,[],size(S,3)) + S(1:3,4,:), pc.tol);
            % get permutation matrix and construct a sparse binary representation
            PM = member_(seitz_apply_(S,pc.tau),pc.tau, pc.tol*1.01); [~,i2p,p2i] = get_connectivity(PM);
            % create structure
            ic = am_cell.define(pc.bas, pc.type, pc.species(i2p), pc.tau(1:3,i2p) ); ic.tol = pc.tol;
            % unassign these parameters which are meaningless
            ic.mass_density = [];
            ic.numb_density = [];
            ic.form_density = [];
            ic.mole_weight  = [];
        end

        function [xc,x2i,i2x]    = get_orbitcellcell(ic, S)
            
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
            xc = am_cell.define(ic.bas, ic.type, ic.species(x2i), tau ); xc.tol = ic.tol;
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
            uc = am_cell.define(pc.bas*B, pc.type, pc.species(u2p), X(2:4,:)); uc.tol = pc.tol;
        end
                
        function [dc,idc]        = get_displaced(pc,bvk,n,kpt,amp,mode,nsteps)
            % Note: can set mode=[] for interactive selection
            % n=[4;4;4]; kpt=[0;0;1/4]; amp=10; mode=6; nsteps=51;
            % [~,md] = get_displaced(pc,bvk,n,kpt,amp,mode,nsteps);
            % clf; [F]=plot_md_cell(md,'view',[0;1;0]); movie(F,3); % write_poscar(md,'POSCAR_test')
            %
            % Q: What effect does the kpt have on vibrational wave?
            % A: Gamma-center and zone boundary points are standing waves.
            %    Positive and negative k-points correspond to waves
            %    traveling in opposite directions.
            %
            % Q: Are the units for dt correct?
            % A: Yes! Check with:
            %
            %       plot(flatten_(diff(u(:,:,:),1,3)./dt), ...
            %               flatten_(v(:,:,2:end)+v(:,:,1:(end-1)))/2,'.');
            %
            % Q: How do traveling waves work?
            % A: Check it out using the code below:
            %
            %         clc
            %         N = 50;
            %         x = [0:(N-1)]/(N-1);
            %         t = [0:(N-1)]/(N-1);
            %
            %         f = 1;
            %         k = 2;
            %         E1 = exp(2i*pi*(+k*x.'+f*t));
            %         E2 = exp(2i*pi*(-k*x.'+f*t));
            %         for i = 1:N
            %         plot(x,real(E1(:,i)),...
            %              x,real(E2(:,i)),...
            %              x,real(E1(:,i)+E2(:,i))/2);
            %         ylim([-1 1]); line([0 1],[0 0]); drawnow;
            %         end
            %
            % Q: What is the reason for producing two phonon excitations
            %    when fourier transforming real displacements?
            % A: Real part of the phonon corresponds to the displacement,
            %    however for a phonon with wavevector k and -k, the real
            %    parts are identical:
            %
            %     +k = [0;0;+1/2]     -k = [0;0;-1/2]
            %    0.0000 + 0.0000i    0.0000 - 0.0000i
            %    0.0000 + 0.0000i    0.0000 - 0.0000i
            %    0.0000 + 0.0000i    0.0000 - 0.0000i
            %   -0.5525 - 0.4015i   -0.5525 + 0.4015i
            %   -0.1480 - 0.1076i   -0.1480 + 0.1076i
            %   -0.0000 - 0.0000i   -0.0000 + 0.0000i
            %    0.5525 + 0.4015i    0.5525 - 0.4015i
            %    0.1480 + 0.1076i    0.1480 - 0.1076i
            %   -0.0000 + 0.0000i   -0.0000 + 0.0000i
            %

            import am_lib.* am_dft.*

            % get a supercell commensurate with the kpoint
            uc = get_supercell(pc,diag(n)); [~,pp] = get_pairs(pc,uc,bvk.cutoff);

            % get phonon energies and eigenvectors at k-point
            bz = get_fbz(pc,[1,1,1]); bz.k = kpt; bz = get_bvk_dispersion(bvk,bz);

            % select a mode
            fprintf('Energies [meV]\n');fprintf('%5.2f \n',bz.hw*1E3);
            if isempty(mode)
                mode = input('Select mode: ');
            else
                fprintf('Mode selected: %i \n',mode);
            end

            % build force constant matrices
            for m = 1:pp.pc_natoms
                phi{m} = zeros(3,3*pp.npairs(m));
            for j = 1:pp.npairs(m)
                % get indicies and irrep force constants
                i = pp.i{m}(j); iq = pp.iq{m}(j); iphi = reshape(bvk.W{i}*bvk.fc{i}(:),3,3);
                % rotate force constants from irrep to orbit
                phi{m}(1:3,[1:3]+3*(j-1)) = pp.Q{1}(1:3,1:3,iq) * permute(iphi,pp.Q{2}(:,iq)) * pp.Q{1}(1:3,1:3,iq).';
            end
            end

            % initialize all arrays [cart]
            u = single(zeros(3,uc.natoms,nsteps));
            v = single(zeros(3,uc.natoms,nsteps));
            f = single(zeros(3,uc.natoms,nsteps));

            % convert phonon energies back to wierd units
            hw = real(bz.hw)./am_lib.units_eV; hw(hw(:)<1E-8)=1E-8;

            % get normal transformations from uc eigenvector
            % (Wallace p 113 10.40, Wallace p 115 eq 10.48)
            % I think there is something wrong with the way the normal mode
            % is being generated. Two modes are getting excited at +/-
            % about gamma. It seems like it's a standing wave rather than a
            % propagating wave ...
            U   = expand_bvk_eigenvectors(bvk,uc,bz);
            q2u = U   ./ repelem(sqrt(uc.mass(uc.species)).',3,1);
            u2q = U'  .* repelem(sqrt(uc.mass(uc.species)).',3,1).';
            v2p = U.' .* repelem(sqrt(uc.mass(uc.species)).',3,1).';

            % build q_sk wave vector
            t = 2*pi*[0:(nsteps-1)]/(nsteps-1)/hw(mode);
            q_sk(bvk.nbranches,nsteps) = 0;
            q_sk(mode,:) = amp * exp(1i*t(:)*hw(mode));

            % displace according to the phonon mode
            shp_ = @(A) reshape(A,3,uc.natoms);
            for i = 1:nsteps
                % get displacement and velocities in [cart : Ang & Ang/fs]
                % (Wallace p 113 10.40)
                u(:,:,i) = real(shp_( q2u * ( q_sk(:,i)            ) ));
                v(:,:,i) = real(shp_( q2u * ( q_sk(:,i).*hw(:).*1i ) ));

                % evaluate forces F on each atom : fc [eV/Ang^2] * u [Ang]
                for m = 1:pp.pc_natoms
                    f(:,pp.c{m},i) = ...
                            - phi{m}*reshape(u(:,pp.o{m},i),size(pp.o{m}).*[3,1]);
                end
            end

            % get normal modes (Wallace p 115 eq 10.49)
            q_sk = (u2q*reshape(u,[],nsteps));
            p_sk = (v2p*reshape(v,[],nsteps));

            % get energies (Wallace p 115 eq 10.53)
            PE(:,1) = -dot(reshape(u,[],nsteps),reshape(f,[],nsteps),1)/2;
            KE(:,1) = reshape(sum(uc.mass(uc.species).*dot(v,v,1),2),1,[])/2;
            PE(:,2) = dot(abs(q_sk).*hw(:),abs(q_sk).*hw(:),1)/2;
            KE(:,2) = dot(abs(p_sk)       ,abs(p_sk)       ,1)/2;

            plot(1:nsteps,KE(:,1),'-',1:nsteps,PE(:,1),'-',1:nsteps,KE(:,1)+PE(:,1),'-',...
                 1:nsteps,KE(:,2),'.',1:nsteps,PE(:,2),'.',1:nsteps,KE(:,2)+PE(:,2),'.');
            legend('KE','PE','KE+PE','nKE','nPE','nKE+nPE'); axis tight;

            % set time step in [fs]
            dt = (t(2)-t(1));

            % convert [cart] to [frac] and u to tau
            tau = matmul_(inv(uc.bas),u)+uc.tau;
            v   = matmul_(inv(uc.bas),v) / sqrt(103.6382);
            f   = matmul_(inv(uc.bas),f);

            % create displaced structure
            dc_ = @(uc,f,tau,v,dt) struct('units','frac',...
                'bas',uc.bas,'tau2pc',uc.tau2pc,'bas2pc',uc.bas2pc,...
                'symb',{{uc.symb{:}}},'mass',uc.mass,'nspecies',uc.nspecies, ...
                'natoms',uc.natoms,'force',f,'tau',tau,'vel',v,'species',uc.species, ...
                'dt',dt,'nsteps',size(tau,3),'u2p',uc.u2p,'p2u',uc.p2u,'u2i',uc.u2i,'i2u',uc.i2u);
            dc = dc_(uc,f,tau,v,dt);

            % get "primitive" irreducible displaced cell
            if nargout==2
                % get primitive cell basis commensurate with the displacement
                [cc,c2d,d2c] = get_primitive(dc); cc.bas = double(cc.bas);

                % reduce dc size
                idc_ = @(dc,cc,pc,c2d,d2c) struct('units',dc.units, ...
                    'bas',cc.bas,'bas2pc',pc.bas/cc.bas,'tau2pc',pc.bas\cc.bas,...
                    'symb',{{dc.symb{:}}},'mass',dc.mass, ...
                    'nspecies',sum(dc.species(c2d).'==[1:max(dc.species(c2d))],1),'natoms',numel(c2d), ...
                    'force',mod_(matmul_(cc.bas\dc.bas,dc.force(:,c2d,:))), ...
                    'tau',  mod_(matmul_(cc.bas\dc.bas,dc.tau(:,c2d,:))),...
                    'vel',       matmul_(cc.bas\dc.bas,dc.vel(:,c2d,:)), ...
                    'species',dc.species(c2d),'dt',dc.dt,'nsteps',dc.nsteps, ...
                    'u2p',dc.u2p(c2d),'p2u',d2c(dc.p2u),'u2i',dc.u2i(c2d),'i2u',d2c(dc.i2u));
                idc = idc_(dc,cc,pc,c2d,d2c);
            end
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
        
        function [eta]           = get_xray_refractive_index(uc,hv,k)
            %
            % n = get_xray_refractive_index(uc,hv,k)
            % 
            % n              [unitless]     n = 1 + delta + i beta
            % hv             [eV]           photon energy
            % k              [1/Ang]        wavevector
            %
            % Eq. 3 in B. L. Henke, E. M. Gullikson, and J. C. Davis,
            % Atomic Data and Nuclear Data Tables 54, 181 (1993). 
            %

            import am_lib.*
            
            % atomic number density per species [atoms/nm^3/species]
            atomic_density_per_species = uc.numb_density .* (uc.nspecies ./ uc.natoms );
            
            % get the effective form factor averaged over species [nhvs,nks]
            nhvs = numel(hv); nks = numel(k); f = zeros(nhvs,nks);
            for i = 1:numel(uc.type)
                f = f + uc.type(i).get_xray_form_factor(hv, k) .* atomic_density_per_species(i);
            end
            
            % refractive index Wormington eq 4.1
            eta = 1 - am_lib.r_0 ./ (2*pi) .* get_photon_wavelength(hv).^2 .* f;
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
            % [~,md] = get_displaced(pc,bvk,n,kpt,amp,mode,nsteps);
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

        function [h]             = plot_charge(uc,v)

            % set isolevels
            if nargin<2
                % fwd_ = @(x) log(x); rev_ = @(x) exp(x);
                fwd_ = @(x) x; rev_ = @(x) x;
                v = rev_(linspace(fwd_(min(uc.chg(:))),fwd_(max(uc.chg(:))),20));
            end

            % generate coordinates [cart]
            mp_ = @(x) [0:(x-1)]./x; [r{1:3}] = ndgrid(mp_(uc.nchg(1)),mp_(uc.nchg(2)),mp_(uc.nchg(3)));
            r = reshape(uc.bas*reshape(permute(cat(4,r{:}),[4,1,2,3]),3,prod(uc.nchg)),[3,uc.nchg]);

            % plot isolevels
            figure(1); clf; set(gcf,'color','w'); hold on;
            sq_ = @(x) reshape(x,uc.nchg); nvs = numel(v); clist=am_lib.colormap_('virdis',nvs);
            for i = 1:nvs
                h=patch(isosurface(sq_(r(1,:)),sq_(r(2,:)),sq_(r(3,:)),uc.chg,v(i)));
                h.FaceColor=clist(i,:); h.FaceAlpha=0.2; h.EdgeColor='none';
            end
            box on; axis tight; daspect([1 1 1]);
            
            % sq_ = @(x) reshape(x,uc.nchg);
            % hold on;
            %     h=slice(permute(sq_(r(1,:)),[2,1,3]),permute(sq_(r(2,:)),[2,1,3]),sq_(r(3,:)),uc.chg,uc.bas(1,1)/2,uc.bas(2,2)/2,uc.bas(3,3)/2); 
            %     for i = 1:numel(h); h(i).EdgeColor='none'; h(i).FaceColor='interp'; end
            % hold off;
        end
        
        function [h]             = plot_sound(uc) % plots sound velocities
            % [pc]=am_cell.load_poscar('poscar');
            % pc.cijkl=am_dft.get_vasp('outcar:elastic');
            % pc.plot_sound();

            density = uc.mass_density; % [g/cm3]
            C = uc.cijkl; % [GPa]
            
            % generate sphere
            n = 100;
            switch 'A'
                case {'M','manually'}
                    [phi,chi] = ndgrid(2*pi*[0:2*n]/(2*n), pi*[0:n]/n); [x,y,z]=sph2cart(chi,phi,1);
                case {'A','auto'}
                    [x,y,z]=sphere(n); [chi,phi,~] = cart2sph(x,y,z); 
                otherwise; error('unknown ssphere generation method');
            end

            % flatten and get sizes
            [n,m] = size(phi); R = [x(:),y(:),z(:)].';

            % loop over Rs
            nRs = size(R,2); v=zeros(3,nRs); 
            for iR = 1:nRs
                % build Christoffel's matrix
                G = zeros(3,3);
                for i = 1:3;for j = 1:3; for k = 1:3; for l = 1:3
                    G(i,k) = G(i,k) + C(i,j,k,l).*R(j,iR).*R(l,iR);
                end; end; end; end
                % solve Christoffel's equation
                v(:,iR) = sort(sqrt(eig(G)/density)); % [km/s]
            end

            % plot the data
            figure(1); clf; set(gcf,'color','w'); vspan = am_lib.minmax_(v(:));
            for ipol = 1:3
                % select velocity
                w = reshape(v(ipol,:),[n,m]);
                % convert velocity to radis
                [x,y,z]=sph2cart(chi,phi,w);
                % plot
                h(ipol) = subplot(1,3,ipol);surf(x,y,z,w,'edgecolor','none');
                daspect([1 1 1]); axis tight; axis off; axis([-1 1 -1 1 -1 1]*vspan(2)); caxis(vspan);    
                % draw axes
                line([0,-3*vspan(2)], [0,0], [0,0], 'LineWidth', 3, 'Color', 'k');
                line([0,0], [0,-3*vspan(2)], [0,0], 'LineWidth', 3, 'Color', 'k');
                line([0,0], [0,0],  [0,3*vspan(2)], 'LineWidth', 3, 'Color', 'k');
                grid on;
            end
            hl = linkprop(h,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
            setappdata(gcf, 'StoreTheLink', hl);

            am_lib.colormap_('virdis',100);  
        end
    end
    
    methods (Static) % [make this access protected] % structure-related symmetries
        
        function taup         = apply_symmetry(S,tau)
            %
            % tau(1:3,:) = atomic coordinates
            % tau(4:end) [optional] = anything associated with tau that you want to transfer 
            
            import am_lib.*

            % get symmetry type from shape
            [s(1),s(2),s(3)] = size(S);
            if     s(1)==3 && s(2)==3 % point symmetry
                apply_ = @(R,tau) cat(1, ...
                    reshape(matmul_(R(1:3,1:3,:),tau),3,[],size(R,3))                  , repmat(tau(4:end,:),1,1,size(S,3)) );
            elseif s(1)==4 && s(2)==4 && all_(eq_(S(4,1:4,:), [0,0,0,1])) % seitz symmetry
                apply_ = @(S,tau) cat(1, ...
                    reshape(matmul_(S(1:3,1:3,:),tau(1:3,:)),3,[],size(S,3))+S(1:3,4,:), repmat(tau(4:end,:),1,1,size(S,3)) );
            elseif s(1)==1 && s(2)==2 && iscell(S) % composit symmetry Q
                apply_ = @(Q,tau) apply_composite_symmetry_(Q,tau);
                [s(1),s(2),s(3)] = size(S{1});
            end
                
            taup = apply_(S,tau);
            taup = reshape(taup,[size(tau),s(3)]);

            function taup = apply_composite_symmetry_(Q,tau)
                % apply symmetry
                for i = 1:size(tau,2)
                    taup(:,i,:,:) = am_dft.apply_symmetry(Q{1},tau(:,i,:,:)); 
                end
                % apply permutation
                for iq = 1:size(Q{1},3)
                    taup(:,Q{2}(:,iq),:,iq) = taup(:,:,:,iq); 
                end
            end
        end

        function [ps_id,t,d,o]= identify_point_symmetries(R,tol)
            import am_lib.* am_dft.*
            
            % number of symmetries
            nsyms=size(R,3); 
            
            % check if double group is input and convert it to SO(3)
            if size(R,1) == 2
                Rp   = SU2_to_SO3(R); Rp = wdv_(Rp); 
                D    = get_wigner(1/2,Rp,'spherical');
                dg_p = ~permute(all(all(eq_(R,+D),1),2),[1,3,2]);
                dg_m = ~permute(all(all(eq_(R,-D),1),2),[1,3,2]);
                dg   = dg_m - dg_p;
                R    = Rp;
            else
                dg = ones(1,nsyms);
            end
            
            % classify
            ps_id = zeros(1,nsyms); t = zeros(1,nsyms); d = zeros(1,nsyms);
            for i = 1:nsyms
                % get trace and determinant (fractional)
                t(i) = trace(R(1:3,1:3,i)); d(i) = det(R(1:3,1:3,i));
                if     (eq_(t(i),+3) && eq_(d(i),+1)); ps_id(i) = 1;  % 'e'
                elseif (eq_(t(i),-1) && eq_(d(i),+1)); ps_id(i) = 2;  % 'c_2'
                elseif (eq_(t(i),+0) && eq_(d(i),+1)); ps_id(i) = 3;  % 'c_3'
                elseif (eq_(t(i),+1) && eq_(d(i),+1)); ps_id(i) = 4;  % 'c_4'
                elseif (eq_(t(i),+2) && eq_(d(i),+1)); ps_id(i) = 5;  % 'c_6'
                elseif (eq_(t(i),-3) && eq_(d(i),-1)); ps_id(i) = 6;  % 'i'
                elseif (eq_(t(i),+1) && eq_(d(i),-1)); ps_id(i) = 7;  % 's_2'
                elseif (eq_(t(i),+0) && eq_(d(i),-1)); ps_id(i) = 8;  % 's_6'
                elseif (eq_(t(i),-1) && eq_(d(i),-1)); ps_id(i) = 9;  % 's_4'
                elseif (eq_(t(i),-2) && eq_(d(i),-1)); ps_id(i) = 10; % 's_3'
                else;                                  ps_id(i) = 0;  % unknown
                end
            end
            
            % make double group values negative
            ps_id = ps_id .* dg;
                
            % finally, get the order
            o = get_order(R,tol);
            
            function order = get_order(S, tol)

                % define comparison function
                check_ = @(x) any(~am_lib.eq_(x(:),0,tol));

                % get sizes
                s = size(S); if numel(s)<3; s(3) = 1; end; nSs = s(3);

                if s(1) == 4 && s(2) == 4       % space symmetry
                    order = ones(1,nSs); 
                    for j = 1:nSs
                        X=S(1:4,1:4,j);
                        while check_(X - eye(4))
                            order(j) = order(j) + 1;
                            X = X*S(:,:,j); X(1:3,4)=mod_(X(1:3,4));
                        end
                    end
                else                            % point symmetry 
                    order = ones(1,nSs); n = size(S,1);
                    for j = 1:nSs
                        X=S(:,:,j);
                        while check_(X - eye(n))
                            order(j) = order(j) + 1;
                            X = X*S(:,:,j);
                        end
                    end
                end
            end
        end

        function bv_code      = identify_bravais(bas, tol, algo)
            % bv_code = identify_bravais_lattice(bas, tol, algo)
            %
            % Identifies the type of crystal system given a basis. 
            %
            % bas     : unit cell column basis vectors 
            % tol     : numerical tolerence
            % algo (1): de Graef p 86
            %      (2): Tinkham  p 61
            % 
            % Antonio Mei Sep 20 2017
            
            import am_lib.*
            
            % set default numerical tolerance
            if nargin < 2; tol = am_dft.tiny; end
            if nargin < 3; algo = 2; end
            
            bv_code = 0; 
            % select the algorithm
            switch algo
                case 1
                    % de Graef p 86
                    M  = bas.'*bas;
                    ii = [M(1,1),M(2,2),M(3,3)];
                    ij = [M(1,2),M(1,3),M(2,3)];
                    t  = numel(uniquetol_(ii, tol));
                    o  = numel(uniquetol_(ij, tol));
                    z  = sum(abs(ij)<am_dft.tiny);

                    if     all([t == 3, o == 3, z == 0]); bv_code = 1; % triclinic
                    elseif all([t == 2, o == 2, z == 2]); bv_code = 2; % hexagonal
                    elseif all([t == 3, o == 2, z == 2]); bv_code = 3; % monoclinic
                    elseif all([t == 2, o == 1, z == 3]); bv_code = 4; % tetragonal
                    elseif all([t == 3, o == 1, z == 3]); bv_code = 5; % orthorhombic
                    elseif all([t == 1, o == 1, z == 3]); bv_code = 6; % cubic
                    elseif all([t == 1, o == 1, z == 0]); bv_code = 7; % rhombahedral
                    else;                                 bv_code = 0; % unknown
                    end
                case 2
                    % Tinkham p 61
                    abc = am_cell.bas2abc(bas);
                    if      eq_(abc(1),abc(2), tol) &&  eq_(abc(1),abc(3), tol) &&  eq_(abc(2),abc(3), tol)
                        if  eq_(abc(4),abc(5), tol) &&  eq_(abc(5),abc(6), tol) &&  eq_(abc(6),abc(4), tol) &&  eq_(abc(4),90, tol)
                            bv_code = 6; % cubic        (a=b=c, alpha=beta=gamma=90)
                        else
                            bv_code = 7; % rhombohedral (a=b=c, alpha=beta=gamma!=90)
                        end
                    elseif ~eq_(abc(1),abc(2), tol) && ~eq_(abc(1),abc(3), tol) && ~eq_(abc(2),abc(3), tol)
                        if ~eq_(abc(4),abc(5), tol) && ~eq_(abc(4),abc(6), tol) && ~eq_(abc(5),abc(6), tol)
                            bv_code = 1; % triclinic    (a!=b!=c, alpha!=beta!=gamma!=90)
                        else
                            if  eq_(abc(4),abc(5), tol) &&  eq_(abc(4),abc(6), tol) &&  eq_(abc(5),abc(6), tol)
                            bv_code = 5; % orthorhombic (a!=b!=c, alpha=beta=gamma=90)
                            else
                            bv_code = 3; % monoclinic   (a!=b!=c, gamma!=(alpha=beta=90))
                            end
                        end
                    elseif  eq_(abc(1),abc(2), tol) && ~eq_(abc(1),abc(3), tol) && ~eq_(abc(2),abc(3), tol) && ...
                            eq_(abc(4),abc(5), tol) &&  eq_(abc(4),abc(6), tol) &&  eq_(abc(4),    90, tol)
                            bv_code = 4; % tetragonal   (a=b!=c, alpha=beta=gamma=90)
                    else
                            bv_code = 2; % hexagonal    (a=b!=c, alpha=beta=90,gamma=120)
                    end
                case 3
                    % run both algorithms and check for match
                    bv_code(1) = am_cell.identify_bravais(bas, tol, 1);
                    bv_code(2) = am_cell.identify_bravais(bas, tol, 2);
                    if bv_code(1) ~= bv_code(2)
                        warning('Algorithm mismatch. \n')
                    end
            end
        end

        function pg_code      = identify_pointgroup(R,tol)
            % pg_code = identify_pointgroup(R)
            %
            % Point symmetries in fractional coordinates so that they are nice integers
            % which can be easily classified.
            %
            %    element:         e    i  c_2  c_3  c_4  c_6  s_2  s_6  s_4  s_3
            %    trace:          +3   -3   -1    0   +1   +2   +1    0   -1   -2
            %    determinant:    +1   -1   +1   +1   +1   +1   -1   -1   -1   -1
            %
            %  The Mathematical Theory of Symmetry in Solids:  Representation Theory for
            %  Point Groups and Space Groups. 1 edition. Oxford?: New York: Oxford University
            %  Press, 2010. page 138, chartab 3.8.
            %
            %  Applied Group Theory: For Physicists and Chemists. Reissue edition.
            %  Mineola, New York: Dover Publications, 2015. page 20.
            %
            %  Casas, Ignasi, and Juan J. Perez. Modification to Flow Chart to
            %  Determine Point Groups. Journal of Chemical Education 69, no. 1
            %  (January 1, 1992): 83. doi:10.1021/ed069p83.2.
            %
            %  Breneman, G. L. ?Crystallographic Symmetry Point Group Notation
            %  Flow Chart.? Journal of Chemical Education 64, no. 3 (March 1, 1987):
            %  216. doi:10.1021/ed064p216.
            %
            %   pg_code --> point_group_name
            %     1 --> c_1       9 --> c_3      17 --> d_4      25 --> c_6v
            %     2 --> s_2      10 --> s_6      18 --> c_4v     26 --> d_3h
            %     3 --> c_2      11 --> d_3      19 --> d_2d     27 --> d_6h
            %     4 --> c_1h     12 --> c_3v     20 --> d_4h     28 --> t
            %     5 --> c_2h     13 --> d_3d     21 --> c_6      29 --> t_h
            %     6 --> d_2      14 --> c_4      22 --> c_3h     30 --> o
            %     7 --> c_2v     15 --> s_4      23 --> c_6h     31 --> t_d
            %     8 --> d_2h     16 --> c_4h     24 --> d_6      32 --> o_h
            %
            import am_lib.*
            %
            nsyms = size(R,3);
            % identify point symmetries
            ps_id = am_cell.identify_point_symmetries(R,tol);
            % count each type of symmetry
            nc2 = sum(ps_id==2);  % 'c_2'
            nc3 = sum(ps_id==3);  % 'c_3'
            nc4 = sum(ps_id==4);  % 'c_4'
            nc6 = sum(ps_id==5);  % 'c_6'
            ni  = sum(ps_id==6);  % 'i'
            ns2 = sum(ps_id==7);  % 's_2'
            ns4 = sum(ps_id==9);  % 's_4'
            % identify point group by comparing number and types of symmetries
            if         nsyms==1         ; pg_code=1 ; % c_1  
            elseif     nsyms==48        ; pg_code=32; % s_2  
            elseif     nsyms==16        ; pg_code=20; % c_2  
            elseif     nsyms==3         ; pg_code=9 ; % c_1h  
            elseif and(nsyms==2 , ni==1); pg_code=2 ; % c_2h  
            elseif and(nsyms==2 ,nc2==1); pg_code=3 ; % d_2  
            elseif and(nsyms==2 ,ns2==1); pg_code=4 ; % c_2v  
            elseif and(nsyms==4 , ni==1); pg_code=5 ; % d_2h  
            elseif and(nsyms==4 ,nc2==3); pg_code=6 ; % c_3  
            elseif and(nsyms==4 ,ns2==2); pg_code=7 ; % s_6  
            elseif and(nsyms==4 ,nc4==2); pg_code=14; % d_3    nc4 == 2 is correct, rather than nc4 == 1.
            elseif and(nsyms==4 ,ns4==2); pg_code=15; % c_3v  
            elseif and(nsyms==6 , ni==1); pg_code=10; % d_3d  
            elseif and(nsyms==6 ,nc2==3); pg_code=11; % c_4  
            elseif and(nsyms==6 ,ns2==3); pg_code=12; % s_4  
            elseif and(nsyms==6 ,nc2==1); pg_code=21; % c_4h  
            elseif and(nsyms==6 ,ns2==1); pg_code=22; % d_4  
            elseif and(nsyms==8 ,ns2==3); pg_code=8 ; % c_4v  
            elseif and(nsyms==8 ,ns2==1); pg_code=16; % d_2d  
            elseif and(nsyms==8 ,ns2==0); pg_code=17; % d_4h  
            elseif and(nsyms==8 ,ns2==4); pg_code=18; % c_6  
            elseif and(nsyms==8 ,ns2==2); pg_code=19; % c_3h  
            elseif and(nsyms==12,ns2==3); pg_code=13; % c_6h  
            elseif and(nsyms==12,ns2==1); pg_code=23; % d_6  
            elseif and(nsyms==12,nc2==7); pg_code=24; % c_6v  
            elseif and(nsyms==12,ns2==6); pg_code=25; % d_3h  
            elseif and(nsyms==12,ns2==4); pg_code=26; % d_6h  
            elseif and(nsyms==12,nc3==8); pg_code=28; % t  
            elseif and(nsyms==24,nc6==2); pg_code=27; % t_h  
            elseif and(nsyms==24, ni==1); pg_code=29; % o  
            elseif and(nsyms==24,nc4==6); pg_code=30; % t_d  
            elseif and(nsyms==24,ns4==6); pg_code=31; % o_h  
            else;                         pg_code=0 ; 
            end
        end

        function [c_code,C]   = identify_centering(pc)
            
            % this checks the lattice for primitive cell vectors which are expected for different centerings. If the cell is of a wierd shape or a supercell, this function will not work.
            
            import am_lib.* am_dft.*
            
            % symmetries
            [~,~,~,R] = get_symmetries(pc); tol=am_dft.tiny;
            
            % get number of point symmetries
            nRs = size(R,3);

            % % basic centring vectors
            T = [  0,1/2,1/2, 1/2, 2/3, 1/3, 2/3, 1/3, 1/3, 2/3;
                 1/2,  0,1/2, 1/2, 1/3, 2/3, 1/3, 2/3, 1/3, 2/3;
                 1/2,1/2,  0, 1/2, 1/3, 2/3,   0,   0, 1/3, 2/3;
                   1   2,  3,   4,   5,   5,   6,   7,   8,   9];

            % substrate centering vecotrs
            X = [ 0, 0,-1,-1,-1, 0,-1, 0;
                  0,-1, 0,-1, 0,-1,-1, 0;
                 -1, 0, 0, 0,-1,-1,-1, 0;
                  0, 0, 0, 0, 0, 0, 0, 0];

            % subtract 1 from nonzero values and append identity
            C = osum_(T,X,2); C = [C(:,~any(eq_(C,-1),1)),[eye(3);0,0,0]]; nCs = size(C,2);

            % eliminate vectors that are collinear
            ex_ = false(nCs,nCs); is_collinear_ = @(X,Y) all(eq_(cross(X(1:3),Y(1:3)),0)); 
            for i = 1:nCs; for j = i:nCs; if i ~= j
                ex_(i,j) = is_collinear_(C(:,i),C(:,j));
            end; end; end
            C = C(:,~any(ex_,1)); nCs = size(C,2);
            
            % find out which subset of vectors V preserve periodic boundary conditions
            check_ = @(A) all(all(abs(A)<tol,1),2); ex_=false(1,nCs); 
            X_ = @(tau,species) sortc_([mod_(tau);species]); X = X_(pc.tau(:,:,1),pc.species);
            for j = 1:nCs; ex_(j) = check_( X_(pc.tau(1:3,:,1)-C(1:3,j),pc.species)-X ); end
            C=C(:,ex_); nCs = size(C,2);

            % rank vectors in increasing length
            C = C(:,rankc_(normc_(C(1:3,:))));

            % get all combinations of vectors
            ijk = am_lib.permn_(1:nCs,3).'; 
            ex_ = ~any([ijk(1,:)==ijk(2,:); 
                        ijk(2,:)==ijk(3,:); 
                        ijk(3,:)==ijk(1,:)],1);
            ijk = ijk(:,ex_); nijks = sum(ex_); 

            % check which combination of three vectors produces symmetries which have only integer components
            for i = 1:nijks         
                Z = C(1:3,ijk(1:3,i)); d = det(Z);
                % check collinearity
                if eq_(d,0); continue; end
                % flip to get right-handed if necessary
                if d<0; Z = fliplr(Z); d = -d; ijk(1:3,i) = flipud(ijk(1:3,i)); end
                % check multiplicity (number of nodes per unit cell)
                if ~any(eq_(1/d,[1,2,3,4])); continue; end 
                % check that symmetries contain only integers
                go = true; iZ=inv(Z);
                for j = 1:nRs
                    Rp = Z*R(:,:,j)*iZ;
                    if any_( ~any(eq_(abs(Rp(:)),[0,1]),2) )
                        go = false; break; 
                    end
                end
                if ~go; continue; end
                % transformation found: separate vector code and centering matrix
                v_code = C(4,ijk(:,i)); C = C(1:3,ijk(:,i)); break;
            end

            % centering matrix
            if     all(any(v_code==[ 1; 2; 3],1)); c_code = 1; % face-centered
            elseif all(any(v_code==[ 1; 0; 0],1)); c_code = 2; % A-centered 
            elseif all(any(v_code==[ 2; 0; 0],1)); c_code = 3; % B-centered 
            elseif all(any(v_code==[ 3; 0; 0],1)); c_code = 4; % C-centered 
            elseif all(v_code==4);                 c_code = 5; % body-centered 
            elseif all(any(v_code==[ 5; 6; 0],1)); c_code = 6; % rhombohedral (hexagonal)
            elseif all(any(v_code==[ 7; 8; 0],1)); c_code = 7; % hexagonal
            elseif all(any(v_code==[ 9;10; 0],1)); c_code = 8; % rhombohedral
            else; c_code = 0; % unknown
            end
        end
        
        function sg_code      = identify_spacegroup(pg_code,tol)
            % NOTE: THIS CODE IS INCOMPLETE. the idea was to narrow down
            % the possible space groups by requiring that 1) the point
            % groups match, 2) the number of symmetries match, and 3) the
            % number of classes match. While matching point groups is good
            % aving classes and symmmetries match do not seem to be that
            % great of criteria...
            %
            % generate databases using:
            % import am_lib.*
            import am_dft.*
            % for i = 1:237
            %     % generate space symmetries and multiplication table
            %     S{i} = generate_sg(i);
            %     % determine pg corresponding to each sgi %4i , ...\n',cell2mat(pg))
            %     pg{i} = identify_pointgroup( reshape(uniquec_( reshape(S{i}(1:3,1:3,:),9,[]) ),3,3,[]) );
            %     % count number of symmetries
            %     nsyms{i} = size(S{i},3);
            %     % indentify number of classes
            %     nclasses{i} = numel(unique(identify_classes(MT{i})));
            % end
            % fprintf('pg_database = [ ... \n'); fprintf('%4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i , ...\n',cell2mat(pg)); fprintf(']; \n');
            % fprintf('nsyms_database = [ ... \n');fprintf('%4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i , ...\n',cell2mat(nsyms)); fprintf(']; \n');
            % fprintf('nclasses_database = [ ... \n');fprintf('%4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i , ...\n',cell2mat(nclasses)); fprintf(']; \n');
            %

            pg_database = [ ...
               1    2    3    3    3    4    4    4    4    5    5    5 , ...
               5    5    5    6    6    6    6    6    6    6    6    6 , ...
               7    7    7    7    7    7    7    7    7    7    7    7 , ...
               7    7    7    7    7    7    7    7    7    7    8    8 , ...
               8    8    8    8    8    8    8    8    8    8    8    8 , ...
               8    8    8    8    8    8    8    8    8    8    8    8 , ...
               8    8   14   14   14   14   14   14   15   15   16   16 , ...
              16   16   16   16   17   17   17   17   17   17   17   17 , ...
              17   17   18   18   18   18   18   18   18   18   18   18 , ...
              18   18   19   19   19   19   19   19   19   19   19   19 , ...
              19   19   20   20   20   20   20   20   20   20   20   20 , ...
              20   20   20   20   20   20   20   20   20   20    9    9 , ...
               9    9   10   10   11   11   11   11   11   11   11   12 , ...
              12   12   12   12   12   13   13   13   13   13   13   21 , ...
              21   21   21   21   21   22   23   23   24   24   24   24 , ...
              24   24   25   25   25   25   26   26   26   26   27   27 , ...
              27   27   28   28   28   28   28   29   29   29   29   29 , ...
              29   29   30   30   30   30   30   30   30   30   31   31 , ...
              31   31   31   31   32   32   32   32   32   32   32   32 , ...
              32   32    9   10   11   12   12   13   13 ];
            nsyms_database = [ ...
               1    2    2    2    4    2    2    4    4    4    4    8 , ...
               4    4    8    4    4    4    4    8    8   16    8    8 , ...
               4    4    4    4    4    4    4    4    4    4    8    8 , ...
               8    8    8    8    8   16   16    8    8    8    8    8 , ...
               8    8    8    8    8    8    8    8    8    8    8    8 , ...
               8    8   16   16   16   16   16   16   32   32   16   16 , ...
              16   16    4    4    4    4    8    8    4    8    8    8 , ...
               8    8   16   16    8    8    8    8    8    8    8    8 , ...
              16   16    8    8    8    8    8    8    8    8   16   16 , ...
              16   16    8    8    8    8    8    8    8    8   16   16 , ...
              16   16   16   16   16   16   16   16   16   16   16   16 , ...
              16   16   16   16   16   16   32   32   32   32    3    3 , ...
               3    9    6   18    6    6    6    6    6    6   18    6 , ...
               6    6    6   18   18   12   12   12   12   36   36    6 , ...
               6    6    6    6    6    6   12   12   12   12   12   12 , ...
              12   12   12   12   12   12   12   12   12   12   24   24 , ...
              24   24   12   48   24   12   24   24   24   96   96   48 , ...
              24   48   24   24   96   96   48   24   24   48   24   96 , ...
              48   24   96   48   48   48   48   48  192  192  192  192 , ...
              96   96    3    6    6    6    6   12   12 ]; %#ok<NASGU>
            nclasses_database = [ ...
               1    2    2    2    4    2    2    4    4    4    4    8 , ...
               4    4    8    4    4    4    4    8    8   16    8    8 , ...
               4    4    4    4    4    4    4    4    4    4    8    8 , ...
               8    8    8    8    8   16   10    8    8    8    8    8 , ...
               8    8    8    8    8    8    8    8    8    8    8    8 , ...
               8    8   16   16   16   16   16   16   32   14   16   16 , ...
              16   16    4    4    4    4    8    8    4    8    8    8 , ...
               8    8   16   10    5    5    5    5    5    5    5    5 , ...
              10   10    5    5    5    5    5    5    5    5   10   10 , ...
              10   10    5    5    5    5    5    5    5    5   10   10 , ...
              10   10   10   10   10   10   10   10   10   10   10   10 , ...
              10   10   10   10   10   10   20   20   14   14    3    3 , ...
               3    9    6    9    3    3    3    3    3    3    6    3 , ...
               3    3    3    9    9    6    6    6    6    9    9    6 , ...
               6    6    6    6    6    6   12   12    6    6    6    6 , ...
               6    6    6    6    6    6    6    6    6    6   12   12 , ...
              12   12    4    8    8    4    8    8    8   16   10   16 , ...
               8   16    5    5   10   10   10    5    5   10    5   10 , ...
              10    5   10   10   10   10   10   10   20   20   14   14 , ...
              20   14    3    6    3    3    3    6    6 ]; %#ok<NASGU>
            sg_code = find( (pg_database==pg_code) ); 
        end

        function lg_code      = identify_laue(pg_code)
            % pg={'c_1' ,'s_2' ,'c_2' ,'c_1h','c_2h','d_2' ,'c_2v','d_2h', ...
            %     'c_3' ,'s_6' ,'d_3' ,'c_3v','d_3d','c_4' ,'s_4' ,'c_4h', ...
            %     'd_4' ,'c_4v','d_2d','d_4h','c_6' ,'c_3h','c_6h','d_6' , ...
            %     'c_6v','d_3h','d_6h','t'   ,'t_h' ,'o'   ,'t_d' ,'o_h'};
            
            % laue group dataset
            lg=[1,1,2,2,2,3,3,3, ...
                4,4,5,5,5,6,6,6, ...
                7,7,7,7,8,8,8,9, ...
                9,9,9,10,10,11,11,11];
            % print point group name
            lg_code = lg(pg_code);
        end
        
        function c_name       = decode_centering(c_code)
            % F face-centered
            % A-centered 
            % B-centered 
            % C-centered 
            % I body-centered 
            % R rhombohedral (hexagonal)
            % H hexagonal
            % D rhombohedral	
            c_database={'F','A','B','C','I','R','H','D'}; 
            c_name = c_database{c_code};
        end
        
        function bv_name      = decode_bravais(bv_code)
            % bravais lattices dataset
            bv={'triclinic','hexagonal','monoclinic','tetragonal',...
                'orthorhombic','cubic','rhombahedral'};
            % print bravais lattice type
            bv_name = bv{bv_code};
        end

        function brv_name     = decode_holohodry(pg_code)
            % point group dataset
            brav={'triclinic','','','','monoclinic','','','orthorhombic', ...
                  '','','','','trigonal','','','','','','','tetragonal',...
                  '','','','','','','hexagonal','','','','','cubic'};
            % print point group name
            brv_name = brav{pg_code};
        end
        
        function lg_name      = decode_laue(lg_code)
            import am_dft.*
            % laue group dataset
            lg = {'-1' ,'2/m' ,'mmm' ,'-3','-3m','4/m' ,'4/mmm','6/m', '6/mmm' ,'m-3' , 'm-3m'};
            % print point group name
            lg_name = lg{lg_code};
        end

        function pg_name      = decode_pg(pg_code)
            % point group dataset
            pg={'c_1' ,'s_2' ,'c_2' ,'c_1h','c_2h','d_2' ,'c_2v','d_2h', ...
                'c_3' ,'s_6' ,'d_3' ,'c_3v','d_3d','c_4' ,'s_4' ,'c_4h', ...
                'd_4' ,'c_4v','d_2d','d_4h','c_6' ,'c_3h','c_6h','d_6' , ...
                'c_6v','d_3h','d_6h','t'   ,'t_h' ,'o'   ,'t_d' ,'o_h'};
            % print point group name
            pg_name = pg{pg_code};
        end

        function sg_name      = decode_sg(sg_code)
            sg_name_database = { ...
                ' P  1      ' ,' P -1      ', ... % monoclinic space groups
                ' P 2       ' ,' P 21      ' ,' C 2       ' ,' P m       ', ...
                ' P c       ' ,' C m       ' ,' C c       ' ,' P 2/m     ', ...
                ' P 21/m    ' ,' C 2/m     ' ,' P 2/c     ' ,' P 21/c    ', ...
                ' C 2/c     ', ...                                              % orthorhombic space groups
                ' P 2 2 2   ' ,' P 2 2 21  ' ,' P 21 21 2 ' ,' P 21 21 21', ...
                ' C 2 2 21  ' ,' C 2 2 2   ' ,' F 2 2 2   ' ,' I 2 2 2   ', ...
                ' I 21 21 21' ,' P m m 2   ' ,' P m c 21  ' ,' P c c 2   ', ...
                ' P m a 2   ' ,' P c a 21  ' ,' P n c 2   ' ,' P m n 21  ', ...
                ' P b a 2   ' ,' P n a 21  ' ,' P n n 2   ' ,' C m m 2   ', ...
                ' C m c 21  ' ,' C c c 2   ' ,' A m m 2   ' ,' A b m 2   ', ...
                ' A m a 2   ' ,' A b a 2   ' ,' F m m 2   ' ,' F d d 2   ', ...
                ' I m m 2   ' ,' I b a 2   ' ,' I m a 2   ' ,' P m m m   ', ...
                ' P n n n   ' ,' P c c m   ' ,' P b a n   ' ,' P m m a   ', ...
                ' P n n a   ' ,' P m n a   ' ,' P c c a   ' ,' P b a m   ', ...
                ' P c c n   ' ,' P b c m   ' ,' P n n m   ' ,' P m m n   ', ...
                ' P b c n   ' ,' P b c a   ' ,' P n m a   ' ,' C m c m   ', ...
                ' C m c a   ' ,' C m m m   ' ,' C c c m   ' ,' C m m a   ', ...
                ' C c c a   ' ,' F m m m   ' ,' F d d d   ' ,' I m m m   ', ...
                ' I b a m   ' ,' I b c a   ' ,' I m m a   ', ...                % tetragonal space groups
                ' P 4       ' ,' P 41      ' ,' P 42      ' ,' P 43      ', ...
                ' I 4       ' ,' I 41      ' ,' P -4      ' ,' I -4      ', ...
                ' P 4/m     ' ,' P 42/m    ' ,' P 4/n     ' ,' P 42/n    ', ...
                ' I 4/m     ' ,' I 41/a    ' ,' P 4 2 2   ' ,' P 4 21 2  ', ...
                ' P 41 2 2  ' ,' P 41 21 2 ' ,' P 42 2 2  ' ,' P 42 21 2 ', ...
                ' P 43 2 2  ' ,' P 43 21 2 ' ,' I 4 2 2   ' ,' I 41 2 2  ', ...
                ' P 4 m m   ' ,' P 4 b m   ' ,' P 42 c m  ' ,' P 42 n m  ', ...
                ' P 4 c c   ' ,' P 4 n c   ' ,' P 42 m c  ' ,' P 42 b c  ', ...
                ' I 4 m m   ' ,' I 4 c m   ' ,' I 41 m d  ' ,' I 41 c d  ', ...
                ' P -4 2 m  ' ,' P -4 2 c  ' ,' P -4 21 m ' ,' P -4 21 c ', ...
                ' P -4 m 2  ' ,' P -4 c 2  ' ,' P -4 b 2  ' ,' P -4 n 2  ', ...
                ' I -4 m 2  ' ,' I -4 c 2  ' ,' I -4 2 m  ' ,' I -4 2 d  ', ...
                ' P 4/m m m ' ,' P 4/m c c ' ,' P 4/n b m ' ,' P 4/n n c ', ...
                ' P 4/m b m ' ,' P 4/m n c ' ,' P 4/n m m ' ,' P 4/n c c ', ...
                ' P 42/m m c' ,' P 42/m c m' ,' P 42/n b c' ,' P 42/n n m', ...
                ' P 42/m b c' ,' P 42/m n m' ,' P 42/n m c' ,' P 42/n c m', ...
                ' I 4/m m m ' ,' I 4/m c m ' ,' I 41/a m d' ,' I 41/a c d', ... % rhombohedral space groups
                ' P 3       ' ,' P 31      ' ,' P 32      ' ,' R 3       ', ...
                ' P -3      ' ,' R -3      ' ,' P 3 1 2   ' ,' P 3 2 1   ', ...
                ' P 31 1 2  ' ,' P 31 2 1  ' ,' P 32 1 2  ' ,' P 32 2 1  ', ...
                ' R 3 2     ' ,' P 3 m 1   ' ,' P 3 1 m   ' ,' P 3 c 1   ', ...
                ' P 3 1 c   ' ,' R 3 m     ' ,' R 3 c     ' ,' P -3 1 m  ', ...
                ' P -3 1 c  ' ,' P -3 m 1  ' ,' P -3 c 1  ' ,' R -3 m    ', ...
                ' R -3 c    ', ...                                              % hexagonal space groups
                ' P 6       ' ,' P 61      ' ,' P 65      ' ,' P 62      ', ...
                ' P 64      ' ,' P 63      ' ,' P -6      ' ,' P 6/m     ', ...
                ' P 63/m    ' ,' P 6 2 2   ' ,' P 61 2 2  ' ,' P 65 2 2  ', ...
                ' P 62 2 2  ' ,' P 64 2 2  ' ,' P 63 2 2  ' ,' P 6 m m   ', ...
                ' P 6 c c   ' ,' P 63 c m  ' ,' P 63 m c  ' ,' P -6 m 2  ', ...
                ' P -6 c 2  ' ,' P -6 2 m  ' ,' P -6 2 c  ' ,' P 6/m m m ', ...
                ' P 6/m c c ' ,' P 63/m c m' ,' P 63/m m c', ...                % cubic space groups
                ' P 2 3     ' ,' F 2 3     ' ,' I 2 3     ' ,' P 21 3    ', ...
                ' I 21 3    ' ,' P m 3     ' ,' P n 3     ' ,' F m 3     ', ...
                ' F d 3     ' ,' I m 3     ' ,' P a 3     ' ,' I a 3     ', ...
                ' P 4 3 2   ' ,' P 42 3 2  ' ,' F 4 3 2   ' ,' F 41 3 2  ', ...
                ' I 4 3 2   ' ,' P 43 3 2  ' ,' P 41 3 2  ' ,' I 41 3 2  ', ...
                ' P -4 3 m  ' ,' F -4 3 m  ' ,' I -4 3 m  ' ,' P -4 3 n  ', ...
                ' F -4 3 c  ' ,' I -4 3 d  ' ,' P m 3 m   ' ,' P n 3 n   ', ...
                ' P m 3 n   ' ,' P n 3 m   ' ,' F m 3 m   ' ,' F m 3 c   ', ...
                ' F d 3 m   ' ,' F d 3 c   ' ,' I m 3 m   ' ,' I a 3 d   ', ... % trigonal groups rhombohedral setting
                ' R 3   |146' ,' R -3  |148' ,' R 3 2 |155' ,' R 3 m |160', ...
                ' R 3 c |161' ,' R -3 m|166' ,' R -3 c|167'};
            for i = 1:numel(sg_code)
                sg_name{i} = sprintf('%s (%i)', strtrim(sg_name_database{sg_code(i)}), sg_code(i)); 
            end
        end

        function ps_name      = decode_ps(ps_code)
            nsyms = numel(ps_code);ps_name=cell(1,nsyms);
            for i = 1:nsyms
                switch abs(ps_code(i))
                    case 1;   ps_name{i}='e';
                    case 2;   ps_name{i}='c_2';
                    case 3;   ps_name{i}='c_3';
                    case 4;   ps_name{i}='c_4';
                    case 5;   ps_name{i}='c_6';
                    case 6;   ps_name{i}='i';
                    case 7;   ps_name{i}='s_2';
                    case 8;   ps_name{i}='s_6';
                    case 9;   ps_name{i}='s_4';
                    case 10;  ps_name{i}='s_3';
                    otherwise;ps_name{i}='';
                end
                % append double group
                if ps_code(i)<0
                    ps_name{i} = ['R',ps_name{i}];
                end
            end
        end

        function ps_name_long = get_long_ps_name(R)
            nsyms=size(R,3); 
            ps_id=am_dft.identify_point_symmetries(R);
            ps_name=am_dft.decode_ps(ps_id);
            ps_axis=am_lib.rnd_(am_lib.R_axis_(R));
            syms x y z
            for i = 1:nsyms
                ps_name_long{i} = sprintf('%s',strtrim(ps_name{i}));
                if abs(ps_id(i))~=1 && abs(ps_id(i))~=6
                    ax = ps_axis(:,i).'; if sum(am_lib.lt_(ax,0))>sum(am_lib.gt_(ax,0)); ax=-ax; end
                    ax_name = sprintf('%s',ax*[x;y;z]); ax_name=strrep(ax_name,' ',''); ax_name=strrep(ax_name,'+',''); 
                    ps_name_long{i} = sprintf('%s(%s)',ps_name_long{i}, ax_name );
                end
            end
        end

        function ss_name_long = get_long_ss_name(S)
            
            import am_dft.* am_lib.*
            
            if isempty(S); ss_name_long=[]; return; end
            
            if iscell(S)
                ss_name_long = get_long_ss_name(S{1});
                for i = 1:size(S{1},3)
                    pp_name = [sprintf('%i',S{2}(1,i)), sprintf(',%i',S{2}(2:end,i))];
                    ss_name_long{i} = sprintf('%s (%s)', ss_name_long{i}, pp_name);
                end
            else
                ss_name_long = get_long_ps_name(S(1:3,1:3,:));
                nsyms=size(S,3); E = find(all(all(eq_(S,eye(4)),1),2));
                for i = 1:nsyms
                    if i==E; continue; end
                    t = S(1:3,4,i);
                    t_name = sprintf('[%3.2f %3.2f %3.2f]',t);
                    ss_name_long{i} = sprintf('%10s %s',ss_name_long{i}, t_name );
                end
            end
        end

        function S            = generate_sg(sg_code,from_memory)

            import am_lib.* am_dft.*

            if nargin<2; from_memory=false; end

            if from_memory
                % ~ 500x faster for large space groups
                S = load_sg_symmetries(sg_code);
                S = reshape(S,4,4,[]);
            else
                
                % load recipe
                recipe = get_recipe(sg_code);
                
                % allocate space
                S = zeros(4,4,sscanf(recipe(2),'%i')+double(recipe(1)=='1'));

                % mix ingredients
                nsyms = 0; k = 0;
                nsyms = nsyms+1; S(1:4,1:4,nsyms) = eye(4);
                %
                k = k+1;
                if recipe(1) == '1'
                    nsyms = nsyms+1; S(1:3,1:3,nsyms) = -eye(3); S(4,4,nsyms)=1;
                end
                k = k+1; ngens = sscanf(recipe(2),'%i');
                for i = 1:ngens
                    nsyms = nsyms+1;
                for j = 1:4
                    k=k+1;
                    switch j
                        case 1
                            S(1:3,1:3,nsyms) = decode_R(recipe(k)); S(4,4,nsyms)=1;
                        otherwise
                            S(j-1,4,nsyms) = mod(decode_T(recipe(k)),1);
                    end
                end
                end

                S = complete_group(S);
            end

            function recipe  = get_recipe(sg_id)
                sg_recipe_database = {...
                    '000                                     ','100                                     ','01cOOO0                                 ', ...
                    '01cODO0                                 ','02aDDOcOOO0                             ','01jOOO0                                 ', ...
                    '01jOOD0                                 ','02aDDOjOOO0                             ','02aDDOjOOD0                             ', ...
                    '11cOOO0                                 ','11cODO0                                 ','12aDDOcOOO0                             ', ...
                    '11cOOD0                                 ','11cODD0                                 ','12aDDOcOOD0                             ', ...
                    '02bOOOcOOO0                             ','02bOODcOOD0                             ','02bOOOcDDO0                             ', ...
                    '02bDODcODD0                             ','03aDDObOODcOOD0                         ','03aDDObOOOcOOO0                         ', ...
                    '04aODDaDODbOOOcOOO0                     ','03aDDDbOOOcOOO0                         ','03aDDDbDODcODD0                         ', ...
                    '02bOOOjOOO0                             ','02bOODjOOD0                             ','02bOOOjOOD0                             ', ...
                    '02bOOOjDOO0                             ','02bOODjDOO0                             ','02bOOOjODD0                             ', ...
                    '02bDODjDOD0                             ','02bOOOjDDO0                             ','02bOODjDDO0                             ', ...
                    '02bOOOjDDD0                             ','03aDDObOOOjOOO0                         ','03aDDObOODjOOD0                         ', ...
                    '03aDDObOOOjOOD0                         ','03aODDbOOOjOOO0                         ','03aODDbOOOjODO0                         ', ...
                    '03aODDbOOOjDOO0                         ','03aODDbOOOjDDO0                         ','04aODDaDODbOOOjOOO0                     ', ...
                    '04aODDaDODbOOOjBBB0                     ','03aDDDbOOOjOOO0                         ','03aDDDbOOOjDDO0                         ', ...
                    '03aDDDbOOOjDOO0                         ','12bOOOcOOO0                             ','03bOOOcOOOhDDD1BBB                      ', ...
                    '12bOOOcOOD0                             ','03bOOOcOOOhDDO1BBO                      ','12bDOOcOOO0                             ', ...
                    '12bDOOcDDD0                             ','12bDODcDOD0                             ','12bDOOcOOD0                             ', ...
                    '12bOOOcDDO0                             ','12bDDOcODD0                             ','12bOODcODD0                             ', ...
                    '12bOOOcDDD0                             ','03bOOOcDDOhDDO1BBO                      ','12bDDDcOOD0                             ', ...
                    '12bDODcODD0                             ','12bDODcODO0                             ','13aDDObOODcOOD0                         ', ...
                    '13aDDObODDcODD0                         ','13aDDObOOOcOOO0                         ','13aDDObOOOcOOD0                         ', ...
                    '13aDDObODOcODO0                         ','04aDDObDDOcOOOhODD1OBB                  ','14aODDaDODbOOOcOOO0                     ', ...
                    '05aODDaDODbOOOcOOOhBBB1ZZZ              ','13aDDDbOOOcOOO0                         ','13aDDDbOOOcDDO0                         ', ...
                    '13aDDDbDODcODD0                         ','13aDDDbODOcODO0                         ','02bOOOgOOO0                             ', ...
                    '02bOODgOOB0                             ','02bOOOgOOD0                             ','02bOODgOOF0                             ', ...
                    '03aDDDbOOOgOOO0                         ','03aDDDbDDDgODB0                         ','02bOOOmOOO0                             ', ...
                    '03aDDDbOOOmOOO0                         ','12bOOOgOOO0                             ','12bOOOgOOD0                             ', ...
                    '03bOOOgDDOhDDO1YBO                      ','03bOOOgDDDhDDD1YYY                      ','13aDDDbOOOgOOO0                         ', ...
                    '04aDDDbDDDgODBhODB1OYZ                  ','03bOOOgOOOcOOO0                         ','03bOOOgDDOcDDO0                         ', ...
                    '03bOODgOOBcOOO0                         ','03bOODgDDBcDDB0                         ','03bOOOgOODcOOO0                         ', ...
                    '03bOOOgDDDcDDD0                         ','03bOODgOOFcOOO0                         ','03bOODgDDFcDDF0                         ', ...
                    '04aDDDbOOOgOOOcOOO0                     ','04aDDDbDDDgODBcDOF0                     ','03bOOOgOOOjOOO0                         ', ...
                    '03bOOOgOOOjDDO0                         ','03bOOOgOODjOOD0                         ','03bOOOgDDDjDDD0                         ', ...
                    '03bOOOgOOOjOOD0                         ','03bOOOgOOOjDDD0                         ','03bOOOgOODjOOO0                         ', ...
                    '03bOOOgOODjDDO0                         ','04aDDDbOOOgOOOjOOO0                     ','04aDDDbOOOgOOOjOOD0                     ', ...
                    '04aDDDbDDDgODBjOOO0                     ','04aDDDbDDDgODBjOOD0                     ','03bOOOmOOOcOOO0                         ', ...
                    '03bOOOmOOOcOOD0                         ','03bOOOmOOOcDDO0                         ','03bOOOmOOOcDDD0                         ', ...
                    '03bOOOmOOOjOOO0                         ','03bOOOmOOOjOOD0                         ','03bOOOmOOOjDDO0                         ', ...
                    '03bOOOmOOOjDDD0                         ','04aDDDbOOOmOOOjOOO0                     ','04aDDDbOOOmOOOjOOD0                     ', ...
                    '04aDDDbOOOmOOOcOOO0                     ','04aDDDbOOOmOOOcDOF0                     ','13bOOOgOOOcOOO0                         ', ...
                    '13bOOOgOOOcOOD0                         ','04bOOOgOOOcOOOhDDO1YYO                  ','04bOOOgOOOcOOOhDDD1YYY                  ', ...
                    '13bOOOgOOOcDDO0                         ','13bOOOgOOOcDDD0                         ','04bOOOgDDOcDDOhDDO1YBO                  ', ...
                    '04bOOOgDDOcDDDhDDO1YBO                  ','13bOOOgOODcOOO0                         ','13bOOOgOODcOOD0                         ', ...
                    '04bOOOgDDDcOODhDDD1YBY                  ','04bOOOgDDDcOOOhDDD1YBY                  ','13bOOOgOODcDDO0                         ', ...
                    '13bOOOgDDDcDDD0                         ','04bOOOgDDDcDDDhDDD1YBY                  ','04bOOOgDDDcDDOhDDD1YBY                  ', ...
                    '14aDDDbOOOgOOOcOOO0                     ','14aDDDbOOOgOOOcOOD0                     ','05aDDDbDDDgODBcDOFhODB1OBZ              ', ...
                    '05aDDDbDDDgODBcDOBhODB1OBZ              ','01nOOO0                                 ','01nOOC0                                 ', ...
                    '01nOOE0                                 ','02aECCnOOO0                             ','11nOOO0                                 ', ...
                    '12aECCnOOO0                             ','02nOOOfOOO0                             ','02nOOOeOOO0                             ', ...
                    '02nOOCfOOE0                             ','02nOOCeOOO0                             ','02nOOEfOOC0                             ', ...
                    '02nOOEeOOO0                             ','03aECCnOOOeOOO0                         ','02nOOOkOOO0                             ', ...
                    '02nOOOlOOO0                             ','02nOOOkOOD0                             ','02nOOOlOOD0                             ', ...
                    '03aECCnOOOkOOO0                         ','03aECCnOOOkOOD0                         ','12nOOOfOOO0                             ', ...
                    '12nOOOfOOD0                             ','12nOOOeOOO0                             ','12nOOOeOOD0                             ', ...
                    '13aECCnOOOeOOO0                         ','13aECCnOOOeOOD0                         ','02nOOObOOO0                             ', ...
                    '02nOOCbOOD0                             ','02nOOEbOOD0                             ','02nOOEbOOO0                             ', ...
                    '02nOOCbOOO0                             ','02nOOObOOD0                             ','02nOOOiOOO0                             ', ...
                    '12nOOObOOO0                             ','12nOOObOOD0                             ','03nOOObOOOeOOO0                         ', ...
                    '03nOOCbOODeOOC0                         ','03nOOEbOODeOOE0                         ','03nOOEbOOOeOOE0                         ', ...
                    '03nOOCbOOOeOOC0                         ','03nOOObOODeOOO0                         ','03nOOObOOOkOOO0                         ', ...
                    '03nOOObOOOkOOD0                         ','03nOOObOODkOOD0                         ','03nOOObOODkOOO0                         ', ...
                    '03nOOOiOOOkOOO0                         ','03nOOOiOODkOOD0                         ','03nOOOiOOOeOOO0                         ', ...
                    '03nOOOiOODeOOO0                         ','13nOOObOOOeOOO0                         ','13nOOObOOOeOOD0                         ', ...
                    '13nOOObOODeOOD0                         ','13nOOObOODeOOO0                         ','03bOOOcOOOdOOO0                         ', ...
                    '05aODDaDODbOOOcOOOdOOO0                 ','04aDDDbOOOcOOOdOOO0                     ','03bDODcODDdOOO0                         ', ...
                    '04aDDDbDODcODDdOOO0                     ','13bOOOcOOOdOOO0                         ','04bOOOcOOOdOOOhDDD1YYY                  ', ...
                    '15aODDaDODbOOOcOOOdOOO0                 ','06aODDaDODbOOOcOOOdOOOhBBB1ZZZ          ','14aDDDbOOOcOOOdOOO0                     ', ...
                    '13bDODcODDdOOO0                         ','14aDDDbDODcODDdOOO0                     ','04bOOOcOOOdOOOeOOO0                     ', ...
                    '04bOOOcOOOdOOOeDDD0                     ','06aODDaDODbOOOcOOOdOOOeOOO0             ','06aODDaDODbODDcDDOdOOOeFBF0             ', ...
                    '05aDDDbOOOcOOOdOOOeOOO0                 ','04bDODcODDdOOOeBFF0                     ','04bDODcODDdOOOeFBB0                     ', ...
                    '05aDDDbDODcODDdOOOeFBB0                 ','04bOOOcOOOdOOOlOOO0                     ','06aODDaDODbOOOcOOOdOOOlOOO0             ', ...
                    '05aDDDbOOOcOOOdOOOlOOO0                 ','04bOOOcOOOdOOOlDDD0                     ','06aODDaDODbOOOcOOOdOOOlDDD0             ', ...
                    '05aDDDbDODcODDdOOOlBBB0                 ','14bOOOcOOOdOOOeOOO0                     ','05bOOOcOOOdOOOeOOOhDDD1YYY              ', ...
                    '14bOOOcOOOdOOOeDDD0                     ','05bOOOcOOOdOOOeDDDhDDD1YYY              ','16aODDaDODbOOOcOOOdOOOeOOO0             ', ...
                    '16aODDaDODbOOOcOOOdOOOeDDD0             ','07aODDaDODbODDcDDOdOOOeFBFhBBB1ZZZ      ','07aODDaDODbODDcDDOdOOOeFBFhFFF1XXX      ', ...
                    '15aDDDbOOOcOOOdOOOeOOO0                 ','15aDDDbDODcODDdOOOeFBB0                 ','01dOOO0                                 ', ...
                    '11dOOO0                                 ','02dOOOfOOO0                             ','02dOOOlOOO0                             ', ...
                    '02dOOOlDDD0                             ','12dOOOfOOO0                             ','12dOOOfDDD0                             '};
                recipe = sg_recipe_database{sg_id};
            end
            function R       = decode_R(str_r)
                switch str_r
                case 'a'; R = [  1  0  0;  0  1  0;  0  0  1]; % a
                case 'b'; R = [ -1  0  0;  0 -1  0;  0  0  1]; % b
                case 'c'; R = [ -1  0  0;  0  1  0;  0  0 -1]; % c
                case 'd'; R = [  0  0  1;  1  0  0;  0  1  0]; % d
                case 'e'; R = [  0  1  0;  1  0  0;  0  0 -1]; % e
                case 'f'; R = [  0 -1  0; -1  0  0;  0  0 -1]; % f
                case 'g'; R = [  0 -1  0;  1  0  0;  0  0  1]; % g
                case 'h'; R = [ -1  0  0;  0 -1  0;  0  0 -1]; % h
                case 'i'; R = [  1  0  0;  0  1  0;  0  0 -1]; % i
                case 'j'; R = [  1  0  0;  0 -1  0;  0  0  1]; % j
                case 'k'; R = [  0 -1  0; -1  0  0;  0  0  1]; % k
                case 'l'; R = [  0  1  0;  1  0  0;  0  0  1]; % l
                case 'm'; R = [  0  1  0; -1  0  0;  0  0 -1]; % m
                case 'n'; R = [  0 -1  0;  1 -1  0;  0  0  1]; % n
                end
            end
            function T       = decode_T(str_t)
                switch str_t
                    case 'A'; T = 1/6; % A
                    case 'B'; T = 1/4; % B
                    case 'C'; T = 1/3; % C
                    case 'D'; T = 1/2; % D
                    case 'E'; T = 2/3; % E
                    case 'F'; T = 3/4; % F
                    case 'G'; T = 5/6; % G
                    case 'O'; T =   0; % O
                    case 'X'; T =-3/8; % X
                    case 'Y'; T =-1/4; % Y
                    case 'Z'; T =-1/8; % Z
                end
            end
        end

        function [R,W]        = generate_pg(pg_code,from_memory)
            %     1.   c_1        9.   c_3        17. d_4       25. c_6v
            %     2.   s_2        10.  s_6        18. c_4v      26. d_3h
            %     3.   c_2        11.  d_3        19. d_2d      27. d_6h
            %     4.   c_1h       12.  c_3v       20. d_4h      28. t
            %     5.   c_2h       13.  d_3d       21. c_6       29. t_h
            %     6.   d_2        14.  c_4        22. c_3h      30. o
            %     7.   c_2v       15.  s_4        23. c_6h      31. t_d
            %     8.   d_2h       16.  c_4h       24. d_6       32. o_h
            
            import am_lib.* am_dft.*
            
            if nargin<2; from_memory=true; end
            
            if ischar(pg_code)
                pg={'c_1' ,'s_2' ,'c_2' ,'c_1h','c_2h','d_2' ,'c_2v','d_2h', ...
                    'c_3' ,'s_6' ,'d_3' ,'c_3v','d_3d','c_4' ,'s_4' ,'c_4h', ...
                    'd_4' ,'c_4v','d_2d','d_4h','c_6' ,'c_3h','c_6h','d_6' , ...
                    'c_6v','d_3h','d_6h','t'   ,'t_h' ,'o'   ,'t_d' ,'o_h'};
                pg_code = find(string(pg)==pg_code); %#ok<STRCLQT>
            end

            if from_memory
                % generated with this code:
                %
                % clear;clc;import am_dft.*
                % for i = 1:32
                %     R{i} = generate_pg(i,false);
                %     x(i) = identify_pointgroup(R{i});
                %     d{i} = decode_pg(x(i));
                %     fprintf('case %i; R = [',i); fprintf('%s',sym(R{i}(1))); fprintf(',%s',sym(R{i}(2:end))); fprintf(']; %% %s\n',d{i});
                % end
                % for i = 1:32
                %     [~,W{i}] = generate_pg(i,false);
                %     fprintf('case %i; W = [',i); fprintf('%s',sym(W{i}(1))); fprintf(',%s',sym(W{i}(2:end))); fprintf('];\n');
                % end
                %
                
                % single-valued representation
                switch pg_code
                case 1;  R = [1,0,0,0,1,0,0,0,1]; % c_1
                case 2;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1]; % s_2
                case 3;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1]; % c_2
                case 4;  R = [1,0,0,0,1,0,0,0,1,1,0,0,0,-1,0,0,0,1]; % c_1h
                case 5;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,-1,0,0,0,-1,0,0,0,-1,1,0,0,0,-1,0,0,0,1]; % c_2h
                case 6;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,1,0,0,0,-1,0,0,0,-1]; % d_2
                case 7;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,1,0,0,0,-1,0,0,0,1,-1,0,0,0,1,0,0,0,1]; % c_2v
                case 8;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,-1,0,0,0,-1,0,0,0,-1,1,0,0,0,-1,0,0,0,-1,1,0,0,0,1,0,0,0,-1,1,0,0,0,-1,0,0,0,1,-1,0,0,0,1,0,0,0,1]; % d_2h
                case 9;  R = [1,0,0,0,1,0,0,0,1,0,0,1,1,0,0,0,1,0,0,1,0,0,0,1,1,0,0]; % c_3
                case 10; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,0,0,1,1,0,0,0,1,0,0,0,-1,-1,0,0,0,-1,0,0,1,0,0,0,1,1,0,0,0,-1,0,0,0,-1,-1,0,0]; % s_6
                case 11; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,0,-1,0,-1,0,0,0,1,1,0,0,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,-1,0,-1,0,-1,0,0,0,1,0,0,0,1,1,0,0]; % d_3
                case 12; R = [1,0,0,0,1,0,0,0,1,1,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1,0,1,0,1,0,0,0,1,0,0,0,1,1,0,0]; % c_3v
                case 13; R = [1,0,0,0,1,0,0,0,1,1/3,-2/3,-2/3,-2/3,-2/3,1/3,-2/3,1/3,-2/3,-1,0,0,0,-1,0,0,0,-1,0,0,1,1,0,0,0,1,0,-1/3,2/3,2/3,2/3,2/3,-1/3,2/3,-1/3,2/3,-2/3,1/3,-2/3,1/3,-2/3,-2/3,-2/3,-2/3,1/3,0,0,-1,-1,0,0,0,-1,0,-2/3,-2/3,1/3,-2/3,1/3,-2/3,1/3,-2/3,-2/3,0,1,0,0,0,1,1,0,0,2/3,-1/3,2/3,-1/3,2/3,2/3,2/3,2/3,-1/3,2/3,2/3,-1/3,2/3,-1/3,2/3,-1/3,2/3,2/3,0,-1,0,0,0,-1,-1,0,0]; % d_3d
                case 14; R = [1,0,0,0,1,0,0,0,1,0,1,0,-1,0,0,0,0,1,-1,0,0,0,-1,0,0,0,1,0,-1,0,1,0,0,0,0,1]; % c_4
                case 15; R = [1,0,0,0,1,0,0,0,1,0,-1,0,1,0,0,0,0,-1,-1,0,0,0,-1,0,0,0,1,0,1,0,-1,0,0,0,0,-1]; % s_4
                case 16; R = [1,0,0,0,1,0,0,0,1,0,1,0,-1,0,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,-1,0,0,0,-1,0,0,0,1,0,-1,0,1,0,0,0,0,-1,0,-1,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,-1,0,1,0,-1,0,0,0,0,-1]; % c_4h
                case 17; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,1,0,-1,0,0,0,0,1,0,1,0,1,0,0,0,0,-1,0,-1,0,-1,0,0,0,0,-1,-1,0,0,0,-1,0,0,0,1,0,-1,0,1,0,0,0,0,1,1,0,0,0,-1,0,0,0,-1]; % d_4
                case 18; R = [1,0,0,0,1,0,0,0,1,0,1,0,-1,0,0,0,0,1,1,0,0,0,-1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,0,1,0,1,0,0,0,0,1,0,-1,0,-1,0,0,0,0,1,0,-1,0,1,0,0,0,0,1,-1,0,0,0,1,0,0,0,1]; % c_4v
                case 19; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,-1,0,1,0,0,0,0,-1,0,-1,0,-1,0,0,0,0,1,0,1,0,1,0,0,0,0,1,-1,0,0,0,-1,0,0,0,1,0,1,0,-1,0,0,0,0,-1,1,0,0,0,-1,0,0,0,-1]; % d_2d
                case 20; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,1,0,-1,0,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,0,1,0,1,0,0,0,0,-1,1,0,0,0,-1,0,0,0,1,0,-1,0,-1,0,0,0,0,-1,-1,0,0,0,-1,0,0,0,1,0,-1,0,1,0,0,0,0,-1,0,-1,0,1,0,0,0,0,1,1,0,0,0,-1,0,0,0,-1,0,-1,0,-1,0,0,0,0,1,0,1,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,-1,0,1,0,-1,0,0,0,0,-1,-1,0,0,0,1,0,0,0,1]; % d_4h
                case 21; R = [1,0,0,0,1,0,0,0,1,-1/3,2/3,2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,0,0,1,1,0,0,0,1,0,2/3,2/3,-1/3,-1/3,2/3,2/3,2/3,-1/3,2/3,0,1,0,0,0,1,1,0,0,2/3,-1/3,2/3,2/3,2/3,-1/3,-1/3,2/3,2/3]; % c_6
                case 22; R = [1,0,0,0,1,0,0,0,1,1/3,-2/3,-2/3,-2/3,1/3,-2/3,-2/3,-2/3,1/3,0,0,1,1,0,0,0,1,0,-2/3,-2/3,1/3,1/3,-2/3,-2/3,-2/3,1/3,-2/3,0,1,0,0,0,1,1,0,0,-2/3,1/3,-2/3,-2/3,-2/3,1/3,1/3,-2/3,-2/3]; % c_3h
                case 23; R = [1,0,0,0,1,0,0,0,1,-1/3,2/3,2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,-1,0,0,0,-1,0,0,0,-1,0,0,1,1,0,0,0,1,0,1/3,-2/3,-2/3,-2/3,1/3,-2/3,-2/3,-2/3,1/3,2/3,2/3,-1/3,-1/3,2/3,2/3,2/3,-1/3,2/3,0,0,-1,-1,0,0,0,-1,0,0,1,0,0,0,1,1,0,0,-2/3,-2/3,1/3,1/3,-2/3,-2/3,-2/3,1/3,-2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,-1/3,2/3,2/3,0,-1,0,0,0,-1,-1,0,0,-2/3,1/3,-2/3,-2/3,-2/3,1/3,1/3,-2/3,-2/3]; % c_6h
                case 24; R = [1,0,0,0,1,0,0,0,1,-1/3,2/3,2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,-1,0,0,0,0,-1,0,-1,0,0,0,1,1,0,0,0,1,0,1/3,-2/3,-2/3,-2/3,-2/3,1/3,-2/3,1/3,-2/3,2/3,2/3,-1/3,-1/3,2/3,2/3,2/3,-1/3,2/3,0,-1,0,-1,0,0,0,0,-1,0,0,-1,0,-1,0,-1,0,0,0,1,0,0,0,1,1,0,0,-2/3,1/3,-2/3,1/3,-2/3,-2/3,-2/3,-2/3,1/3,-2/3,-2/3,1/3,-2/3,1/3,-2/3,1/3,-2/3,-2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,-1/3,2/3,2/3]; % d_6
                case 25; R = [1,0,0,0,1,0,0,0,1,-1/3,2/3,2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,1,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,1,0,-1/3,2/3,2/3,2/3,2/3,-1/3,2/3,-1/3,2/3,2/3,2/3,-1/3,-1/3,2/3,2/3,2/3,-1/3,2/3,0,1,0,1,0,0,0,0,1,0,0,1,0,1,0,1,0,0,0,1,0,0,0,1,1,0,0,2/3,-1/3,2/3,-1/3,2/3,2/3,2/3,2/3,-1/3,2/3,2/3,-1/3,2/3,-1/3,2/3,-1/3,2/3,2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,-1/3,2/3,2/3]; % c_6v
                case 26; R = [1,0,0,0,1,0,0,0,1,1/3,-2/3,-2/3,-2/3,1/3,-2/3,-2/3,-2/3,1/3,1,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,1,0,1/3,-2/3,-2/3,-2/3,-2/3,1/3,-2/3,1/3,-2/3,-2/3,-2/3,1/3,1/3,-2/3,-2/3,-2/3,1/3,-2/3,0,1,0,1,0,0,0,0,1,0,0,1,0,1,0,1,0,0,0,1,0,0,0,1,1,0,0,-2/3,1/3,-2/3,1/3,-2/3,-2/3,-2/3,-2/3,1/3,-2/3,-2/3,1/3,-2/3,1/3,-2/3,1/3,-2/3,-2/3,-2/3,1/3,-2/3,-2/3,-2/3,1/3,1/3,-2/3,-2/3]; % d_3h
                case 27; R = [1,0,0,0,1,0,0,0,1,-1/3,2/3,2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,-1,0,0,0,0,-1,0,-1,0,0,0,1,1,0,0,0,1,0,-1,0,0,0,-1,0,0,0,-1,1/3,-2/3,-2/3,-2/3,-2/3,1/3,-2/3,1/3,-2/3,2/3,2/3,-1/3,-1/3,2/3,2/3,2/3,-1/3,2/3,1/3,-2/3,-2/3,-2/3,1/3,-2/3,-2/3,-2/3,1/3,0,-1,0,-1,0,0,0,0,-1,1,0,0,0,0,1,0,1,0,0,0,-1,0,-1,0,-1,0,0,0,1,0,0,0,1,1,0,0,0,0,-1,-1,0,0,0,-1,0,-2/3,1/3,-2/3,1/3,-2/3,-2/3,-2/3,-2/3,1/3,-1/3,2/3,2/3,2/3,2/3,-1/3,2/3,-1/3,2/3,-2/3,-2/3,1/3,-2/3,1/3,-2/3,1/3,-2/3,-2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,-1/3,2/3,2/3,-2/3,-2/3,1/3,1/3,-2/3,-2/3,-2/3,1/3,-2/3,0,1,0,1,0,0,0,0,1,0,0,1,0,1,0,1,0,0,0,-1,0,0,0,-1,-1,0,0,2/3,-1/3,2/3,-1/3,2/3,2/3,2/3,2/3,-1/3,2/3,2/3,-1/3,2/3,-1/3,2/3,-1/3,2/3,2/3,-2/3,1/3,-2/3,-2/3,-2/3,1/3,1/3,-2/3,-2/3]; % d_6h
                case 28; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,1,0,0,0,1,1,0,0,0,1,0,0,0,-1,-1,0,0,0,-1,0,0,0,1,-1,0,0,0,0,1,1,0,0,0,1,0,0,-1,0,0,0,-1,1,0,0,0,0,-1,-1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,1,0,0,0,-1,0,-1,0,0,0,-1,0,0,0,1,1,0,0,0,-1,0,0,0,-1]; % t
                case 29; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,1,0,0,0,1,1,0,0,-1,0,0,0,-1,0,0,0,-1,0,1,0,0,0,-1,-1,0,0,1,0,0,0,-1,0,0,0,1,0,-1,0,0,0,1,-1,0,0,0,0,1,1,0,0,0,1,0,0,-1,0,0,0,-1,-1,0,0,0,-1,0,0,0,-1,1,0,0,0,0,-1,-1,0,0,0,1,0,0,-1,0,0,0,1,1,0,0,0,0,1,-1,0,0,0,-1,0,0,1,0,0,0,-1,1,0,0,0,0,-1,1,0,0,0,-1,0,0,0,-1,-1,0,0,0,-1,0,0,1,0,0,0,1,-1,0,0,0,0,1,1,0,0,0,-1,0,-1,0,0,0,-1,0,0,0,1,0,0,-1,1,0,0,0,1,0,1,0,0,0,-1,0,0,0,-1,0,0,1,-1,0,0,0,1,0,1,0,0,0,1,0,0,0,-1,-1,0,0,0,1,0,0,0,1]; % t_h
                case 30; R = [1,0,0,0,1,0,0,0,1,0,1,0,0,0,1,1,0,0,0,1,0,-1,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,0,1,0,-1,0,1,0,0,-1,0,0,0,0,1,0,1,0,-1,0,0,0,-1,0,0,0,1,1,0,0,0,0,-1,0,1,0,0,-1,0,1,0,0,0,0,1,0,-1,0,0,0,-1,1,0,0,0,0,1,0,1,0,-1,0,0,0,-1,0,0,0,1,-1,0,0,0,0,-1,0,1,0,1,0,0,0,0,-1,-1,0,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,-1,1,0,0,0,-1,0,0,1,0,0,0,-1,-1,0,0,0,0,1,-1,0,0,0,-1,0,0,1,0,1,0,0,0,0,-1,-1,0,0,0,1,0,0,0,-1,1,0,0,0,-1,0,0,0,-1,0,0,-1,0,-1,0,-1,0,0,-1,0,0,0,0,-1,0,-1,0,0,-1,0,-1,0,0,0,0,-1]; % o
                case 31; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1,1,0,0,1,0,0,0,-1,0,0,0,-1,0,-1,0,-1,0,0,0,0,1,0,-1,0,0,0,1,-1,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,-1,-1,0,0,0,-1,0,1,0,0,0,0,-1,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,-1,1,0,0,0,0,1,0,1,0,1,0,0,0,0,1,1,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1,0,-1,0,-1,0,0,0,0,1,-1,0,0,0,-1,0,-1,0,0,0,0,-1,0,1,0,0,0,-1,0,1,0,-1,0,0,0,0,-1,-1,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,-1,0,-1,0,1,0,0,0,0,-1,1,0,0,0,-1,0]; % t_d
                case 32; R = [1,0,0,0,1,0,0,0,1,0,1,0,0,0,1,1,0,0,0,1,0,-1,0,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,0,0,1,1,0,0,0,1,0,0,0,1,0,-1,0,1,0,0,0,-1,0,0,0,-1,-1,0,0,-1,0,0,0,0,1,0,1,0,-1,0,0,0,-1,0,0,0,1,0,-1,0,1,0,0,0,0,-1,1,0,0,0,0,-1,0,1,0,0,0,-1,-1,0,0,0,-1,0,0,-1,0,1,0,0,0,0,1,0,-1,0,0,0,-1,1,0,0,0,0,-1,0,1,0,-1,0,0,0,0,1,0,1,0,-1,0,0,1,0,0,0,0,-1,0,-1,0,0,-1,0,0,0,1,-1,0,0,1,0,0,0,1,0,0,0,-1,0,0,-1,0,1,0,1,0,0,0,0,-1,-1,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,1,0,0,0,0,1,0,-1,0,0,1,0,-1,0,0,0,0,-1,0,0,-1,1,0,0,0,-1,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,-1,0,0,0,0,-1,0,-1,0,1,0,0,0,0,1,-1,0,0,0,-1,0,0,1,0,0,0,-1,1,0,0,0,1,0,1,0,0,0,0,-1,0,0,1,0,-1,0,-1,0,0,-1,0,0,0,1,0,0,0,-1,0,0,1,1,0,0,0,-1,0,-1,0,0,0,0,-1,0,1,0,1,0,0,0,-1,0,0,0,-1,0,0,1,-1,0,0,0,1,0,0,0,-1,0,-1,0,-1,0,0,0,-1,0,0,0,1,1,0,0,-1,0,0,0,0,-1,0,-1,0,0,0,-1,1,0,0,0,1,0,0,-1,0,-1,0,0,0,0,1,1,0,0,0,-1,0,0,0,1,0,-1,0,-1,0,0,0,0,-1,-1,0,0,0,1,0,0,0,1,0,0,1,0,1,0,1,0,0,1,0,0,0,0,1,0,1,0,0,1,0,1,0,0,0,0,1]; % o_h
                end
                R = reshape(R,3,3,[]);
                
                % double-valued representation
                switch pg_code
                case 1; W = [1,0,0,1,-1,0,0,-1];
                case 2; W = [1,0,0,1,1i,0,0,1i,-1,0,0,-1,-1i,0,0,-1i];
                case 3; W = [1,0,0,1,0,-1,1,0,-1,0,0,-1,0,1,-1,0];
                case 4; W = [1,0,0,1,0,-1i,1i,0,-1,0,0,-1,0,1i,-1i,0];
                case 5; W = [1,0,0,1,0,-1,1,0,1i,0,0,1i,0,-1i,1i,0,-1,0,0,-1,0,1,-1,0,-1i,0,0,-1i,0,1i,-1i,0];
                case 6; W = [1,0,0,1,1i,0,0,-1i,0,-1,1,0,0,1i,1i,0,-1,0,0,-1,-1i,0,0,1i,0,1,-1,0,0,-1i,-1i,0];
                case 7; W = [1,0,0,1,1i,0,0,-1i,0,-1i,1i,0,0,-1,-1,0,-1,0,0,-1,-1i,0,0,1i,0,1i,-1i,0,0,1,1,0];
                case 8; W = [1,0,0,1,1i,0,0,-1i,0,-1,1,0,1i,0,0,1i,0,1i,1i,0,-1,0,0,1,0,-1i,1i,0,0,-1,-1,0,-1,0,0,-1,-1i,0,0,1i,0,1,-1,0,-1i,0,0,-1i,0,-1i,-1i,0,1,0,0,-1,0,1i,-1i,0,0,1,1,0];
                case 9; W = [1,0,0,1,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,-1,0,0,-1,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2];
                case 10; W = [1,0,0,1,1i,0,0,1i,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,-1,0,0,-1,-1i,0,0,-1i,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2];
                case 11; W = [1,0,0,1,(2^(1/2)*1i)/2,2^(1/2)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(1/2 - 1i/2),0,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,-1,0,0,-1,-(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)/2,(2^(1/2)*1i)/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2];
                case 12; W = [1,0,0,1,-2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,-1,0,0,-1,2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2];
                case 13; W = [1,0,0,1,(6^(1/2)*1i)/6,- (2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,6^(1/2)/6 - (2^(1/2)*3^(1/2)*1i)/3,-(6^(1/2)*1i)/6,1i,0,0,1i,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,-6^(1/2)/6,(2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 + (6^(1/2)*1i)/6,6^(1/2)/6,-(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)*(- 1/6 + 1i/6),6^(1/2)*(1/6 + 1i/6),(2^(1/2)*3^(1/2)*1i)/3,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,-(6^(1/2)*1i)/6,- (2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(6^(1/2)*1i)/6,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,(2^(1/2)*3^(1/2))/3,6^(1/2)*(- 1/6 - 1i/6),6^(1/2)*(- 1/6 + 1i/6),-(2^(1/2)*3^(1/2))/3,6^(1/2)/6,6^(1/2)/6 - (2^(1/2)*3^(1/2)*1i)/3,(2^(1/2)*3^(1/2)*1i)/3 + 6^(1/2)/6,-6^(1/2)/6,1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,-1,0,0,-1,-(6^(1/2)*1i)/6,(2^(1/2)*3^(1/2)*1i)/3 + 6^(1/2)/6,(2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,(6^(1/2)*1i)/6,-1i,0,0,-1i,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,6^(1/2)/6,(6^(1/2)*1i)/6 - (2^(1/2)*3^(1/2))/3,- (2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,-6^(1/2)/6,(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)*(1/6 - 1i/6),6^(1/2)*(- 1/6 - 1i/6),-(2^(1/2)*3^(1/2)*1i)/3,1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,(6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 + (6^(1/2)*1i)/6,(6^(1/2)*1i)/6 - (2^(1/2)*3^(1/2))/3,-(6^(1/2)*1i)/6,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,-(2^(1/2)*3^(1/2))/3,6^(1/2)*(1/6 + 1i/6),6^(1/2)*(1/6 - 1i/6),(2^(1/2)*3^(1/2))/3,-6^(1/2)/6,(2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,- (2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,6^(1/2)/6,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2];
                case 14; W = [1,0,0,1,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),1i,0,0,-1i,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),-1,0,0,-1,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),-1i,0,0,1i,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2)];
                case 15; W = [1,0,0,1,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),1i,0,0,-1i,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),-1,0,0,-1,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),-1i,0,0,1i,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2)];
                case 16; W = [1,0,0,1,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),1i,0,0,1i,1i,0,0,-1i,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),-1,0,0,1,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),-1,0,0,-1,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),-1i,0,0,-1i,-1i,0,0,1i,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),1,0,0,-1,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2)];
                case 17; W = [1,0,0,1,0,-1,1,0,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(1/2 - 1i/2),0,1i,0,0,-1i,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),0,1i,1i,0,-1,0,0,-1,0,1,-1,0,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,-1i,0,0,1i,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),0,-1i,-1i,0];
                case 18; W = [1,0,0,1,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),0,-1i,1i,0,1i,0,0,-1i,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),0,-1,-1,0,-1,0,0,-1,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),0,1i,-1i,0,-1i,0,0,1i,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(1/2 - 1i/2),0,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),0,1,1,0];
                case 19; W = [1,0,0,1,0,-1,1,0,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,1i,0,0,-1i,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),0,1i,1i,0,-1,0,0,-1,0,1,-1,0,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,-1i,0,0,1i,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),0,-1i,-1i,0];
                case 20; W = [1,0,0,1,0,-1,1,0,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),1i,0,0,1i,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(1/2 + 1i/2),0,0,-1i,1i,0,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(1/2 - 1i/2),0,1i,0,0,-1i,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),0,1i,1i,0,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,-1,0,0,1,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),0,-1,-1,0,-1,0,0,-1,0,1,-1,0,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),-1i,0,0,-1i,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,0,1i,-1i,0,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,-1i,0,0,1i,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),0,-1i,-1i,0,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,1,0,0,-1,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),0,1,1,0];
                case 21; W = [1,0,0,1,-(3^(1/2)*1i)/3,3^(1/2)*(1/3 - 1i/3),3^(1/2)*(- 1/3 - 1i/3),(3^(1/2)*1i)/3,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,3^(1/2)*(1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(1/2 + 1i/6),1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,3^(1/2)*(1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/2 - 1i/6),-1,0,0,-1,(3^(1/2)*1i)/3,3^(1/2)*(- 1/3 + 1i/3),3^(1/2)*(1/3 + 1i/3),-(3^(1/2)*1i)/3,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,3^(1/2)*(- 1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(- 1/2 - 1i/6),- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,3^(1/2)*(- 1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/2 + 1i/6)];
                case 22; W = [1,0,0,1,3^(1/2)/3,3^(1/2)*(1/3 + 1i/3),3^(1/2)*(1/3 - 1i/3),-3^(1/2)/3,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,3^(1/2)*(1/6 + 1i/2),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/2),1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,3^(1/2)*(- 1/6 + 1i/2),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/2),-1,0,0,-1,-3^(1/2)/3,3^(1/2)*(- 1/3 - 1i/3),3^(1/2)*(- 1/3 + 1i/3),3^(1/2)/3,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,3^(1/2)*(- 1/6 - 1i/2),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 - 1i/2),- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,3^(1/2)*(1/6 - 1i/2),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/2)];
                case 23; W = [1,0,0,1,-(3^(1/2)*1i)/3,3^(1/2)*(1/3 - 1i/3),3^(1/2)*(- 1/3 - 1i/3),(3^(1/2)*1i)/3,1i,0,0,1i,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,3^(1/2)/3,3^(1/2)*(1/3 + 1i/3),3^(1/2)*(1/3 - 1i/3),-3^(1/2)/3,3^(1/2)*(1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(1/2 + 1i/6),- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,3^(1/2)*(1/6 + 1i/2),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/2),3^(1/2)*(1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/2 - 1i/6),1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,3^(1/2)*(- 1/6 + 1i/2),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/2),-1,0,0,-1,(3^(1/2)*1i)/3,3^(1/2)*(- 1/3 + 1i/3),3^(1/2)*(1/3 + 1i/3),-(3^(1/2)*1i)/3,-1i,0,0,-1i,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,-3^(1/2)/3,3^(1/2)*(- 1/3 - 1i/3),3^(1/2)*(- 1/3 + 1i/3),3^(1/2)/3,3^(1/2)*(- 1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(- 1/2 - 1i/6),1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,3^(1/2)*(- 1/6 - 1i/2),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 - 1i/2),3^(1/2)*(- 1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/2 + 1i/6),- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,3^(1/2)*(1/6 - 1i/2),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/2)];
                case 24; W = [1,0,0,1,-(3^(1/2)*1i)/3,3^(1/2)*(1/3 - 1i/3),3^(1/2)*(- 1/3 - 1i/3),(3^(1/2)*1i)/3,(2^(1/2)*1i)/2,2^(1/2)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,(6^(1/2)*1i)/6,- (2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,6^(1/2)/6 - (2^(1/2)*3^(1/2)*1i)/3,-(6^(1/2)*1i)/6,3^(1/2)*(1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(1/2 + 1i/6),0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(1/2 - 1i/2),0,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,-(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)*(- 1/6 + 1i/6),6^(1/2)*(1/6 + 1i/6),(2^(1/2)*3^(1/2)*1i)/3,-(6^(1/2)*1i)/6,- (2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(6^(1/2)*1i)/6,3^(1/2)*(1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/2 - 1i/6),-1,0,0,-1,(3^(1/2)*1i)/3,3^(1/2)*(- 1/3 + 1i/3),3^(1/2)*(1/3 + 1i/3),-(3^(1/2)*1i)/3,-(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)/2,(2^(1/2)*1i)/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,-(6^(1/2)*1i)/6,(2^(1/2)*3^(1/2)*1i)/3 + 6^(1/2)/6,(2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,(6^(1/2)*1i)/6,3^(1/2)*(- 1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(- 1/2 - 1i/6),0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)*(1/6 - 1i/6),6^(1/2)*(- 1/6 - 1i/6),-(2^(1/2)*3^(1/2)*1i)/3,(6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 + (6^(1/2)*1i)/6,(6^(1/2)*1i)/6 - (2^(1/2)*3^(1/2))/3,-(6^(1/2)*1i)/6,3^(1/2)*(- 1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/2 + 1i/6)];
                case 25; W = [1,0,0,1,-(3^(1/2)*1i)/3,3^(1/2)*(1/3 - 1i/3),3^(1/2)*(- 1/3 - 1i/3),(3^(1/2)*1i)/3,-2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,-6^(1/2)/6,(2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 + (6^(1/2)*1i)/6,6^(1/2)/6,3^(1/2)*(1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(1/2 + 1i/6),0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,(2^(1/2)*3^(1/2))/3,6^(1/2)*(- 1/6 - 1i/6),6^(1/2)*(- 1/6 + 1i/6),-(2^(1/2)*3^(1/2))/3,6^(1/2)/6,6^(1/2)/6 - (2^(1/2)*3^(1/2)*1i)/3,(2^(1/2)*3^(1/2)*1i)/3 + 6^(1/2)/6,-6^(1/2)/6,3^(1/2)*(1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/2 - 1i/6),-1,0,0,-1,(3^(1/2)*1i)/3,3^(1/2)*(- 1/3 + 1i/3),3^(1/2)*(1/3 + 1i/3),-(3^(1/2)*1i)/3,2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,6^(1/2)/6,(6^(1/2)*1i)/6 - (2^(1/2)*3^(1/2))/3,- (2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,-6^(1/2)/6,3^(1/2)*(- 1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(- 1/2 - 1i/6),0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,-(2^(1/2)*3^(1/2))/3,6^(1/2)*(1/6 + 1i/6),6^(1/2)*(1/6 - 1i/6),(2^(1/2)*3^(1/2))/3,-6^(1/2)/6,(2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,- (2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,6^(1/2)/6,3^(1/2)*(- 1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/2 + 1i/6)];
                case 26; W = [1,0,0,1,3^(1/2)/3,3^(1/2)*(1/3 + 1i/3),3^(1/2)*(1/3 - 1i/3),-3^(1/2)/3,-2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,(6^(1/2)*1i)/6,- (2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,6^(1/2)/6 - (2^(1/2)*3^(1/2)*1i)/3,-(6^(1/2)*1i)/6,3^(1/2)*(1/6 + 1i/2),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/2),0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,-(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)*(- 1/6 + 1i/6),6^(1/2)*(1/6 + 1i/6),(2^(1/2)*3^(1/2)*1i)/3,-(6^(1/2)*1i)/6,- (2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(6^(1/2)*1i)/6,3^(1/2)*(- 1/6 + 1i/2),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/2),-1,0,0,-1,-3^(1/2)/3,3^(1/2)*(- 1/3 - 1i/3),3^(1/2)*(- 1/3 + 1i/3),3^(1/2)/3,2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,-(6^(1/2)*1i)/6,(2^(1/2)*3^(1/2)*1i)/3 + 6^(1/2)/6,(2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,(6^(1/2)*1i)/6,3^(1/2)*(- 1/6 - 1i/2),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 - 1i/2),0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)*(1/6 - 1i/6),6^(1/2)*(- 1/6 - 1i/6),-(2^(1/2)*3^(1/2)*1i)/3,(6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 + (6^(1/2)*1i)/6,(6^(1/2)*1i)/6 - (2^(1/2)*3^(1/2))/3,-(6^(1/2)*1i)/6,3^(1/2)*(1/6 - 1i/2),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/2)];
                case 27; W = [1,0,0,1,-(3^(1/2)*1i)/3,3^(1/2)*(1/3 - 1i/3),3^(1/2)*(- 1/3 - 1i/3),(3^(1/2)*1i)/3,(2^(1/2)*1i)/2,2^(1/2)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,1i,0,0,1i,(6^(1/2)*1i)/6,- (2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,6^(1/2)/6 - (2^(1/2)*3^(1/2)*1i)/3,-(6^(1/2)*1i)/6,3^(1/2)*(1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(1/2 + 1i/6),3^(1/2)/3,3^(1/2)*(1/3 + 1i/3),3^(1/2)*(1/3 - 1i/3),-3^(1/2)/3,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(1/2 - 1i/2),0,-2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,-(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)*(- 1/6 + 1i/6),6^(1/2)*(1/6 + 1i/6),(2^(1/2)*3^(1/2)*1i)/3,-6^(1/2)/6,(2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 + (6^(1/2)*1i)/6,6^(1/2)/6,-(6^(1/2)*1i)/6,- (2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(6^(1/2)*1i)/6,3^(1/2)*(1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/2 - 1i/6),3^(1/2)*(1/6 + 1i/2),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/2),0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,(2^(1/2)*3^(1/2))/3,6^(1/2)*(- 1/6 - 1i/6),6^(1/2)*(- 1/6 + 1i/6),-(2^(1/2)*3^(1/2))/3,6^(1/2)/6,6^(1/2)/6 - (2^(1/2)*3^(1/2)*1i)/3,(2^(1/2)*3^(1/2)*1i)/3 + 6^(1/2)/6,-6^(1/2)/6,3^(1/2)*(- 1/6 + 1i/2),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/2),-1,0,0,-1,(3^(1/2)*1i)/3,3^(1/2)*(- 1/3 + 1i/3),3^(1/2)*(1/3 + 1i/3),-(3^(1/2)*1i)/3,-(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)/2,(2^(1/2)*1i)/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,-1i,0,0,-1i,-(6^(1/2)*1i)/6,(2^(1/2)*3^(1/2)*1i)/3 + 6^(1/2)/6,(2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,(6^(1/2)*1i)/6,3^(1/2)*(- 1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(- 1/2 - 1i/6),-3^(1/2)/3,3^(1/2)*(- 1/3 - 1i/3),3^(1/2)*(- 1/3 + 1i/3),3^(1/2)/3,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)*(1/6 - 1i/6),6^(1/2)*(- 1/6 - 1i/6),-(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)/6,(6^(1/2)*1i)/6 - (2^(1/2)*3^(1/2))/3,- (2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,-6^(1/2)/6,(6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 + (6^(1/2)*1i)/6,(6^(1/2)*1i)/6 - (2^(1/2)*3^(1/2))/3,-(6^(1/2)*1i)/6,3^(1/2)*(- 1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/2 + 1i/6),3^(1/2)*(- 1/6 - 1i/2),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 - 1i/2),0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,-(2^(1/2)*3^(1/2))/3,6^(1/2)*(1/6 + 1i/6),6^(1/2)*(1/6 - 1i/6),(2^(1/2)*3^(1/2))/3,-6^(1/2)/6,(2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,- (2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,6^(1/2)/6,3^(1/2)*(1/6 - 1i/2),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/2)];
                case 28; W = [1,0,0,1,0,-1,1,0,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1i,0,0,-1i,0,1i,1i,0,-1,0,0,-1,0,1,-1,0,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,-1i,0,0,1i,0,-1i,-1i,0];
                case 29; W = [1,0,0,1,0,-1,1,0,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,1i,0,0,1i,1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,0,-1i,1i,0,1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1i,0,0,-1i,1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,0,1i,1i,0,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,-1,0,0,1,0,-1,-1,0,-1,0,0,-1,0,1,-1,0,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,-1i,0,0,-1i,- 1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,0,1i,-1i,0,- 1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,-1i,0,0,1i,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,0,-1i,-1i,0,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1,0,0,-1,0,1,1,0];
                case 30; W = [1,0,0,1,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)/2,-(2^(1/2)*1i)/2,1i,0,0,-1i,2^(1/2)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,2^(1/2)/2,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 - 1i/2,2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,2^(1/2)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(1/2 + 1i/2),0,0,-1,1,0,0,1i,1i,0,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,2^(1/2)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(1/2 - 1i/2),0,-1,0,0,-1,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,-2^(1/2)/2,(2^(1/2)*1i)/2,-1i,0,0,1i,-2^(1/2)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),- 1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,-2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,- 1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,-2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,-2^(1/2)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,0,1,-1,0,0,-1i,-1i,0,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)/2,(2^(1/2)*1i)/2,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(- 1/2 + 1i/2),0];
                case 31; W = [1,0,0,1,1i,0,0,-1i,0,-1,1,0,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,0,1i,1i,0,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 - 1i/2,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),-2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,(2^(1/2)*1i)/2,-2^(1/2)/2,-2^(1/2)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,(2^(1/2)*1i)/2,2^(1/2)/2,2^(1/2)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,-1,0,0,-1,-1i,0,0,1i,0,1,-1,0,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,0,-1i,-1i,0,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(1/2 - 1i/2),0,- 1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),- 1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,- 1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,-(2^(1/2)*1i)/2,2^(1/2)/2,2^(1/2)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2];
                case 32; W = [1,0,0,1,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),1i,0,0,1i,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)/2,-(2^(1/2)*1i)/2,1i,0,0,-1i,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,2^(1/2)/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,-2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,2^(1/2)/2,1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 - 1i/2,-1,0,0,1,2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,(2^(1/2)*1i)/2,-2^(1/2)/2,-2^(1/2)/2,(2^(1/2)*1i)/2,2^(1/2)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 + 1i/2,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(1/2 + 1i/2),0,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,0,-1,1,0,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,(2^(1/2)*1i)/2,2^(1/2)/2,2^(1/2)/2,(2^(1/2)*1i)/2,0,1i,1i,0,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,(2^(1/2)*1i)/2,2^(1/2)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,0,-1i,1i,0,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(1/2 - 1i/2),0,0,-1,-1,0,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,-1,0,0,-1,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),-1i,0,0,-1i,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,-(2^(1/2)*1i)/2,2^(1/2)/2,-2^(1/2)/2,(2^(1/2)*1i)/2,-1i,0,0,1i,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),-2^(1/2)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),- 1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,- 1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1,0,0,-1,-2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,-(2^(1/2)*1i)/2,2^(1/2)/2,2^(1/2)/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,0,1,-1,0,- 1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,0,-1i,-1i,0,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,- 1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 - 1i/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)/2,(2^(1/2)*1i)/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(1/2 - 1i/2),0,0,1i,-1i,0,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,0,1,1,0,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0];
                end
                W = reshape(W,2,2,[]);
                
            else
                % load recipe
                recipe = get_recipe(pg_code);
                % initialize
                R = zeros(3,3,numel(recipe)+1);
                % make generators
                nsyms = 0;
                nsyms = nsyms+1; R(:,:,nsyms) = eye(3);
                for k = 1:numel(recipe)
                    nsyms=nsyms+1; R(:,:,nsyms) = decode_R(recipe(k)); 
                end
                
                % add elements until no new elements are generated
                R = complete_group(R);
                
                % [IMPORTANT] use rhombohedral rather than hexagonal setting for hexagonal point groups
                if any(pg_code==[9:13,21:27])
                    T = [-1 1 1; 2 1 1; -1 -2 1].'/3;
                    R = matmul_(matmul_(inv(T),R),(T));
                end
                
                % generate double group
                if nargout > 1
                    % get double group
                    j=1/2; [W] = get_wigner(j,R,'spherical');                     
                    % Add plus and minus in accordance with 
                    % V. Heine, Group Theory in Quantum Mechanics (Elsevier, 2014), p 62, eq 8.24.
                    W = cat(3,W,-W); % W = complete_group(reshape(uniquec_(reshape(wdv_(W),4,[])),2,2,[]));
                    % remove numerical noise
                    W = wdv_(W);
                end
            end

            function recipe  = get_recipe(sg_id)
                pg_recipe_database = {...
                     '','h','c','j','ch','bc','bj','bch','n','hn','en','kn','fhn','g','m','gh','cg',...
                    'gj','cm','cgh','bn','in','bhn','ben','bkn','ikn','benh','cd','cdh','dg','bcld','dgh'};
                recipe = pg_recipe_database{sg_id};
            end
            function R       = decode_R(str_r)
                switch str_r
                case 'a'; R = [  1  0  0;  0  1  0;  0  0  1]; % a
                case 'b'; R = [ -1  0  0;  0 -1  0;  0  0  1]; % b
                case 'c'; R = [ -1  0  0;  0  1  0;  0  0 -1]; % c
                case 'd'; R = [  0  0  1;  1  0  0;  0  1  0]; % d
                case 'e'; R = [  0  1  0;  1  0  0;  0  0 -1]; % e
                case 'f'; R = [  0 -1  0; -1  0  0;  0  0 -1]; % f
                case 'g'; R = [  0 -1  0;  1  0  0;  0  0  1]; % g
                case 'h'; R = [ -1  0  0;  0 -1  0;  0  0 -1]; % h
                case 'i'; R = [  1  0  0;  0  1  0;  0  0 -1]; % i
                case 'j'; R = [  1  0  0;  0 -1  0;  0  0  1]; % j
                case 'k'; R = [  0 -1  0; -1  0  0;  0  0  1]; % k
                case 'l'; R = [  0  1  0;  1  0  0;  0  0  1]; % l
                case 'm'; R = [  0  1  0; -1  0  0;  0  0 -1]; % m
                case 'n'; R = [  0 -1  0;  1 -1  0;  0  0  1]; % n
                end
            end 
        end

        function [S,MT]       = complete_group(S,tol)

            import am_lib.*
            
            % set tolernece
            if nargin<2; tol=am_dft.tiny; end

            % get size
            s=size(S); d=s(1)*s(2);
            
            % switch
            if s(1)==4 && s(2)==4 && all_(eq_(S(4,1:4,:), [0,0,0,1], tol))
                algo = 2; % seitz symmetry
            else
                algo = 1; % point symmetry
            end

            % exclude empty symmetries
            S = S(:,:,any(reshape(S,d,[])~=0,1)); 
            
            % add identity if not present 
            E = eye(s(1),s(2)); if ~any(all(all(S==E))); S = cat(3,E,S); end
            
            % get number of symmetries at this point
            nsyms=size(S,3);

            % allocate space
            S=cat(3,S,zeros(size(S,1),size(S,1),2*192-nsyms));

            % add symmetry until no new symmetry is generated
            nsyms_last=0; MT = zeros(2*192,2*192); 
            while nsyms ~= nsyms_last
                nsyms_last=nsyms;
                for i = 1:nsyms_last
                for j = 1:nsyms_last
                    if MT(i,j)==0
                        A = symmul_(S(:,:,i),S(:,:,j),algo,tol);
                        A_id = member_(A(:),reshape(S(:,:,1:nsyms),d,[]),tol);
                        if A_id == 0
                            nsyms = nsyms+1; S(:,:,nsyms) = A; MT(i,j) = nsyms;
                        else
                            MT(i,j) = A_id;
                        end
                    end
                end
                end
            end
            
            % trim output
            S  = S(:,:,1:nsyms);
            MT = MT(1:nsyms,1:nsyms);
            
            function C = symmul_(A,B,algo,tol)
                switch algo
                    case 1
                        % point symmetry
                        C = A*B;
                    case 2
                        % seitz symmetry
                        C = A*B; C(1:3,4)=mod(C(1:3,4)+tol,1)-tol;
                end
            end
        end
 
        function abc          = bas2abc(bas) % convert from column basis vectors to [a,b,c,alpha,beta,gamma] 
            % [a,b,c,alpha,beta,gamma] = bas2abc(bas)
            M = bas.' * bas;
            % de Graef p 86
            a = sqrt(M(1,1));
            b = sqrt(M(2,2));
            c = sqrt(M(3,3));
            alpha= acosd(M(2,3)/(b*c));
            beta = acosd(M(1,3)/(a*c));
            gamma= acosd(M(1,2)/(a*b));
            abc = [a,b,c,alpha,beta,gamma];
        end

        function bas          = abc2bas(abc,brav) % convert from [a,b,c,alpha,beta,gamma] to column basis vectors 
            import am_cell.abc2bas
            switch nargin
                case 1
                    % bas = abc2bas([a,b,c,alpha,beta,gamma])
                    a=abc(1); alpha=abc(4);
                    b=abc(2); beta =abc(5);
                    c=abc(3); gamma=abc(6);
                    % correct rounding errors
                    abs_cap = @(x) max(min(x,180), -180);
                    val = abs_cap( (cosd(alpha)*cosd(beta)-cosd(gamma))/(sind(alpha)*sind(beta)) ); gamma_star = acosd(val);
                    % generate basis (column vectors)
                    bas(:,1) = [ a * sind(beta)                    ,                                0.0, a * cosd(beta) ];
                    bas(:,2) = [-b * sind(alpha) * cosd(gamma_star), b * sind(alpha) * sind(gamma_star), b * cosd(alpha)];
                    bas(:,3) = [                                0.0,                                0.0, c              ];

                case 2
                    if      contains(brav,'cubic')
                        % [a, a, a, 90, 90, 90]
                        bas = abc2bas([abc(1),abc(1),abc(1),90,90,90]);
                    elseif  contains(brav,'tetra')
                        % [a, a, c, 90, 90, 90]
                        bas = abc2bas([abc(1),abc(1),abc(2),90,90,90]);
                    elseif  contains(brav,'orth')
                        % [a, b, c, 90, 90, 90]
                        bas = abc2bas([abc(1),abc(2),abc(3),90,90,90]);
                    elseif  contains(brav,'mon')
                        % [a, b, c, 90, beta, 90]
                        bas = abc2bas([abc(1),abc(2),abc(3),90,abc(4),90]);
                    elseif  contains(brav,'hex')
                        % [a, a, c, 90, 90, 120]
                        bas = abc2bas([abc(1),abc(1),abc(2),90,90,120]);
                    elseif  contains(brav,'rhomb')
                        % [a, a, a, alpha, alpha, alpha]
                        bas = abc2bas([abc(1),abc(1),abc(1),abc(2),abc(2),abc(2)]);
                    else
                        error('abc2bas: unknown bravis lattice');
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

            import am_lib.*
            
            % stupid matlab requires these symbolic variables to be initialized
            x=[];y=[];z=[]; %#ok<NASGU>
            
            % do not get conventional cell. it is full if bugs
            do_cc = false; cc=[];
            
            % time
            fprintf(' ... getting cell'); tic
            
            % validate input
            validatestring(flag,{'poscar','cif','material'});
            switch flag
                case {'poscar','cif'};  if exist(arg,'file')~=2; fprintf('\n'); error('File does not exist: %s',arg); end
            end
            
            % convert bas from [Ang] to [nm] when loading from cif/poscar
            switch flag
                case 'poscar';   [uc]     = am_cell.load_poscar(arg); 
                case 'cif';      [uc,str] = load_cif(arg);      uc.bas = uc.bas*0.1;
                case 'material'; [uc]     = load_material(arg);
                case 'create';   [uc]     = am_dft.create_cell(arg{:});
            end
            
            % set tolerance
            if nargin>2; uc.tol = tol; end

            % get primitive cell
            [pc,p2u,u2p] = get_primitive(uc);

            % get irreducible cell
            [ic,i2p,p2i] = get_irreducible(pc);

            % get conventional cell
            if do_cc
            [cc,c2p,p2c] = get_conventional_cell(pc); 
            end

            % complete mapping
            u2i = p2i(u2p); i2u = p2u(i2p);
            if do_cc
            u2c = p2c(u2p); c2u = p2u(c2p);
            c2i = p2i(c2p); i2c = p2c(i2p); 
            end
            
            % augment properties
            for plist = {'u2p','p2u','u2i','i2u','bas2pc','tau2pc'}; addprop(uc,plist{:}); end
            for plist = {'p2u','u2p','p2i','i2p'}; addprop(pc,plist{:}); end
            for plist = {'i2p','p2i','i2u','u2i'}; addprop(ic,plist{:}); end

            % sace mapping to cells
            pc.p2i = p2i; pc.p2u = p2u; if do_cc; pc.p2c = p2c; end
            pc.i2p = i2p; pc.u2p = u2p; if do_cc; pc.c2p = c2p; end

            ic.i2p = i2p; ic.i2u = i2u; if do_cc; ic.i2c = i2c; end
            ic.p2i = p2i; ic.u2i = u2i; if do_cc; ic.u2i = c2i; end

            uc.u2p = u2p; uc.u2i = u2i; if do_cc; uc.u2c = u2c; end
            uc.p2u = p2u; uc.i2u = i2u; if do_cc; uc.c2u = c2u; end
           
%             if do_cc; cc.c2i = c2i; cc.c2u = c2u; cc.c2p = c2p; end
%             if do_cc; cc.i2c = i2c; cc.u2c = u2c; cc.p2c = p2c; end
%             
%             % save bas2pc and tau2pc to convert [uc/cc-frac] to [pc-frac]
%             uc.bas2pc = pc.bas/uc.bas; uc.tau2pc = pc.bas\uc.bas;
%             if do_cc
%             cc.bas2pc = pc.bas/cc.bas; cc.tau2pc = pc.bas\cc.bas;
%             end

            % print basic symmetry info
            [~,H,~,R] = get_symmetries(pc);
            bv_code = am_cell.identify_bravais(pc.bas, pc.tol);

            % holohodry should give same info as bravais lattice
            hg_code = am_cell.identify_pointgroup(H,pc.tol); 
            pg_code = am_cell.identify_pointgroup(R,pc.tol); 
            sg_code = am_cell.identify_spacegroup(pg_code,pc.tol); % BETA

            % print relevant information
            verbose = true;
            if verbose
                fprintf(' (%.3f s) \n',toc);
                fprintf('     %-16s = %s\n','formula',uc.formula);
                fprintf('     %-16s = %s\n','primitive',am_cell.decode_bravais(bv_code));
                fprintf('     %-16s = %s\n','holohodry',am_cell.decode_holohodry(hg_code));
                fprintf('     %-16s = %s\n','point group',am_cell.decode_pg(pg_code));
                fprintf('     %-16s = %s\n','laue group',am_cell.decode_laue(am_cell.identify_laue(pg_code)));
                fprintf('     %-16s = %s\n','space group',strrep(cell2mat(join(am_cell.decode_sg(sg_code),',')),' ',''));
                fprintf('     %-16s = %-8.3f [g/cm3] \n','mass density',uc.mass_density);
                fprintf('     %-16s = %-8.3f [atoms/nm3]\n','number density',uc.numb_density);
                fprintf('     %-16s = %-8.3f [f.u./nm3]\n','formula density',uc.form_density);
                fprintf('     %-16s = %-8.3f [amu/f.u.]\n','molecular weight',uc.mole_weight);
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
                    'symb',{symb},'mass',am_dft.get_atomic_mass(symb),'nspecies',sum(1:max(species(:))==species(:),1), ...
                    'natoms',numel(species),'tau',tau,'species',species);
                uc = uc_(bas,symb,species,tau);
                
                if numel(uc.mass) ~= numel(uc.symb); error('mass and symb mismatch'); end 
                if numel(uc.mass) ~= numel(uc.nspecies); error('mass and nspecies mismatch'); end 
                
                uc_gen = [];
                eval(['uc_gen = am_dft.',str]);
                if uc_gen.natoms ~= uc.natoms || ...
                   am_lib.any_(~am_lib.eq_(am_lib.sortc_(uc_gen.tau),am_lib.sortc_(uc.tau)))
                    str = 'Mismatch';
                end
            end

        end
        
        function [uc]         = load_poscar(fposcar) % can also be used to load the charge density 
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
                if contains(fposcar,'CHGCAR') || contains(fposcar,'CHG')
                    fgetl(fid);
                    n = sscanf(fgetl(fid),'%i');
                    t = textscan(fid,'%f');
                    uc.nchg = n(:).';
                    uc.chg  = reshape(t{:},uc.nchg);
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
