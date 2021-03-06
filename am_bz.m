classdef am_bz %  < (dynamicprops is a terrible idea)
    
    properties
        % basic properties
        units    = []; % frac/cart
        recbas   = []; % reciprocal basis 
        n        = []; % monkhorst-pack dimensions or zone path divisions
        nks      = []; % number of kpoints
        k        = []; % kpoint coordinates [3,nks]
        w        = []; % kpoint weights [1,nks]
        x        = []; % linear coordinate for plotting path
        xt       = []; % linear coordinate ticks
        xl       = []; % linear coordinate labels
        tol      = []; % numerical tolerance
        % dispersions
        E        = []; % electron dispersion [nbands,nks,nmds] % molecular dynamics times steps if available
        hw       = []; % phonon dispersion [nbands,nks]
        % projection information
        nbands   = []; % number of bands
        norbitals= []; % number of orbitals
        natoms   = []; % number of atoms
        nspins   = []; % number of spins
        lmproj   = []; % electronic projections [nspins,norbitals,natoms,nbands,nks]
        eigenv   = []; % phonon eigenvectors :: need to change this to: [1,3,natoms,nbands,nks]
        optmat   = []; % squared optical matrix elements [nspins,3,nbands,nbands,nks]
        % oribtal labels
        orbital  = []; % orbital labels
        % tetrahedron information
        ntets    = []; % number of tetrahedra 
        tet      = []; % tetrahedra coordination [4,ntets]
        tetw     = []; % tetrahedra weight [1,ntets]
        tetv     = []; % tetrahedron volume
        % mapping
        f2i      = []; % full to irreducible
        i2f      = []; % irreducible to full
        i2b      = []; % bragg: irreducible to full
        b2i      = []; % bragg: full to irreducible
        % x-ray information
        hv       = []; % photon energy
        F        = []; % structure factor [1,nks]
        L        = []; % lorentz factor [1,nks]
        P        = []; % polarization factor [1,nks]
        I        = []; % x-ray diffracted intensity (includes lorentz polarization and weights, i.e. multiplicities) [1,nks]
    end
    
    methods (Static)

        function [bz]             = define(bas,k,varargin)
            % recbas   = []; % reciprocal basis
            % n        = []; % dimensions for monkhorst-pack mesh
            % k        = []; % kpoint coordinates
            % w        = []; % kpoint weights
            bz          = am_bz;
            bz.units    = 'frac-recp';
            bz.recbas   = inv(bas).';
            bz.k        = k;
            bz.tol      = 1E-8;
            % parse remaining input
            bz = bz.set(varargin{:});
            % populate
            bz.nks = size(bz.k,2);
            if isempty(bz.w); bz.w = ones(1,bz.nks); end
        end

        function [fbz,ibz]        = get_zone(pc,n)

            import am_lib.* am_dft.*

            % get full brillouin zone
            [fbz] = am_bz.get_fbz(pc,n);

            % get irreducible zone
            [ibz,i2f,f2i] = fbz.get_ibz(pc,'tetra');

            % save mapping to zones
            % for plist = {'f2i','i2f'}; addprop(fbz,plist{:}); end
            % for plist = {'f2i','i2f'}; addprop(ibz,plist{:}); end
            fbz.f2i = f2i; fbz.i2f = i2f;
            ibz.i2f = i2f; ibz.f2i = f2i;
        end
        
        function [fbs,ibs]        = get_bragg(uc,hv,k_max,threshold)
            % bragg brillouin zone
            import am_lib.* am_dft.*

            % k max is the largest wavevector magntiude considered
            if nargin < 2 || isempty(hv); Cu = am_atom.define('Cu'); hv = Cu.get_emission_line('kalpha1'); end
            if nargin < 3 || isempty(k_max); th2_max=180; lambda=0.15406; k_max = 2*sind(th2_max/2)/lambda; end
            if nargin < 4 || isempty(threshold); threshold = am_lib.tiny; end

            fbs = am_bz.get_fbs(uc,hv,k_max);
            
            % get structure factors and scattering intensity
                [fbs.F,fbs.L,fbs.P] = uc.get_structure_factor(fbs,hv);
                fbs.I = abs(fbs.F).^2.*fbs.L.*fbs.P.*fbs.w; fbs.I = fbs.I ./ max(fbs.I(:))*100;
                % exlcude stuff above threshold
                ex_ = fbs.I > threshold;
                [fbs.F,fbs.I,fbs.L,fbs.P,fbs.w,fbs.k,fbs.nks] = ...
                    deal(fbs.F(ex_),fbs.I(ex_),fbs.L(ex_),fbs.P(ex_),fbs.w(ex_),fbs.k(:,ex_),sum(ex_));
            
            % get symmetrically equivalent bragg spots
            % save "irreducible bragg spots" structure
            [ibs,i2b,b2i] = get_ibz(fbs,uc,'nomod,addinv');
            % for plist = {'i2b','b2i'}; addprop(ibs,plist{:}); end
            % for plist = {'i2b','b2i'}; addprop(fbs,plist{:}); end
            ibs.i2b = i2b; ibs.b2i = b2i;
            fbs.i2b = i2b; fbs.b2i = b2i; 
        end

        function [bzp]            = get_path(pc,n,brav)

            import am_lib.* am_dft.*
            
            if isnumeric(pc); bas = pc; pc=[]; pc.bas=bas; end

            % define kpoint path
            if     contains( lower(brav), 'fcc-short' )
                G=[0;0;0];  X1=[0;1;1]/2; X2=[2;1;1]/2;
                L=[1;1;1]/2; K=[6;3;3]/8;
                xl={'G','X','K','G','L'};
                ks=[G,X2,K,G];
                ke=[X1,K,G,L];
            elseif contains( lower(brav), 'fcc-long' )
                G=[0;0;0];  X1=[0;1;1]/2; W=[1;3;2]/4;
                U=[2;5;5]/8; L=[1;1;1]/2; K=[3;6;3]/8;
                xl={'G','X','W','K','G','L','U','W','L','K'};
                ks=[G,X1,W,K,G,L,U,W,L];
                ke=[X1,W,K,G,L,U,W,L,K];
            elseif contains( lower(brav), 'tetra' )
                G=[0;0;0];   Z=[0;0;1]/2; A=[1;1;1]/2; 
                M=[1;1;0]/2; X=[0;1;0]/2; R=[0;1;1]/2;
                xl={'G','X','M','G','Z','R','A','Z'};
                ks=[G,X,M,G,Z,R,A];
                ke=[X,M,G,Z,R,A,Z];
            elseif contains( lower(brav), 'sc' )
                G=[0;0;0];   X=[0;1;0]/2;
                M=[1;1;0]/2; R=[1;1;1]/2;
                xl={'G','X','M','G','R','X'};
                ks=[G,X,M,G,R];
                ke=[X,M,G,R,X];
            elseif contains( lower(brav), 'hex' )
                % for a pc.bas ordered like so:
                %     3.0531   -1.5266         0
                %          0    2.6441         0
                %          0         0    3.4526
                G=[0;0;0];   K=[1/3;1/3;0];   M=[1/2;0;0];
                A=[0;0;1/2]; H=[1/3;1/3;1/2]; L=[1/2;0;1/2];
                xl={'G','K','M','G','A','H','L','A'};
                ks=[G,K,M,G,A,H,L];
                ke=[K,M,G,A,H,L,A];
            else
                error('invalid bravais lattice');
            end

            % get number of kpoints
            nqs=size(ks,2); recbas = inv(pc.bas).';

            % get path: convert to [cart-recp] to get x spacings right then convert back to [frac-recp]
            [k,x,xt] = get_path(recbas*ks,recbas*ke,nqs,n); k=recbas\k;

            % create path object
            bzp_ = @(recbas,xl,xt,nks,x,k) struct('units','frac', ...
                'recbas',recbas,'xl',{{xl{:}}},'xt',xt,'nks',nks,'x',x,'k',k);
            bzp = bzp_(recbas,xl,xt,size(k,2),x,k);

            function [k,x,xt] = get_path(qs,qe,nqs,N)
              % define path (includes both boundaries)
                path_ = @(k,q,N) cumsum([zeros(3,1),repmat((k-q)/(N-1),1,N-1)],2)+repmat(q,1,N);
                x_    = @(k,q,N) [0, repmat(norm((k-q)/(N-1)),1,N-1) ];

                % build path
                nks = N*nqs; k = zeros(3,nks); x = zeros(1,nks);
                for i = 1:nqs
                    k(1:3,[1:N]+(i-1)*N) = path_(qe(:,i),qs(:,i),N);
                    x(    [1:N]+(i-1)*N) =    x_(qe(:,i),qs(:,i),N);
                end
                x = cumsum(x);

                % set labels coordinates
                xt = x([1,N*[1:nqs]]);
            end
        end

        function [bzs]            = get_surf(pc,n,vy,vx)
            % bzs=get_surf(pc,[101,101],[1;0;0],[0;1;0]);
            % bzs=get_bvk_dispersion(bvk,bzs);
            % plot_surf(bzs,1)

            import am_lib.* am_dft.*

            % get number of kpoints
            nks=prod(n); recbas = inv(pc.bas).';

            % get surface in [cart-recp] then convert back to [frac-recp]
            Q_ = @(n) [0:(n-1)]./(n-1); [Y{1:2}]=meshgrid(Q_(n(2)),Q_(n(1)));
            k = Y{1}(:).'.*vx + Y{2}(:).'.*vy; k = recbas\k;

            % create path object
            bzs_ = @(recbas,n,nks,k) struct('units','frac', ...
                'recbas',recbas,'nks',nks,'n',n,'k',k);
            bzs = bzs_(recbas,n,nks,k);
        end

        function [bzl]            = get_line(pc,n,vs,ve)
            % surf: vs and ve are the start and end vectors in cart

            import am_lib.* am_dft.*

            % get number of kpoints
            nks=n; recbas = inv(pc.bas).';

            % define path (includes both boundaries)
            path_ = @(k,q,N) cumsum([zeros(3,1),repmat((k-q)/(N-1),1,N-1)],2)+repmat(q,1,N);
            x_    = @(k,q,N) [0, repmat(norm((k-q)/(N-1)),1,N-1) ];

            % build path in recp.-frac
            k(1:3,:) = path_(ve-vs,[0;0;0],n); k = recbas\k;
            x(    :) =    x_(ve-vs,[0;0;0],n); x = cumsum(x);

            % create path object
            bzl_ = @(recbas,n,nks,x,k) struct('units','frac',...
                'recbas',recbas,'nks',nks,'n',n,'x',x,'k',k);
            bzl = bzl_(recbas,n,nks,x,k);
        end

        function [bza]            = get_angles(pc,hv,th2)

            import am_lib.* am_dft.*

            % get number of kpoints
            nks=numel(th2); recbas = inv(pc.bas).';

            % convert 2-theta into k points
            [kx,kz]=angle2kxkz(th2/2,th2,hv,'w2th'); k = [kx(:).';zeros(1,nks);kz(:).']; 

            % build path in recp.-frac
            k = recbas\k;

            bza = am_bz.define(pc.bas,k);
            bza.x = th2(:).';
        end

    end
    
    methods 
        
        function [bz]             = set(bz,varargin) % set bz properties 
            if isempty(varargin); return; end
            % parse remaining input
            nvs = numel(varargin);
            if ~am_lib.iseven_(nvs); error('{key, value} inputs required'); end; i=1; 
            while true; bz.(varargin{i}) = varargin{i+1}; i=i+2; if i>nvs; break; end; end
        end
        
        function [ibz,i2f,f2i]    = get_ibz(fbz,pc,flag) % get irreducible zone 

            import am_lib.* 

            if nargin<3; flag=''; end

            % get point symmetries [real-frac --> rec-frac] by transposing R
            [~,~,~,R] = pc.get_symmetries(); R = permute(R,[2,1,3]);

            % add inversion for time-reversal?
            if contains(flag,'addinv') 
            if ~any(all(all(R==-eye(3),1),2))
                R = am_cell.complete_group( cat(3,R,-eye(3)) );
            end
            end
            
            % build permutation matrix for kpoints related by point symmetries
            % nomod is used for bragg reflections, otherwise all points would go to gamma
            if contains(flag,'nomod') 
                PM = member_(    (matmul_(R,fbz.k)),fbz.k);
            else
                PM = member_(mod_(matmul_(R,fbz.k)),fbz.k);
            end
            A = get_connectivity(PM);

            % set identifiers
            i2f = round(am_lib.findrow_(A)).'; f2i = round(([1:size(A,1)]*A)); w=sum(A,2).';
            if abs(sum(w)-fbz.nks)>am_lib.eps; error('mismatch: kpoint mesh and point group symmetry'); end

            % create structure
            ibz = am_bz.define(pc.bas,fbz.k(:,i2f),'n',fbz.n,'w',w);
            
            % if requested
            if contains(flag,'tetra')
                % get all tetrahedron
                tet = get_tetrahedra(fbz.recbas,fbz.n); 
                vol = abs(det([fbz.k(:,tet(:,1));[1,1,1,1]]))/factorial(3);
                % get irreducible tetrahedra
                [tet,~,tet_f2i] = unique(sort(f2i(tet)).','rows'); 
                tet=tet.'; tetw = hist(tet_f2i,[1:size(tet,2)].'-.5);
                % augment structure with tetrahedron
                ibz = ibz.set('ntets',size(tet,2),'tet',tet,'tetw',tetw,'tetv',vol);
            end
            
            % finally, pass along additional parameters if they exist
            field = fieldnames(fbz); isfield_ = @(x) any(strcmp(field,x));
            for f = {'n','nbands','norbitals','natoms','hv'}
                if isfield_(f) && ~isempty(fbz.(f{:})); ibz.(f{:}) = fbz.(f{:}); end
            end
            for f = {'F','L','P','I','E','hw'}
                if isfield_(f) && ~isempty(fbz.(f{:})); ibz.(f{:}) = fbz.(f{:})(:,i2f); end
            end
            for f = {'lmproj','eigenv'}
                if isfield_(f) && ~isempty(fbz.(f{:})); ibz.(f{:}) = fbz.(f{:})(:,:,:,:,i2f); end
            end

            % subfunctions
            function tet           = get_tetrahedra(recbas,n)
                % divide mesh into boxes
                box = grid2box(n); nboxes = size(box,2);
                % divide a single box into six tetrahedron
                tetrahedron = box2tetrahedron(recbas);
                % loop over boxes
                tet = zeros(4,6*nboxes); t = 0;
                for b = 1:nboxes
                    % loop over tetrahedron/box
                    for j = 1:6
                        % augment tetrahedron counter
                        t = t + 1;
                        % define tetrahedra corners using indices of kpoints
                        tet(:,t) = box(tetrahedron(:,j),b);
                    end
                end

                % subfunctions
                function box           = grid2box(n)
                    % get mesh
                    [Z{1:3}]=ndgrid([1:n(1)],[1:n(2)],[1:n(3)]); ki = reshape(cat(3+1,Z{:}),[],3).';
                    % get box vertices
                    boxv = [0,1,0,1,0,1,0,1;0,0,1,1,0,0,1,1;0,0,0,0,1,1,1,1];
                    % there will be 1 box per kpoint and 8 vertices per box
                    nks = prod(n); box = zeros(8,nks);
                    % get boxes for each kpoint
                    box_ = @(d,i) mod(boxv(d,:)+ki(d,i)-1,n(d))+1;
                    for m = 1:nks; box(:,m) = sub2ind(n,box_(1,m),box_(2,m),box_(3,m)); end
                end
                function tetrahedron   = box2tetrahedron(recbas)
                    %     7-------8
                    %    /|      /|
                    %   / |     / |
                    %  5-------6  |
                    %  |  3----|--4
                    %  | /     | /
                    %  |/      |/
                    %  1-------2
                    %
                    boxvc = recbas*[0,1,0,1,0,1,0,1;0,0,1,1,0,0,1,1;0,0,0,0,1,1,1,1];
                    % get indices of diagonal pairs
                    diags=[1,2,3,4;8,7,6,5];
                    % get distances across diagonals
                    d=zeros(1,4); for m = 1:4; d(m) = norm(boxvc(:,diags(2,m))-boxvc(:,diags(1,m))); end
                    % record smallest diagonal
                    [~,si]=min(d);
                    % create connectivity list defining tetrahedra
                    switch si
                        case (1)
                        tetrahedron(:,1) = [1,8,2,4];
                        tetrahedron(:,2) = [1,8,2,6];
                        tetrahedron(:,3) = [1,8,3,4];
                        tetrahedron(:,4) = [1,8,3,7];
                        tetrahedron(:,5) = [1,8,5,6];
                        tetrahedron(:,6) = [1,8,5,7];
                        case (2)
                        tetrahedron(:,1) = [2,7,1,3];
                        tetrahedron(:,2) = [2,7,1,5];
                        tetrahedron(:,3) = [2,7,3,4];
                        tetrahedron(:,4) = [2,7,4,8];
                        tetrahedron(:,5) = [2,7,5,6];
                        tetrahedron(:,6) = [2,7,6,8];
                        case (3)
                        tetrahedron(:,1) = [3,6,1,2];
                        tetrahedron(:,2) = [3,6,1,5];
                        tetrahedron(:,3) = [3,6,2,4];
                        tetrahedron(:,4) = [3,6,4,8];
                        tetrahedron(:,5) = [3,6,5,7];
                        tetrahedron(:,6) = [3,6,7,8];
                        case (4)
                        tetrahedron(:,1) = [4,5,1,2];
                        tetrahedron(:,2) = [4,5,1,3];
                        tetrahedron(:,3) = [4,5,2,6];
                        tetrahedron(:,4) = [4,5,3,7];
                        tetrahedron(:,5) = [4,5,6,8];
                        tetrahedron(:,6) = [4,5,7,8];
                    end
                end
            end
        end

        
    end
    
    
    methods % cell properties
        
        function [nesting]        = get_nesting(fbz,ibz,degauss,Ep,flag)
            % get nesting function on ibz at probing energies Ep using smearing degauss
            % degauss, degauss = 0.04 61x61x61 kpoint mesh
            % note that these probing energies are not JDOS energies! They are simply EF energies!
            
            import am_lib.* am_dft.*

            % number of probing energies
            nEps = numel(Ep);

            % get momentum conserving q-point triplets
            qqq = get_qqq(fbz,ibz);

            % copy irreducible energies onto fbz mesh
            E = ibz2fbz(fbz,ibz,ibz.E);

            % compute spectral function A on the full mesh
            dE = (reshape(E,[1,size(E)]) - Ep(:))/degauss;
            if     contains(flag,'fermi')
                A = am_lib.fermi_dirac_dydx_(dE);
            elseif contains(flag,'mp')
                A = am_lib.methfessel_paxton_dydx_(dE,1);
            elseif contains(flag,'mv')
                A = am_lib.marzari_vanderbilt_dydx_(dE);
            elseif contains(flag,'gauss')
                A = am_lib.gauss_(dE);
            elseif contains(flag,'lorentz')
                A = am_lib.lorentz_(dE);
            else
                error('ERROR [get_dos_quick]: flag not recognized')
            end
            A = reshape( sum(A./degauss,2), [nEps,fbz.n] );

            % compute nesting on ibz and transfer to ibz
            nesting = zeros(nEps,ibz.nks); access_ = @(x,i) reshape(x(i),[],1); m = 2;
            for j = 1:nEps
                nesting(j,:) = accumarray( fbz.f2i( qqq(1,:,m)).' , ...
                           access_(A(j,:,:,:),qqq(2,:,m))     ...
                        .* access_(A(j,:,:,:),qqq(3,:,m))   , [], @sum  )./prod(ibz.n);
            end
        end
        
        function [nesting]        = get_nesting_jdos(fbz,ibz,degauss,Ep,flag)
            % get jdos nesting function on ibz at probing energies Ep (This Ep is the one that conserved energy!)
            % degauss, degauss = 0.04 61x61x61 kpoint mesh
            %
            
            import am_lib.* am_dft.*
            
            if nargin>5; flag=[]; error('only gaussian is implemented'); end
            
            % number of probing energies
            nEps = numel(Ep);

            % get momentum conserving q-point triplets
            % q1 + q2 = q3 (absorption)  q1 = q2 + q3 (emission)      : m = 1
            % q  + k  = k' (e scatters)  q  = k  + k' (e+h recombine) : m = 2
            qqq = am_dft.get_qqq(fbz,ibz);

            % copy irreducible energies onto fbz mesh
            E = reshape(am_bz.ibz2fbz(fbz,ibz,ibz.E),[ibz.nbands,fbz.n]);

            % LOL! I am trying to integrate over a ((40x40),((15x15x11)x(15x15x11))) matrix! lol... that's 9E9 numbers... LOL!!!
            % How to approach this? Loop over pairs of bands which has fewer elements than pairs of wavevectors
            % Also, some pairs of bands can be quickly ignored 
            nesting = zeros(nEps,ibz.nks); access_ = @(x,i) reshape(x(i),[],1); m = 2; nEps = numel(Ep); t=0;
            % estimate time
            for c = 1:ibz.nbands
            for v = 1:ibz.nbands
                if max(E(c,:))<0 && max(E(v,:))<0; continue; end
                if min(E(c,:))>0 && min(E(v,:))>0; continue; end
                if min(E(v,:))>max(E(c,:)); continue; end
                t=t+1;
            end
            end
            h = waitbar(0,'Integrating...'); tot=t; t=0;
            for c = 1:ibz.nbands
            for v = 1:ibz.nbands
                % skip if both Ec and Ev are in the valence band
                if max(E(c,:))<0 && max(E(v,:))<0; continue; end
                % skip if both Ec and Ev are in the conduction band
                if min(E(c,:))>0 && min(E(v,:))>0; continue; end
                % skip if Ev is strictly above Ec
                if min(E(v,:))>max(E(c,:)); continue; end
                % timer
                t=t+1; waitbar(t./tot,h,'Integrating...'); 
                % if some Ec and some Ev are conduction/valence bands, consider the situation
                Ec = access_(E(c,:,:,:),qqq(2,:,m));
                Ev = access_(E(v,:,:,:),qqq(3,:,m));
                % keep only points where Ec and Ev are not together (conduction vs valence)
                ex_ = xor(Ec<0, Ev<0); ex_(ex_) = Ec(ex_)>Ev(ex_);
                % evaluate nesting
                % loop over energies so that accumarray can be applied
                for n = 1:nEps
                    nesting(n,:) = nesting(n,:) + ...
                        accumarray( fbz.f2i( qqq(1,ex_,m)).' , ...
                                am_lib.gauss_((  (Ec(ex_)-Ev(ex_)).'-Ep(n)   )./degauss)./degauss , [] , @sum ).';
                end
                % another approach (didn't work out too well. i don't know what is wrong... )
            %     nesting = nesting + ...
            %         accumarray( flatten_( repmat(fbz.f2i(qqq(1,ex_,m)), nEps, 1 )) , ...
            %                 flatten_(am_lib.gauss_((  (Ec(ex_)-Ev(ex_)).'-Ep(:) )./degauss)./degauss) , [] , @sum ).';
            end
            end
            close(h);
            nesting = nesting./prod(fbz.n);
            
        end

        function [dos]            = get_dos(ibz,Ep,flag,degauss)
            % Ep=linspace(-6,6,1000); degauss= 0.08;
            % hold on;
            % for method = {'tet,'fermi','mp','mv','gauss','lorentz'}
            %     [D] = am_dft.get_dos(ibz,Ep,degauss,method{:});
            %     plot(Ep,D);
            % end
            % hold off;
            if nargin<3 || isempty(flag); flag = 'tet,dos'; end
            if nargin<4 || isempty(degauss); degauss = 0.05; end
            
            % get energies E and occupation weights occw
            if     contains(flag,'ojdos')
                % occupation-weighted jdos: jdos properly weighed by fermi function instead of theta function
                % and divided by 1/E.^2
                kT = 0.025852; % 300K
                % E(valence,conduction,kpoints) 
                Ev = permute(ibz.E,[1,3,2]); Ec = permute(ibz.E,[3,1,2]); E = Ec - Ev;
                % get occupational weights 
                occw = am_lib.fermi_dirac_(Ec./kT) .* ( 1 - am_lib.fermi_dirac_(Ev./kT)); 
                E = reshape(E,ibz.nbands^2,ibz.nks); occw = reshape(occw,ibz.nbands^2,ibz.nks);
            elseif contains(flag,'jdos')
                % checked against vasp. confirmed.
                % only consider transitions from the valence band (E<0) to the conduction band (E>0)
                Ev = permute(ibz.E,[1,3,2]); Ec = permute(ibz.E,[3,1,2]); E = Ec - Ev;
                % E(valence,conduction,kpoints) 
                % make all valence    - valence    states E => 1E8 (something large enough to take it out of range)
                % make all conduction - conduction states E => 1E8 (something large enough to take it out of range)
                % proper way would be to compute fermi functions and use that as the projection weight
                % but just doing it quick and dirty here (equivalent to assuming kt -> 0, fermi fucntion -> step function)
                ex_ = (Ec<0) & (Ev<0); E(ex_) = 1E8;
                ex_ = (Ec>0) & (Ev>0); E(ex_) = 1E8;
                E = reshape(E,ibz.nbands^2,ibz.nks); E = sort(E); occw = 1;
            elseif contains(flag,'dos')
                E = ibz.E; occw = 1;
            else
                error('ERROR [get_dos]: flag must include one of [ojdos, jdos, dos]')
            end
            
            % perform integration
            if contains(flag,'tet')
                % tetrahedron method
                if     contains(flag,'ojdos')
                    % sort E and occw together
                    [E,I] = sort(E); [m,n]=size(occw); occw = occw(sub2ind([m n],I,repmat(1:n,m,1))); 
                    % prepare occw for pdos
                    occw = permute(occw,[3,1,2]);
                    % compute optical jdos
                    D = am_dft.get_pdos_tet(Ep,E,ibz.tet,ibz.tetw,ibz.tetv,occw,'mex');
                elseif contains(flag,'jdos') || contains(flag,'dos')
                    D = am_dft.get_dos_tet(Ep,E,ibz.tet,ibz.tetw,ibz.tetv,'mex');
                end
            else
                % discrete integration
                Ep = permute(Ep(:),[3,2,1]); dE = ((E-Ep)./degauss);
                if     contains(flag,'fermi')
                    D = am_lib.sum_(ibz.w .* am_lib.fermi_dirac_dydx_(dE) .* occw,[1,2]);
                elseif contains(flag,'mp')
                    D = am_lib.sum_(ibz.w .* am_lib.methfessel_paxton_dydx_(dE,1) .* occw,[1,2]);
                elseif contains(flag,'mv')
                    D = am_lib.sum_(ibz.w .* am_lib.marzari_vanderbilt_dydx_(dE) .* occw,[1,2]);
                elseif contains(flag,'gauss')
                    D = am_lib.sum_(ibz.w .* am_lib.gauss_(dE) .* occw,[1,2]);
                elseif contains(flag,'lorentz')
                    D = am_lib.sum_(ibz.w .* am_lib.lorentz_(dE) .* occw,[1,2]);
                else
                    error('ERROR [get_dos]: flag not recognized')
                end
                D=D./degauss./sum(ibz.w(:));
            end
            % create structure
            dos.E=Ep(:).';
            dos.D=D(:).';
        end
 
    end
    
    methods % plotting

        function                    plot_interpolated(fbz,bzp,x, varargin)
            % interpolate x, defined on the monkhorst-pack fbz, on [frac] path
            %
            % plot_interpolated(fbz,bzp, ibz2fbz(fbz,ibz,ibz.E) ,'-');
            % hold on; plot_dispersion(tb,bzp,'electron','--');
            %

            import am_lib.fftinterp_

            % define figure properties
            fig_ = @(h)       set(h,'color','white');
            axs_ = @(h,xt,xl) set(h,'Box','on','XTick',xt,'Xticklabel',xl);
            fig_(gcf);

            % plot results
            plot(bzp.x,fftinterp_(x,bzp.k,fbz.n), varargin{:});
            axs_(gca,bzp.xt,bzp.xl); axis tight; xlabel('Wavevector k');
        end

        function                    plot_dispersion_projected(bzp,pc,flag)

            % get eigenvalues and band character weights
            c = zeros(bzp.nbands,bzp.nks); w = [];
            
            if     contains(flag,'orbital') % projected on orbitals
                character_labels = {'s','p','d','f'}; c = projections2color_(bzp,pc,'orbital');
                if contains(flag,'atom')
                    if numel(pc.i2p)~=2; error('cannot broaden width in more than 2 dimensions'); end
                    w = projections2color_(bzp,pc,'atom')/2;
                end
            elseif contains(flag,'atom') % projected on atoms
                character_labels = pc.symb; c = projections2color_(bzp,pc,'atom');
            else % just plot bands straight up
                flag = [flag,'none'];
                character_labels=''; c = ones(bzp.nbands,bzp.nks); 
            end

            % define figure properties
            fig_ = @(h)       set(h,'color','white');
            axs_ = @(h,xt,xl) set(h,'Box','on','XTick',xt,'Xticklabel',xl);
            fig_(gcf);

            % plot band structure
            fig_(gcf); 
            if contains(flag,'none')
                plot(bzp.x,bzp.E,'-','color','b');
            else
                hold on;
                for m = 1:bzp.nbands
                    if isempty(w)
                        am_lib.plotc_(bzp.x,bzp.E(m,:),c(m,:));
                    else
                        am_lib.plotc_(bzp.x,bzp.E(m,:),c(m,:),w(m,:));
                        plot(bzp.x,bzp.E,'-w','linewidth',0.1);
                    end
                end
                hold off;
                % apply color map and label axes
                colormap( am_lib.colormap_('spectral',100).^(2) ); h = colorbar; caxis([0,1]);
                cticks = am_lib.assign_cmap_(eye(numel(character_labels))); [~,inds] = sort( cticks );
                set(h,'Ticks',cticks(inds),'TickLabels',character_labels(inds));
            end
            
            % label axes
            axs_(gca,bzp.xt,bzp.xl); axis tight; ylabel('Energy [eV]'); xlabel('Wavevector k');

            % plot fermi level
            line([0,bzp.x(end)],[0,0],'linewidth',2,'color',[1,1,1]*0.5,'linestyle',':');
            
            function c = projections2color_(bz,pc,flag)
                % normalize columns of matrixs
                normc_ = @(m) ones(size(m,1),1)*sqrt(1./sum(m.*m)).*m;
                % colorize
                if     contains(flag,'orbital')
                for i = 1:bz.nks
                    %    s     py     pz     px    dxy    dyz    dz2    dxz    dx2    f-3    f-2    f-1     f0     f1     f2     f3
                    % l-PROJECTION (sum over spin, atoms, m-quantum numbers)
                    % lmproj(nspins,norbitals,natoms,nbands,nkpts)
                    Vp(1,:) = squeeze(sum(sum(sum(bz.lmproj(:,   [1],:,:,i),1),2),3)); % s
                    Vp(2,:) = squeeze(sum(sum(sum(bz.lmproj(:, [2:4],:,:,i),1),2),3)); % p
                    Vp(3,:) = squeeze(sum(sum(sum(bz.lmproj(:, [5:9],:,:,i),1),2),3)); % d
                    Vp(4,:) = squeeze(sum(sum(sum(bz.lmproj(:,[9:16],:,:,i),1),2),3)); % f
                    c(:,i) = am_lib.assign_cmap_(normc_(Vp));
                end
                elseif contains(flag,'atom')
                    for i = 1:bz.nks
                        % l-PROJECTION (sum over spin, atoms, m-quantum numbers)
                        % lmproj(nspins,norbitals,natoms,nbands,nkpts)
                        for j = 1:numel(unique(pc.p2i))
                            Vp(j,:) = squeeze(sum(sum(sum(bz.lmproj(:,:,pc.p2i==j,:,i),1),2),3));
                        end
                        c(:,i) = am_lib.assign_cmap_(normc_(Vp));
                    end
                end
            
            end
        end

        function                    plot_nesting(ibz,fbz,bzp,degauss,Ep,varargin)
            import am_lib.*
            % plot results
            plot_interpolated(fbz,bzp, am_bz.ibz2fbz(fbz,ibz,am_bz.get_nesting(fbz,ibz,degauss,Ep)) , varargin{:})
        end
        
        function [Aibz]          =  plot_nesting_jdos(fbz,ibz,bzp,degauss,Ep)
            
            Aibz = am_dft.get_nesting_jdos(fbz,ibz,degauss,Ep);
            Afbz = am_dft.ibz2fbz(fbz,ibz,Aibz);
            Abzp = am_lib.fftinterp_(Afbz,bzp.k,fbz.n);
            
            subplot(3,1,1:2);
            [EE,XX]=ndgrid(Ep,bzp.x);
            % ex_ = ones(size(Abzp)); ex_(Abzp<0.01) = 0;
            % surf(XX,EE*scale_energies_,(abs(Abzp)./Ep.^2).^(0.1).*ex_,'edgecolor','none'); view([0 0 1]); axis tight;
            surf(XX,EE,abs(Abzp),'edgecolor','none'); view([0 0 1]); axis tight;

            axs_ = @(h,xt,xl) set(h,'Box','on','XTick',xt,'Xticklabel',xl);
            axs_(gca,bzp.xt,bzp.xl); axis tight; ylabel('Energy [eV]'); xlabel('Wavevector k');

        end

        function                    plot_bz(fbz)

            import am_lib.* am_dft.*

            % initialize figure
            set(gcf,'color','w'); hold on;

            % plot points
            % h = scatter3_(uc2ws(fbz.recbas*fbz.k,fbz.recbas),'.');

            % generate all possible points halfway between reciprocal lattice vectors
            P = [-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1;...
                 -1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1;...
                 -1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1]/2;
            P = uc2ws(fbz.recbas*P,fbz.recbas); P = uniquec_(P);

            % get the wigner-seitz planes corresponding to each point
            N = P./normc_(P); nplanes=size(P,2); NP=zeros(4,nplanes);
            for i = 1:nplanes; NP(:,i) = [N(:,i);N(:,i).'*P(:,i)]; end

            % find all combinations of three wigner-seitz planes which intersect at a point
            ijk = nchoosek_(size(P,2),3); nijks = size(ijk,2); ex_ = false(1,nijks);
            for i = 1:nijks; ex_(i) = (rank(NP(1:3,ijk(:,i)).')==3); end

            % find points at the intersection of three wigner-seitz planes
            X=zeros(3,sum(ex_));
            for i = find(ex_); X(:,i) = NP(1:3,ijk(:,i)).'\NP(4,ijk(:,i)).'; end

            % shift to wigner seitz cell and get unique values [relax the edge by ~0.99999x]
            X = uc2ws(X*0.999,fbz.recbas)/0.999; X = uniquec_(X);

            % plot convex hull
            plothull_(X);

            hold off; daspect([1 1 1]); box on;
        end

        function                    plot_path(bzp)

            import am_lib.* am_dft.*

            % initialize figure
            set(gcf,'color','w'); hold on;

            % plot brillouin zone path
            k=bzp.recbas*bzp.k;
            hold on; plot3(k(1,:),k(2,:),k(3,:),':r','linewidth',2);

            % plot brillouin zone boundary
            hold on; plot_bz(bzp);

            % plot reciprocal lattice vectors
            hold on; plotv3_(bzp.recbas,'o-','linewidth',2);

            hold off; daspect([1 1 1]); box on;
        end

        function                    plot_surf(bzs,band)
            % bzs=get_surf(pc,[101,101],[1;0;0],[0;1;0]);
            % bzs=get_bvk_dispersion(bvk,bzs);
            % plot_surf(bzs,1)

            import am_lib.* am_dft.*

            % initialize figure
            set(gcf,'color','w'); hold on;

            % plot bz surface
            k=bzs.recbas*bzs.k; k(abs(k)<am_lib.eps)=0; s_ = @(i) reshape(k(i,:),bzs.n);
            surf(s_(1),s_(2),s_(3),reshape(real(bzs.hw(band,:)),bzs.n),'facecolor','interp');
            shading interp;

            % plot brillouin zone boundary
            plot_bz(bzs);

            % fix axes
            minmax = @(x) [min(x(:)),max(x(:))];
            caxis(minmax(bzs.hw(band,:)))

            hold off; daspect([1 1 1]); box on;
        end

        function h                = plot_fermi_surface(fbz,ibz,Ef,flag,C)
            if nargin < 4 || isempty(flag); flag=''; end
            if nargin < 5 || isempty(C); C=[]; end % C is the colorcode
            
            % % color based on curvature? (each band is a scalar field: hessian is diagonal)
            % r = ndgrid(1:fbz.n(1),1:fbz.n(2),1:fbz.n(3)); r = permute(r,[4,1,2,3]);
            % for i = 1:dft.nbands; C(i,:,:,:) = ifftn( fftn(E(i,:,:,:)).*(1i*r).^2 ); end; C = abs(real(C));
            
            % if C is a string, color based on field
            if contains(flag,'optmat')
                if isempty(ibz.optmat); error('field is empty'); end
                C = am_lib.sum_(ibz.optmat,[1,2])/3; % sum over spins and average over directions
                C = reshape( am_bz.ibz2fbz(fbz,ibz,reshape(C,[ibz.nbands.^2,ibz.nks])),[ibz.nbands.^2,fbz.n]);
            end
            
            if contains(flag,'jdos')
                % checked against vasp. confirmed.
                % only consider transitions from the valence band (E<0) to the conduction band (E>0)
                Ev = permute(ibz.E,[1,3,2]); Ec = permute(ibz.E,[3,1,2]); E = Ec - Ev;
                % E(valence,conduction,kpoints) 
                % make all valence    - valence    states E => 1E8 (something large enough to take it out of range)
                % make all conduction - conduction states E => 1E8 (something large enough to take it out of range)
                ex_ = (Ec<0) & (Ev<0); E(ex_) = 1E8;
                ex_ = (Ec>0) & (Ev>0); E(ex_) = 1E8;
                E = reshape( E , [ibz.nbands^2, ibz.nks] ); E = sort(E);
                E = reshape( am_bz.ibz2fbz(fbz,ibz,E), [ibz.nbands.^2, fbz.n]); 
            else % just regular dos
                E = reshape( am_bz.ibz2fbz(fbz,ibz,ibz.E), [ibz.nbands, fbz.n]); 
            end

            % plot                                                                                 [4,4,1]  2
            h = am_lib.plot_isosurface_(reshape(fbz.recbas*fbz.k,[3,fbz.n]), fbz.recbas, E, Ef, C, [1,1,1], 3, 'cspline,center');
            % plot bz boundary
            if contains(flag,'tetragonal'); fbz.draw_boundary('tetragonal'); end
            % set lighting settings
            daspect([1 1 1]); axis tight; axis off; campos([10.7135,23.0003,9.4299]); camlight('right','infinite'); camproj('perspective');
        end
        
        function h                = draw_boundary(fbz,flag)
            switch flag
                case 'tetragonal'; h = draw_tetragonal_boundary(fbz);
                otherwise; error('unknown geometry');
            end
            
            function h = draw_tetragonal_boundary(fbz)
                R3 = [ 0,-1, 0;-1, 0, 0; 0, 0,-1];
                s = [1,1,1];
                xy = fbz.recbas*[[0;0;0],[s(1);0;0],[s(1);s(2);0],[0;s(2);0]];
                xz = fbz.recbas*[[0;0;0],[0;0;s(3)],[s(1);0;s(3)],[s(1);0;0]];
                alpha = 0.015; i=0;
                i=i+1; h(i) = patch('Faces',[1:4],'Vertices',(xy)','EdgeColor','k','FaceColor','k','facealpha',alpha,'LineWidth',1);
                i=i+1; h(i) = patch('Faces',[1:4],'Vertices',(xz)','EdgeColor','k','FaceColor','k','facealpha',alpha,'LineWidth',1);
                i=i+1; h(i) = patch('Faces',[1:4],'Vertices',(xy+fbz.recbas*[0;0;s(3)])','EdgeColor','k','FaceColor','k','facealpha',alpha,'LineWidth',1);
                i=i+1; h(i) = patch('Faces',[1:4],'Vertices',(xz+fbz.recbas*[0;s(2);0])','EdgeColor','k','FaceColor','k','facealpha',alpha,'LineWidth',1);
                i=i+1; h(i) = patch('Faces',[1:4],'Vertices',(R3*xz+fbz.recbas*[0;s(2);s(3)])','EdgeColor','k','FaceColor','k','facealpha',alpha,'LineWidth',1);
                i=i+1; h(i) = patch('Faces',[1:4],'Vertices',(R3*xz+fbz.recbas*[s(1);s(2);s(3)])','EdgeColor','k','FaceColor','k','facealpha',alpha,'LineWidth',1);
            end
        end
        
        % plot x-ray stuff

        function                    print_bragg_table(fbs,hv)
            th_ = @(hv,q) asind(am_lib.get_photon_wavelength(hv)*q/2);
            % print results
            fprintf('     %5s %5s %5s %10s %5s %10s %10s %10s %15s\n','h','k','l','2th [deg]','w','L','P','Fhkl^2 [%]','w*L*P*Fhkl^2 [%]');
            fprintf('     %5s %5s %5s %10s %5s %10s %10s %10s %15s\n','-----','-----','-----','----------','-----','----------','----------','----------','---------------');
            for j = 1:fbs.nks
                % exclude everything with peak height smaller than threshold
                th2 = 2 * th_( hv, norm(fbs.recbas*fbs.k(:,j)) );
                fprintf('     %5i %5i %5i %10.3f %5i %10.3f %10.3f %10.3f %15.3f\n',fbs.k(:,j),th2,fbs.w(j),fbs.L(j),fbs.P(j),fbs.Fk2(j),fbs.intensity(j));
            end
        end
        
        function [h]              = plot_ibs_1D(ibs,labels,hv,label_threshold,varargin)
            import am_lib.* am_dft.*
            % set default maximum range in plot
            max_th2_range = 110;
            %
            % yscaler_ = @(x) x.^0.1;
            yscaler_ = @(x) x;
            % set threshold
            if nargin<4 || isempty(label_threshold); label_threshold = 0; end
            if nargin<2 || isempty(labels); labels = cell(1,numel(ibs)); end
            % number of bragg structures (one for each cell)
            nibss=numel(ibs);
            % plot results
            if nibss>1
                clist = am_lib.cmap_('spectral',nibss).';
                for j = 1:nibss
                    ax(j) = axes('position',[0.025 (0.1+0.87*(j-1)/nibss) 0.95 0.85/nibss]);
                    h=plot_ibs_1D(ibs(j),[],label_threshold,'color',clist(:,j)); 
                    if j~=1; set(gca,'XTickLabel',[]); xlabel(''); end
                    % set y axis label properties
                    ax(j).YLabel.String=labels{j};
                    ax(j).YLabel.Color=clist(:,j);
                end
                linkaxes(ax);
            else
                % plot Bragg peaks
                set(gcf,'color','w'); 
                h=hggroup;
                for i = 1:ibs.nks
                    th2 = 2*get_th(norm(ibs.recbas*ibs.k(:,i)),hv);
                    if th2 < max_th2_range
                        line([th2,th2],[0,yscaler_(abs(ibs.F(i)).^2)],'linewidth',2,varargin{:},'Parent',h);
                        if abs(ibs.F(i)).^2 > label_threshold
                            % text(th2,ibs.Fk2(i),sprintf('  [%i%i%i]  %.2f^\\circ  %i', ibs.k(:,i),th2,ibs.w(i)),'Parent',h);
                            text(th2,abs(ibs.F(i)).^2+2,sprintf('[%i%i%i]  %.2f^\\circ', ibs.k(:,i),th2),'EdgeColor','k','BackgroundColor','w','Parent',h);
                            % text(th2,yscaler_(ibs.Fk2(i)),sprintf('  [%i%i%i]', ibs.k(:,i)),varargin{:},'Parent',h);
                        end
                    end
                end
                box on; xlabel('2\theta [deg]'); xlim([5 max_th2_range]); ylim([0 yscaler_(130)]); %xlabel('intensity [a.u.]');
                set(gca,'XTick',[0:10:max_th2_range]); set(gca,'YTick',[]); grid on; set(gca,'YMinorGrid','on');
            end
            H=findobj(gca,'Type','text');
            set(H,'Rotation',90);
        end
        
        function [h]              = plot_ibs_2D(ibs)
            % uc=am_dft.load_cell('material','TiN');
            % [fbs,ibs] = get_bragg(uc);
            % hv = get_atomic_emission_line_energy(get_atomic_number('Cu'),'kalpha1');
            % [ibs] = get_structure_factor(uc,ibs,hv); print_bragg_table(ibs,hv);
            % [fbs] = get_structure_factor(uc,fbs,hv);
            % % h = plot_ibs_1D(ibs,[],0.1)
            % % [h]       = plot_ibs_2D(ibs)
            % v1=[0 0 1]; v2=[1 0 0];  
            % plot_fbs_2D(fbs,v1,v2)
            % include edge points 
            N = 201; euler = exp(2i*pi*[0:N]/N);
            x = real(euler);y = imag(euler); ibs.Fk2 = abs(ibs.F).^2;
            figure(1); set(gcf,'color','w'); h=hggroup;
            for i = 1:ibs.nks
                r = norm(ibs.recbas*ibs.k(:,i));
                line(r*x,r*y,'color','r','linewidth',log10(ibs.Fk2(i)+1)/1.5+0.5,'Parent',h)
                text(0,r,sprintf('%i %i %i',ibs.k(:,i)),'HorizontalAlignment','center','BackgroundColor','w','EdgeColor','r','Parent',h);
            end
            daspect([1 1 1]); axis tight; box on;
        end
        
        function [h]              = plot_fbs_2D(fbs,v1,v2,flag)
            import am_lib.*
            %
            if nargin<4; flag=''; end
            % check for orthogonality
            if ~am_lib.isorthogonal_(v1,v2); error('v1 and v2 must be orthogonal'); end
            % normalize v1 and v2
            v1 = v1(:)./norm(v1); v2 = v2(:)./norm(v2);
            % see which points lie on the plane
            for i = 1:fbs.nks
                if eq_(det([v1(:),v2(:),fbs.k(:,i)]),0)
                    ex_(i) = true;
                else
                    ex_(i) = false;
                end
            end
            % exclude points not on the plane
            fbs.nks = sum(ex_); [fbs.k,fbs.w,fbs.F,fbs.L,fbs.P,fbs.I,fbs.b2i] = ...
                deal(fbs.k(:,ex_),fbs.w(ex_),fbs.F(ex_),fbs.L(ex_),fbs.P(ex_),fbs.I(ex_),fbs.b2i(ex_)); 
            % convert to cartesian coordinates
            v1_cart = fbs.recbas*v1(:); v1_cart = v1_cart./normc_(v1_cart); 
            v2_cart = fbs.recbas*v2(:); v2_cart = v2_cart./normc_(v2_cart); 
            k_cart  = fbs.recbas*fbs.k; 
            % project points onto the diffraction plane
            fbs.x = [v1_cart.'*k_cart; v2_cart.'*k_cart];
            % align points v2 to point up?
            if contains(flag,'align')
                % get tilt
                % v1_cart * v2 is correct!
                tilt = atan2d(v1_cart.'*fbs.recbas*v2,v2_cart.'*fbs.recbas*v2);
                % get rep
                R = am_lib.rotzd_( tilt ); R = R(1:2,1:2);
                % transform
                fbs.x = R * fbs.x;
            end
            % plot only upper half?
            if contains(flag,'xrd')
                ex_ = isaccessible_(fbs.x(1,:),fbs.x(2,:),fbs.hv); 
                % ex_ = fbs.x(2,:)>-fbs.tol; % plot only upper half
                fbs.nks = sum(ex_); [fbs.k,fbs.w,fbs.x,fbs.F,fbs.L,fbs.P,fbs.I,fbs.b2i] = ...
                    deal(fbs.k(:,ex_),fbs.w(ex_),fbs.x(:,ex_),fbs.F(ex_),fbs.L(ex_),fbs.P(ex_),fbs.I(ex_),fbs.b2i(ex_));
            end
            % plot plane with points
            figure(1); set(gcf,'color','w'); % h=hggroup; 
            hold on; 
                h = scatter(fbs.x(1,:),fbs.x(2,:),rescale(log10(abs(fbs.F).^2),20,100)/2,'filled','linewidth',1.5);%,'Parent',h); 
                scatter(0,0,log10(100)*100,'k','filled','linewidth',2);%,'Parent',h); 
            hold off;
            if contains(flag,'legend')
                for i = 1:fbs.nks
                    text(fbs.x(1,i),fbs.x(2,i)+0.1,sprintf('%i %i %i',fbs.k(:,i)), ...
                        'HorizontalAlignment','center','BackgroundColor','w','EdgeColor',h.CData);
                end
            end
            daspect([1 1 1]); axis tight; box on;
            
            function [ex_]            = isaccessible_(kx,kz,hv)
                    [alpha_i,alpha_f] = am_lib.kxkz2angle(kx,kz,hv,'in/exit');
                    ex_ = true(size(alpha_i));
                    ex_(ex_) = alpha_i(ex_) < 180;
                    ex_(ex_) = alpha_i(ex_) >   0;
                    ex_(ex_) = alpha_f(ex_) < 180 - alpha_i(ex_);
                    ex_(ex_) = alpha_f(ex_) >   0;
            end
        end
        
        function [h]              = plot_xrd_accessible(fbs,flag)
            if nargin < 2; flag='xrd'; end 
            hv = fbs.hv; 
            hold all;
                switch flag
                    case 'xrr'; N = 50; m = 5;
                        alpha = 0:0.1:5; 
                        for i = 1:numel(alpha)
                            [kx,kz] = am_lib.angle2kxkz(ones(1,N)*alpha(i),linspace(0,m,N),hv,'in/exit'); line(kx,kz,'linestyle',':');
                            [kx,kz] = am_lib.angle2kxkz(linspace(0,m,N),ones(1,N)*alpha(i),hv,'in/exit'); line(kx,kz);
                        end
                        grid off; view([0 0 1]); box on;
                    case 'xrd'; N = 50; m = 180;
                        % generate accessible-region
                        alpha_i=zeros(N); alpha_f=zeros(N);
                        alpha_i = repmat(linspace(0,m,N),N,1); 
                        for i = 1:N; alpha_f(:,i) = linspace(0,m-alpha_i(1,i),N); end
                        w = alpha_i; th2 = alpha_i+alpha_f;
                        % convert to reciprocal coordinates
                        [kx,kz] = am_lib.angle2kxkz(w,th2,hv,'w2th');
                        % plot 
                        % h = mesh(kx,kz,log(w*1E-8),'EdgeColor',[1 1 1],'linewidth',0.2);
                        h = mesh(kx,kz,log(w*1E-8),'facecolor','w'); %,'EdgeColor',[1 1 1]*0.85);
                        grid off; view([0 0 1]); box on; colormap('gray'); 
                    otherwise
                        error('unknown option');
                end
            hold off;
        end
        
    end
    
    methods (Static, Access = protected) % internal

        function [fbz]            = get_fbz(pc,n) % n = either kpoint mesh [n1,n2,n3] or density [n]
            switch numel(n)
                case 1 % convert kpoint density to odd grid
                    n_ = @(bas,s) 2*ceil(am_lib.normc_(inv(bas).')*s/2)+1; 
                    fbz = am_bz.get_fbz(pc,n_(pc.bas,n)); 
                    return;
                case 3
                    if any(mod(n,1)~=0); error('n must be integers'); end
                    % generate primitive lattice vectors
                    Q_ = @(i) [0:(n(i)-1)]./n(i); [Y{1:3}]=ndgrid(Q_(1),Q_(2),Q_(3)); k=reshape(cat(3+1,Y{:}),[],3).';
                    % create structure
                    fbz = am_bz.define(pc.bas,k,'n',n);
                otherwise
                    error('unknown n');
            end
        end

        function [fbs]            = get_fbs(pc,hv,k_max,N)
            
            import am_lib.* 
            
            if nargin < 2 || isempty(hv); Cu = am_atom.define('Cu'); hv = Cu.get_emission_line('kalpha1'); end
            if nargin < 3 || isempty(k_max); th2_max=180; lambda=0.15406; k_max = 2*sind(th2_max/2)/lambda; end
            if nargin < 4 || isempty(N); N = ceil((max(k_max./(1./am_lib.normc_(pc.bas))))); end 
            
            % get miller indicies [hkl] excluding gamma
            k=permn_([N:-1:-N],3).'; k=k(:,~all(k==0,1)); 
            % get reciprocal basis vectors [1/nm]
            recbas = inv(pc.bas).';
            % sort by distance
            k = k(:,rankc_(normc_(recbas*k)));               
            % identify values which cannot be reached by diffractometer
            k = k(:,normc_(recbas*k)<k_max); 
            % save "full bragg spots" structure 
            fbs = am_bz.define(pc.bas,k); fbs.hv = hv;
        end
        
    end

    methods (Static, Access = protected) % aux brillouin zones

        function [y]              = ibz2fbz(fbz,ibz,x) % copy ibz values on fbz grid 
            % copy ibz values on fbz grid
            %
            %       x [ n , ibz.nks ] ==> y [ n , fbz.nks ]
            %
            % n can be any number of rows
            %
            y = zeros(size(x,1),fbz.nks);
            for i = [1:ibz.nks]; y(:,fbz.f2i==i) = repmat( x(:,i) , [1,sum(fbz.f2i==i)] ); end
        end
        
        function [qqq]            = get_qqq(fbz,ibz) % get_qqq(fbz,ibz) get momentum-conserving vectors for three-scattering processes 
            % get all possible wavevector triplets which conserve momentum
            %        q1 + q2 = q3 (absorbtion) and q1 = q2 + q3 (emission)
            %        q1 are ibz points
            %
            %        qqq( 1:3 , ibz.nks , [1(abs):2(ems)])
            %
            % Check with:
            %     % NOTE: umklapp triplets are given by values which a reciprocal lattice
            %     % vector not equal to [0;0;0], i.e. [0;0;1].
            %     % fbz.k(:,qqq(1,:))+fbz.k(:,qqq(2,:))-fbz.k(:,qqq(3,:))
            %     check_ = @(x) all(mod_(x(:))<am_lib.eps);
            %     check_( fbz.k(:,qqq(1,:,1))+fbz.k(:,qqq(2,:,1))-fbz.k(:,qqq(3,:,1)) )
            %     check_( fbz.k(:,qqq(1,:,2))-fbz.k(:,qqq(2,:,2))-fbz.k(:,qqq(3,:,2)) )

            import am_lib.mod_

            sub2inds_ = @(n,q) reshape( ...
                               sub2ind(n,round(q(1,:).*n(1)+1), ...
                                         round(q(2,:).*n(2)+1), ...
                                         round(q(3,:).*n(3)+1)), 1, size(q,2)*size(q,3)*size(q,4) );

            % get wavevector triplets
            ex_ = repelem([1:ibz.nks],fbz.nks);
            qqq = zeros(3,ibz.nks.*fbz.nks,2);
            for i = [1:ibz.nks]
                % Procedure: get q1, q2, and q3 vector satisfying:
                %     q1 + q2 = q3
                %     q1 + reshape(fbz.k,[3,fbz.n]) - circshift(reshape(fbz.k,[3,fbz.n]),-[0,[q{:}]-1])
                %     q1 = q2 + (-q3)
                %     q1 - reshape(fbz.k,[3,fbz.n]) + circshift(reshape(fbz.k,[3,fbz.n]),+[0,[q{:}]-1])
                %
                % get shift corresponding to ibz point
                [q{1:3}]=ind2sub(fbz.n,fbz.i2f(i));
                % get q1
                q1 = ([q{:}].'-1)./fbz.n(:);
                % set q2
                q2 = reshape(fbz.k,[3,fbz.n]);
                % get q3: q1 + q2 = q3 (absorption)  q1 = q2 + q3 (emission)
                %         q  + k  = k' (e scatters)  q  = k  + k' (e+h recombine)
                q3_abs =       circshift(reshape(fbz.k,[3,fbz.n]),-[0,[q{:}]-1]);
                q3_ems = mod_(-circshift(reshape(fbz.k,[3,fbz.n]),+[0,[q{:}]-1]));
                % save indicies
                qqq(1,ex_==i,1:2) = repmat(sub2inds_(fbz.n,q1),[1,fbz.nks,2]);
                qqq(2,ex_==i,1:2) = repmat(sub2inds_(fbz.n,q2),[1,1,2]      );
                qqq(3,ex_==i,1  ) =        sub2inds_(fbz.n,q3_abs)           ;
                qqq(3,ex_==i,2  ) =        sub2inds_(fbz.n,q3_ems)           ;
            end
        end

    end
    
end
