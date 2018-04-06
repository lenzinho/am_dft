classdef am_bz
    
    properties
        % basic properties
        units    = []; % frac/cart
        recbas   = []; % reciprocal basis
        n        = []; % dimensions for monkhorst-pack mesh
        nks      = []; % number of kpoints
        k        = []; % kpoint coordinates
        w        = []; % kpoint weights
        E        = []; % electron dispersion
        hw       = []; % phonon dispersion
        tol      = []; % numerical tolerance
        % tetrahedron information
        ntets    = []; % number of tetrahedra
        tet      = []; % tetrahedra coordination
        tetw     = []; % tetrahedra weight
        tetv     = []; % tetrahedra volume
        % x-ray related info
    end
    
    methods (Static)

        function [uc]            = define(pc,n,k,w) % define(pc,n,k,w)
            % create brillouin zone
            if nargin<4; w=ones(1,size(k,2)); end
            % recbas   = []; % reciprocal basis
            % n        = []; % dimensions for monkhorst-pack mesh
            % k        = []; % kpoint coordinates
            % w        = []; % kpoint weights
            uc          = am_bz;
            uc.units    = 'frac';
            uc.recbas   = inv(pc.bas).';
            uc.n        = n;
            uc.k        = k;
            uc.w        = w;
        end
       
    end
    
    methods  % convert between cells

        function [fbz,ibz]        = get_zones(pc,n)

            import am_lib.* am_dft.*

            % get full brillouin zone
            [fbz] = get_fbz(pc,n);

            % get irreducible zone
            [ibz,i2f,f2i] = get_ibz(fbz,pc,'tetra');

            % save mapping to zones
            fbz.f2i = f2i; fbz.i2f = i2f;
            ibz.i2f = i2f; ibz.f2i = f2i;
        end
        
        function [fbz]            = get_fbz(pc,n) % n = either kpoint mesh [n1,n2,n3] or density [n]
            switch numel(n)
                case 1 % convert kpoint density to odd grid
                    n_ = @(bas,s) 2*ceil(am_lib.normc_(inv(bas).')*s/2)+1; fbz = get_fbz(pc,n_(bz.bas,n));
                case 3
                    if any(mod(n,1)~=0); error('n must be integers'); end
                    % generate primitive lattice vectors
                    Q_ = @(i) [0:(n(i)-1)]./n(i); [Y{1:3}]=ndgrid(Q_(1),Q_(2),Q_(3)); k=reshape(cat(3+1,Y{:}),[],3).';
                    % create structure
                    fbz = define(pc,n,k);
                otherwise
                    error('unknown n');
            end
        end

        function [ibz,i2f,f2i]    = get_ibz(fbz,pc,flag)

            import am_lib.* 

            if nargin<3; flag=''; end

            % get point symmetries [real-frac --> rec-frac] by transposing R
            [~,~,~,R] = pc.get_symmetries(); R = permute(R,[2,1,3]);

            % add inversion for time-reversal?
            if contains(flag,'addinv') 
            if ~any(all(all(R==-eye(3),1),2))
                R = complete_group( cat(3,R,-eye(3)) );
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
            fbz = define(pc,fbz.n,fbz.k(:,i2f),w);
            
            % pass along additional parameters if they exist
            for f = {'n','hv','F','F2'}
                if isfield(fbz,f) && ~isempty(fbz.(f{:}))
                    ibz.(f{:}) = fbz.(f{:}); 
                end
            end
            
            % if requested
            if contains(flag,'tetra')
                % get all tetrahedron
                tet = get_tetrahedra(fbz.recbas,fbz.n); 
                vol = abs(det([fbz.k(:,tet(:,1));[1,1,1,1]]))/factorial(3);
                % get irreducible tetrahedra
                [tet,~,tet_f2i] = unique(sort(f2i(tet)).','rows'); 
                tet=tet.'; tetw = hist(tet_f2i,[1:size(tet,2)].'-.5);
                % augment structure with tetrahedron
                ibz.ntets = size(tet,2);
                ibz.tet   = tet; 
                ibz.tetw  = tetw;
                ibz.tetv  = vol;
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

        function [bzp]            = get_bz_path(pc,n,brav)

            import am_lib.* am_dft.*

            % define kpoint path
            if     contains( lower(brav), 'fcc-short' )
                G=[0;0;0];  X1=[0;1;1]/2; X2=[2;1;1]/2;
                L=[1;1;1]/2; K=[6;3;3]/8;
                ql={'G','X','K','G','L'};
                qs=[G,X2,K,G];
                qe=[X1,K,G,L];
            elseif contains( lower(brav), 'fcc-long' )
                G=[0;0;0];  X1=[0;1;1]/2; W=[1;3;2]/4;
                U=[2;5;5]/8; L=[1;1;1]/2; K=[3;6;3]/8;
                ql={'G','X','W','K','G','L','U','W','L','K'};
                qs=[G,X1,W,K,G,L,U,W,L];
                qe=[X1,W,K,G,L,U,W,L,K];
            elseif contains( lower(brav), 'tetra' )
                G=[0;0;0];   Z=[0;0;1]/2; A=[1;1;1]/2; 
                M=[1;1;0]/2; X=[0;1;0]/2; R=[0;1;1]/2;
                ql={'G','X','M','G','Z','R','A','Z'};
                qs=[G,X,M,G,Z,R,A];
                qe=[X,M,G,Z,R,A,Z];
            elseif contains( lower(brav), 'sc' )
                G=[0;0;0];   X=[0;1;0]/2;
                M=[1;1;0]/2; R=[1;1;1]/2;
                ql={'G','X','M','G','R','X'};
                qs=[G,X,M,G,R];
                qe=[X,M,G,R,X];
            elseif contains( lower(brav), 'hex' )
                % for a pc.bas ordered like so:
                %     3.0531   -1.5266         0
                %          0    2.6441         0
                %          0         0    3.4526
                G=[0;0;0];   K=[1/3;1/3;0];   M=[1/2;0;0];
                A=[0;0;1/2]; H=[1/3;1/3;1/2]; L=[1/2;0;1/2];
                ql={'G','K','M','G','A','H','L','A'};
                qs=[G,K,M,G,A,H,L];
                qe=[K,M,G,A,H,L,A];
            else
                error('invalid bravais lattice');
            end

            % get number of kpoints
            nqs=size(qs,2); recbas = inv(pc.bas).';

            % get path: convert to [cart-recp] to get x spacings right then convert back to [frac-recp]
            [k,x,qt] = get_path(recbas*qs,recbas*qe,nqs,n); k=recbas\k;

            % create path object
            bzp_ = @(recbas,ql,qt,nks,x,k) struct('units','frac', ...
                'recbas',recbas,'ql',{{ql{:}}},'qt',qt,'nks',nks,'x',x,'k',k);
            bzp = bzp_(recbas,ql,qt,size(k,2),x,k);

            function [k,x,qt] = get_path(qs,qe,nqs,N)
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
                qt = x([1,N*[1:nqs]]);
            end
        end

        function [bzs]            = get_bz_surf(pc,n,vy,vx)
            % bzs=get_bz_surf(pc,[101,101],[1;0;0],[0;1;0]);
            % bzs=get_bvk_dispersion(bvk,bzs);
            % plot_bz_surf(bzs,1)

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

        function [bzl]            = get_bz_line(pc,n,vs,ve)
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

        function [bza]            = get_bz_angles(pc,hv,th2)

            import am_lib.* am_dft.*

            % get number of kpoints
            nks=numel(th2); recbas = inv(pc.bas).';

            % convert 2-theta into k points
            [kx,kz]=angle2kxkz(th2/2,th2,hv); k = [kx(:).';zeros(1,nks);kz(:).']; 

            % build path in recp.-frac
            k = recbas\k;

            % create path object
            bza_ = @(recbas,nks,x,k) struct('units','frac',...
                'recbas',recbas,'nks',nks,'n',[],'x',x,'k',k);
            bza = bza_(recbas,nks,th2(:).',k);
        end

        function [fbs,ibs]        = get_bz_bragg(uc,k_max,hv,threshold)
            % bragg brillouin zone
            import am_lib.* am_dft.*

            % k max is the largest wavevector magntiude considered
            if nargin < 2 || isempty(k_max); th2_max=180; lambda=0.15406; k_max = 2*sind(th2_max/2)/lambda; end
            if nargin < 3 || isempty(hv); hv = get_atomic_emission_line_energy(get_atomic_number('Cu'),'kalpha1'); end
            if nargin < 4 || isempty(threshold); threshold = am_lib.tiny; end

            fbs = get_fbs(uc,k_max);
            
            % get structure factors and scattering intensity
                [fbs.F,fbs.L,fbs.P] = get_structure_factor(uc,fbs,hv);
                fbs.I = abs(fbs.F).^2.*fbs.L.*fbs.P.*fbs.w; fbs.I = fbs.I ./ max(fbs.I(:))*100;
                % exlcude stuff above threshold
                ex_ = fbs.I > threshold;
                [fbs.F,fbs.L,fbs.P,fbs.I,fbs.w,fbs.k,fbs.nks] = ...
                    deal(fbs.F(ex_),fbs.L(ex_),fbs.P(ex_),fbs.I(ex_),fbs.w(ex_),fbs.k(:,ex_),sum(ex_));
            
            % get symmetrically equivalent bragg spots
            % save "irreducible bragg spots" structure
            [ibs,i2b,b2i] = get_ibz(fbs,uc,'nomod,addinv');
            ibs.i2b = i2b; ibs.b2i = b2i;
            fbs.i2b = i2b; fbs.b2i = b2i; 
        end
        
        function [fbs]            = get_fbs(uc,k_max,N)
            
            import am_lib.* 
            
            if nargin < 2 || isempty(k_max); th2_max=180; lambda=0.15406; k_max = 2*sind(th2_max/2)/lambda; end
            if nargin < 3 || isempty(N); N = ceil((max(k_max./(1./am_lib.normc_(uc.bas))))); end 
            
            % get miller indicies [hkl] excluding gamma
            k=permn_([N:-1:-N],3).'; k=k(:,~all(k==0,1)); 
            % get reciprocal basis vectors [1/nm]
            recbas = inv(uc.bas).';
            % sort by distance
            k = k(:,rankc_(normc_(recbas*k)));               
            % identify values which cannot be reached by diffractometer
            k = k(:,normc_(recbas*k)<k_max); 
            % save "full bragg spots" structure 
            bb_ = @(recbas,k) struct('units','frac-recp',...
                'recbas',recbas,'nks',size(k,2),'k',k,'w',ones(1,size(k,2)));
            fbs = bb_(recbas,k);
            
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
            axs_ = @(h,qt,ql) set(h,'Box','on','XTick',qt,'Xticklabel',ql);
            fig_(gcf);

            % plot results
            plot(bzp.x,fftinterp_(x,bzp.k,fbz.n), varargin{:});
            axs_(gca,bzp.qt,bzp.ql); axis tight; xlabel('Wavevector k');
        end

        function                    plot_dispersion_projected(dft,bzp,pc,flag)
            % FPOSCAR = 'POSCAR'; Ef = 5.0740;
            % [~,pc] = load_cells(FPOSCAR);
            % [dft]   = load_procar('evk/PROCAR',Ef);
            % [bzp]   = get_bz_path(pc,40,'sc');

            % get eigenvalues and band character weights
            c = zeros(dft.nbands,dft.nks); w = [];
            
            if     contains(flag,'orbital') % projected on orbitals
                character_labels = {'s','p','d','f'}; c = projections2color_(dft,pc,'orbital');
                if contains(flag,'atom')
                    if numel(pc.i2p)~=2; error('cannot broaden width in more than 2 dimensions'); end
                    w = projections2color_(dft,pc,'atom')/2;
                end
            elseif contains(flag,'atom') % projected on atoms
                character_labels = pc.symb; c = projections2color_(dft,pc,'atom');
            else % just plot bands straight up
                flag = [flag,'none'];
                character_labels=''; c = ones(dft.nbands,dft.nks); 
            end

            % define figure properties
            fig_ = @(h)       set(h,'color','white');
            axs_ = @(h,qt,ql) set(h,'Box','on','XTick',qt,'Xticklabel',ql);
            fig_(gcf);

            % plot band structure
            fig_(gcf); 
            if contains(flag,'none')
                plot(bzp.x,dft.E,'-','color','b');
            else
                hold on;
                for m = 1:dft.nbands
                    if isempty(w)
                        am_lib.plotc_(bzp.x,dft.E(m,:),c(m,:));
                    else
                        am_lib.plotc_(bzp.x,dft.E(m,:),c(m,:),w(m,:));
                        plot(bzp.x,dft.E,'-w','linewidth',0.1);
                    end
                end
                hold off;
                % apply color map and label axes
                colormap( am_lib.colormap_('spectral',100).^(2) ); h = colorbar; caxis([0,1]);
                cticks = am_lib.assign_cmap_(eye(numel(character_labels))); [~,inds] = sort( cticks );
                set(h,'Ticks',cticks(inds),'TickLabels',character_labels(inds));
            end
            
            % label axes
            axs_(gca,bzp.qt,bzp.ql); axis tight; ylabel('Energy [eV]'); xlabel('Wavevector k');

            % plot fermi level
            line([0,bzp.x(end)],[0,0],'linewidth',2,'color',[1,1,1]*0.5,'linestyle',':');
            
            
            
            function c = projections2color_(dft,pc,flag)
                % normalize columns of matrix
                normc_ = @(m) ones(size(m,1),1)*sqrt(1./sum(m.*m)).*m;
                % colorize
                if     contains(flag,'orbital')
                for i = 1:dft.nks
                    %    s     py     pz     px    dxy    dyz    dz2    dxz    dx2    f-3    f-2    f-1     f0     f1     f2     f3
                    % l-PROJECTION (sum over spin, atoms, m-quantum numbers)
                    % lmproj(nspins,norbitals,nions,nbands,nkpts)
                    Vp(1,:) = squeeze(sum(sum(sum(dft.lmproj(:,   [1],:,:,i),1),2),3)); % s
                    Vp(2,:) = squeeze(sum(sum(sum(dft.lmproj(:, [2:4],:,:,i),1),2),3)); % p
                    Vp(3,:) = squeeze(sum(sum(sum(dft.lmproj(:, [5:9],:,:,i),1),2),3)); % d
                    Vp(4,:) = squeeze(sum(sum(sum(dft.lmproj(:,[9:16],:,:,i),1),2),3)); % f
                    c(:,i) = am_lib.assign_cmap_(normc_(Vp));
                end
                elseif contains(flag,'atom')
                    for i = 1:dft.nks
                        % l-PROJECTION (sum over spin, atoms, m-quantum numbers)
                        % lmproj(nspins,norbitals,nions,nbands,nkpts)
                        for j = 1:numel(unique(pc.p2i))
                            Vp(j,:) = squeeze(sum(sum(sum(dft.lmproj(:,:,pc.p2i==j,:,i),1),2),3));
                        end
                        c(:,i) = am_lib.assign_cmap_(normc_(Vp));
                    end
                end
            
            end
        end

        function                    plot_nesting(ibz,fbz,bzp,degauss,Ep, varargin)

            import am_lib.* am_dft.*

            % plot results
            plot_interpolated(fbz,bzp, ibz2fbz(fbz,ibz,get_nesting(fbz,ibz,degauss,Ep)) , varargin{:})
        end
        
        function [Aibz]          =  plot_nesting_jdos(fbz,ibz,bzp,dft,degauss,Ep)
            
            Aibz = am_dft.get_nesting_jdos(fbz,ibz,dft,degauss,Ep);
            Afbz = am_dft.ibz2fbz(fbz,ibz,Aibz);
            Abzp = am_lib.fftinterp_(Afbz,bzp.k,fbz.n);
            
            subplot(3,1,1:2);
            [EE,XX]=ndgrid(Ep,bzp.x);
            % ex_ = ones(size(Abzp)); ex_(Abzp<0.01) = 0;
            % surf(XX,EE*scale_energies_,(abs(Abzp)./Ep.^2).^(0.1).*ex_,'edgecolor','none'); view([0 0 1]); axis tight;
            surf(XX,EE,abs(Abzp),'edgecolor','none'); view([0 0 1]); axis tight;

            axs_ = @(h,qt,ql) set(h,'Box','on','XTick',qt,'Xticklabel',ql);
            axs_(gca,bzp.qt,bzp.ql); axis tight; ylabel('Energy [eV]'); xlabel('Wavevector k');

        end

        function                   plot_bz(fbz)

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

        function                   plot_bz_path(bzp)

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

        function                   plot_bz_surf(bzs,band)
            % bzs=get_bz_surf(pc,[101,101],[1;0;0],[0;1;0]);
            % bzs=get_bvk_dispersion(bvk,bzs);
            % plot_bz_surf(bzs,1)

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

        function h               = plot_fermi_surface(ibz,fbz,dft,Ef,flag,C)
            if nargin < 5; flag=''; end
            if nargin < 6; C=[]; end % C is the colorcode
            
            % % color based on curvature? (each band is a scalar field: hessian is diagonal)
            % r = ndgrid(1:fbz.n(1),1:fbz.n(2),1:fbz.n(3)); r = permute(r,[4,1,2,3]);
            % for i = 1:dft.nbands; C(i,:,:,:) = ifftn( fftn(E(i,:,:,:)).*(1i*r).^2 ); end; C = abs(real(C));
            
            if contains(flag,'jdos')
                % checked against vasp. confirmed.
                % only consider transitions from the valence band (E<0) to the conduction band (E>0)
                Ev = permute(dft.E,[1,3,2]); Ec = permute(dft.E,[3,1,2]); E = Ec - Ev;
                % E(valence,conduction,kpoints) 
                % make all valence    - valence    states E => 1E8 (something large enough to take it out of range)
                % make all conduction - conduction states E => 1E8 (something large enough to take it out of range)
                % proper way would be to compute fermi functions and use that as the projection weight
                % but just doing it quick and dirty here (equivalent to assuming kt -> 0, fermi fucntion -> step function)
                ex_ = (Ec<0) & (Ev<0); E(ex_) = 1E8;
                ex_ = (Ec>0) & (Ev>0); E(ex_) = 1E8;
                E = reshape(E,dft.nbands^2,dft.nks); E = sort(E);
                E = reshape( am_dft.ibz2fbz(fbz,ibz,E) ,[dft.nbands.^2,fbz.n]); 
            else % just regular dos
                E = reshape( am_dft.ibz2fbz(fbz,ibz,dft.E) ,[dft.nbands,fbz.n]); 
            end


            % plot                                                                                 [4,4,1]  2
            h = am_lib.plot_isosurface_(reshape(fbz.recbas*fbz.k,[3,fbz.n]), fbz.recbas, E, Ef, C, [1,1,1], 3, 'cubic,center');
            % plot bz boundary
            if contains(flag,'tetra'); am_dft.draw_bz_boundary(fbz,'tetra'); end
            % set lighting settings
            daspect([1 1 1]); axis tight; axis off;
            campos([10.7135,23.0003,9.4299]); 
            camlight('right','infinite');
            camproj('perspective');
        end
        
        function h               = draw_bz_boundary(fbz,flag)
            switch flag
                case 'tetra'; h = draw_tetragonal_bz_boundary(fbz);
                otherwise; error('unknown geometry');
            end
            
            function draw_tetragonal_bz_boundary(fbz)
                R3 = [ 0,-1, 0;-1, 0, 0; 0, 0,-1];
                s = [1,1,1];
                xy = fbz.recbas*[[0;0;0],[s(1);0;0],[s(1);s(2);0],[0;s(2);0]];
                xz = fbz.recbas*[[0;0;0],[0;0;s(3)],[s(1);0;s(3)],[s(1);0;0]];
                alpha = 0.015;
                patch('Faces',[1:4],'Vertices',(xy)','EdgeColor','k','FaceColor','k','facealpha',alpha,'LineWidth',1);
                patch('Faces',[1:4],'Vertices',(xz)','EdgeColor','k','FaceColor','k','facealpha',alpha,'LineWidth',1);
                patch('Faces',[1:4],'Vertices',(xy+fbz.recbas*[0;0;s(3)])','EdgeColor','k','FaceColor','k','facealpha',alpha,'LineWidth',1);
                patch('Faces',[1:4],'Vertices',(xz+fbz.recbas*[0;s(2);0])','EdgeColor','k','FaceColor','k','facealpha',alpha,'LineWidth',1);
                patch('Faces',[1:4],'Vertices',(R3*xz+fbz.recbas*[0;s(2);s(3)])','EdgeColor','k','FaceColor','k','facealpha',alpha,'LineWidth',1);
                patch('Faces',[1:4],'Vertices',(R3*xz+fbz.recbas*[s(1);s(2);s(3)])','EdgeColor','k','FaceColor','k','facealpha',alpha,'LineWidth',1);
            end
        end
        
    end
    
    % aux library

    methods (Static)
    end
    
end
