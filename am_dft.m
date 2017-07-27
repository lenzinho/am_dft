classdef am_dft

    
    properties (Constant)
        tiny      = 1E-3; % precision of atomic coordinates
        eps       = 1E-8; % numerical precision
        units_eV  = 0.06465555; % sqrt( [eV/Ang^2] * [1/amu] ) --> 0.06465555 [eV]
        units_THz = 98.22906;   % sqrt( [eV/Ang^2] * [1/amu] ) --> 98.22906 [THz=1/ps]
        units_GHz = 98229.06;   % sqrt( [eV/Ang^2] * [1/amu] ) --> 98229.06 [GHz=1/fs]
        
        % mex compiler parameters
        usemex    = false;
        FC        = 'ifort'; 
        FFLAGS    = '-O3 -parallel -fpp -fPIC -lmx -lmex -lmat -nofor_main -bundle -implicitnone -assume realloc_lhs';
        MPATH     = '/Applications/MATLAB_R2016b.app';
        LIBS      = ['-L',am_dft.MPATH,'/bin/maci64 -I',am_dft.MPATH,'/extern/include'];
        EXT       = '.mexmaci64';
        DEBUG     = '-debug'
    end
    
    % program level
    
    methods (Static)
        
        function [uc,pc,md,bvk,pp,bvt,pt] = get_phonons(opts)
            %
            % clear; clc;
            %
            % opts.continue=true;
            % opts.fposcar='POSCAR';
            % opts.fforce_position='infile.force_position';
            % opts.cutoff2=3;
            % opts.cutoff3=3; % ignored if set to zero
            % opts.dt=2; %fs 
            %
            % [uc,pc,md,bvk,pp,bvt,pt] = get_phonons(opts)
            %
            
            import am_lib.*
            import am_dft.*

            % check inputs
            for f = {opts.fposcar,opts.fforce_position}
            if ~exist(f{:},'file'); error('File not found: %s',f{:}); end
            end

            % file name to save to/load from
            sname = sprintf('%s_%s_%0.2f_%0.2f.mat',...
                opts.fposcar, opts.fforce_position, opts.cutoff2, opts.cutoff3);
            if and(exist(sname,'file'),opts.continue); load(sname); else

                % get cells
                [uc,pc] = get_cells(opts.fposcar); 
                % load md
                [md] = load_md(uc,opts.fforce_position,opts.dt);
                % get both pairs and triplets?
                if or(opts.cutoff3==0,isempty(opts.cutoff3))
                    % get pair shells
                    [bvk,pp] = get_bvk(pc,uc,md,opts.cutoff2);
                    % get triplet shells
                    bvt = []; pt = [];
                else
                    % get pair shells ( do not fit for pair force constants here )
                    [bvk,pp] = get_bvk(pc,uc,md,opts.cutoff2,'-identify -model');
                    % get triplet shells ( fit for pair and triplets simultaneously here )
                    [bvt,pt,bvk] = get_bvt(pc,uc,md,opts.cutoff3,bvk,pp,'-identify -model -fit');
                end
                % save results
                save(sname,'uc','pc','md','bvk','pp','bvt','pt');
            end
            
            % plot results
                % plot correlation for dft vs bvk forces on atoms
                figure('color','white'); plot_bvk_vs_aimd(uc,md,bvk,pp,bvt,pt); drawnow;
        end

        function [uc,pc,dft,tb,pp] = get_tightbinding(opts)
            % clear;clc
            % 
            % import am_lib.*
            import am_dft.*
            % 
            % opts.continue=true;
            % opts.fposcar='POSCAR';
            % opts.feigenval='EIGENVAL.scf';
            % opts.spdf={'p','d'};
            % opts.nskips=5; 
            % opts.Ef=12.8174; % grep OUTCAR for E-fermi  
            % opts.cutoff2=3;
            % 
            % [uc,pc,dft,tb,pp] = get_tightbinding(opts);
            %
            % NOTE: correlation coefficient for a good fit should exceed
            % 0.99. For VN with pd orbitals and a pair cutoff of 3 Ang
            % (second neighbor), R^2 = 0.998.

            import am_lib.*
            import am_dft.*

            % check inputs
            for f = {opts.fposcar}
            if ~exist(f{:},'file'); error('File not found: %s',f{:}); end
            end

            % file name to save to/load from
            sname = sprintf('%s_%s_%02i_%0.2f_%0.2f_%s.mat',...
                opts.fposcar, opts.feigenval, opts.nskips, ...
                opts.Ef, opts.cutoff2, strjoin(opts.spdf,'+') );
            if and(exist(sname,'file'),opts.continue); load(sname); else

                % get cells
                [uc,pc] = get_cells(opts.fposcar); 
                % load dft
                [dft]   = load_eigenval(opts.feigenval,opts.Ef);
                % get tb
                [tb,pp] = get_tb(pc,get_supercell(pc,diag([5,5,5])),dft,opts.cutoff2,opts.spdf,opts.nskips);
                % save results
                save(sname,'uc','pc','dft','tb','pp');
                
            end
            
            % plot results
                % plot correlation for dft vs bvk forces on atoms
                figure('color','white'); plot_tb_vs_dft(tb,dft); drawnow;
        end
     
    end
    
    % core

    methods (Static)

        % vasp

        function           save_poscar(uc,fposcar)
            n = size(uc.tau,3);
            for i = 1:n
                if n == 1; fid=fopen(sprintf('%s'     ,fposcar,i),'w'); else
                           fid=fopen(sprintf('%s_%06i',fposcar,i),'w');
                end
                fprintf(fid,' %i ',uc.nspecies); 
                fprintf(fid,'%s \n',sprintf('POSCAR %i of %i',i,n)); 
                fprintf(fid,'%12.8f \n',1.0);  % latpar
                fprintf(fid,'%12.8f %12.8f %12.8f \n',uc.bas(:,1)); 
                fprintf(fid,'%12.8f %12.8f %12.8f \n',uc.bas(:,2)); 
                fprintf(fid,'%12.8f %12.8f %12.8f \n',uc.bas(:,3)); 
                fprintf(fid,' %s ',uc.symb{:}); fprintf(fid,'\n');
                fprintf(fid,' %i ',uc.nspecies); fprintf(fid,'\n');
                fprintf(fid,'Direct \n'); 
                fprintf(fid,'%12.8f %12.8f %12.8f \n',uc.tau(:,:,i));
                fclose(fid);
            end
        end

        function [uc]    = load_poscar(fposcar)
            import am_lib.*
            import am_dft.*
            fid=fopen(fposcar,'r');            % open file
            header=fgetl(fid); uc.units='frac';% read header (often times contains atomic symbols)
            latpar=sscanf(fgetl(fid),'%f');    % read lattice parameter
            a1=sscanf(fgetl(fid),'%f %f %f');  % first basis vector
            a2=sscanf(fgetl(fid),'%f %f %f');  % second basis vector
            a3=sscanf(fgetl(fid),'%f %f %f');  % third basis vector
            uc.bas=latpar*[a1,a2,a3];          % construct the basis (column vectors)
            buffer=fgetl(fid);                 % check vasp format
            if ~isempty(sscanf(buffer,'%f'))
                % vasp 4 format (symbols are missing, jumps straight to species)
                uc.symb=regexp(header, '([^ \s][^\s]*)', 'match');
                uc.nspecies=sscanf(buffer,repmat('%f' ,1,length(uc.symb)))';
            else
                % vasp 5 format (symbols are present)
                uc.symb=regexp(buffer, '([^ \s][^\s]*)', 'match');
                uc.nspecies=sscanf(fgetl(fid),repmat('%f' ,1,length(uc.symb)))';
            end
            for i=1:length(uc.nspecies); uc.mass(i)=Z2mass(symb2Z(uc.symb{i})); end
            uc.natoms=sum(uc.nspecies);
            coordtype=lower(strtrim(fgetl(fid)));
            l=0;
            for i=1:length(uc.nspecies)
                for j=1:uc.nspecies(i); l=l+1;
                    uc.tau(:,l)=sscanf(fgetl(fid),'%f %f %f');
                    uc.species(l)=i;
                end
            end
            if ~strcmp(coordtype(1),'d'); uc.tau=uc.bas\uc.tau*latpar; end
            fclose(fid);
            %
            function [Z] = symb2Z(symb)
            s = {'h'  ,'he' ,'li' ,'be' ,'b'  ,'c'  ,'n'  ,'o'  ,'f'  ,'ne' ,'na' ,'mg' ,'al' ,'si' ,'p'  ,'s'  , ...
                 'cl' ,'ar' ,'k'  ,'ca' ,'sc' ,'ti' ,'v'  ,'cr' ,'mn' ,'fe' ,'co' ,'ni' ,'cu' ,'zn' ,'ga' ,'ge' , ...
                 'as' ,'se' ,'br' ,'kr' ,'rb' ,'sr' ,'y'  ,'zr' ,'nb' ,'mo' ,'tc' ,'ru' ,'rh' ,'pd' ,'ag' ,'cd' , ...
                 'in' ,'sn' ,'sb' ,'te' ,'i'  ,'xe' ,'cs' ,'ba' ,'la' ,'ce' ,'pr' ,'nd' ,'pm' ,'sm' ,'eu' ,'gd' , ...
                 'tb' ,'dy' ,'ho' ,'er' ,'tm' ,'yb' ,'lu' ,'hf' ,'ta' ,'w'  ,'re' ,'os' ,'ir' ,'pt' ,'au' ,'hg' , ...
                 'tl' ,'pb' ,'bi' ,'po' ,'at' ,'rn' ,'fr' ,'ra' ,'ac' ,'th' ,'pa' ,'u'  ,'np' ,'pu' ,'am' ,'cm' , ...
                 'bk' ,'cf' ,'es' ,'fm' ,'md' ,'no' ,'lr' ,'rf' ,'db' ,'sg' ,'bh' ,'hs' ,'mt' ,'ds' ,'rg' ,'uub', ...
                 'uut','uuq','uup','uuh'}; Z = find(strcmp(strtrim(lower(symb)),s));
            end
            function [mass] = Z2mass(Z)
            m = [   1.007947000,     4.002602000,     6.941200000,     9.012182000,    10.811500000, ...
                   12.011100000,    14.006747000,    15.999430000,    18.998403000,    20.179760000, ...
                   22.989769000,    24.305060000,    26.981540000,    28.085530000,    30.973762000, ...
                   32.066600000,    35.452790000,    39.948100000,    39.098310000,    40.078900000, ...
                   44.955911000,    47.883000000,    50.941510000,    51.996160000,    54.938051000, ...
                   55.847300000,    58.933201000,    58.693400000,    63.546300000,    65.392000000, ...
                   69.723100000,    72.612000000,    74.921592000,    78.963000000,    79.904000000, ...
                   83.801000000,    85.467830000,    87.621000000,    88.905852000,    91.224200000, ...
                   92.906382000,    95.941000000,    98.000000000,   101.072000000,   102.905503000, ...
                  106.421000000,   107.868220000,   112.411800000,   114.821000000,   118.710700000, ...
                  121.757000000,   127.603000000,   126.904473000,   131.292000000,   132.905435000, ...
                  137.327700000,   138.905520000,   140.115400000,   140.907653000,   144.243000000, ...
                  145.000000000,   150.363000000,   151.965900000,   157.253000000,   158.925343000, ...
                  162.503000000,   164.930323000,   167.263000000,   168.934213000,   173.043000000, ...
                  174.967100000,   178.492000000,   180.947910000,   183.853000000,   186.207100000, ...
                  190.210000000,   192.223000000,   195.083000000,   196.966543000,   200.593000000, ...
                  204.383320000,   207.210000000,   208.980373000,   209.000000000,   210.000000000, ...
                  222.000000000,   223.000000000,   226.025000000,   227.028000000,   232.038110000, ...
                  231.035900000,   238.028910000,   237.048000000,   244.000000000,   243.000000000, ...
                  247.000000000,   247.000000000,   251.000000000,   252.000000000,   257.000000000, ...
                  258.000000000,   259.000000000,   262.000000000,   261.000000000,   262.000000000, ...
                  263.000000000,   262.000000000,   265.000000000,   266.000000000]; mass = m(Z);
            end
        end

        function [dft]   = load_eigenval(feigenval,Ef)
            % load dispersion [frac-recp] and shift Fermi energy to zero
            
            fprintf(' ... loading dft band structure from %s ',feigenval); tic;
            
            fid=fopen(feigenval);
                % skip first five lines
                for i = 1:5; fgetl(fid); end
                buffer = strsplit(strtrim(fgetl(fid)));
                dft.nelecs = sscanf(buffer{1},'%i');
                dft.nks    = sscanf(buffer{2},'%i');
                dft.nbands = sscanf(buffer{3},'%i');
                for i = 1:dft.nks
                    % skip line
                    fgetl(fid);
                    % get kpnts
                    buffer = strsplit(strtrim(fgetl(fid)));
                    dft.k(1,i) = sscanf(buffer{1},'%f');
                    dft.k(2,i) = sscanf(buffer{2},'%f');
                    dft.k(3,i) = sscanf(buffer{3},'%f');
                    % loop over bands
                    for j = 1:dft.nbands
                        buffer = strsplit(strtrim(fgetl(fid)));
                        dft.E(j,i)  = sscanf(buffer{2},'%f');
                    end
                    dft.E(:,i) = sort(dft.E(:,i));
                end
                % zero fermi and sort bands in increasing order
                dft.E = sort(dft.E - Ef); 
                % close file
            fclose(fid);
            
            fprintf('(%.f secs)\n',toc);
        end

        function [md]    = load_md(uc,fforces,dt)
            %
            % Loads outcar preprocessed with outcar2fp.sh (generate_script).
            %

            import am_lib.*
            import am_dft.*
            
            fprintf(' ... loading displacements vs forces'); tic;

            % count number of lines in file and check that all runs completed properly
            nlines = count_lines_(fforces); if mod(nlines,uc.natoms)~=0; error('lines appear to be missing.'); end

            % open file and parse: use single precision here, solves for force constants much faster
            fid = fopen(fforces); fd = reshape(single(fscanf(fid,'%f')),6,uc.natoms,nlines/uc.natoms); fclose(fid);
            
            % convert to [uc-frac]
            fd(1:3,:,:) = matmul_(inv(uc.bas),fd(1:3,:,:));
            fd(4:6,:,:) = matmul_(inv(uc.bas),fd(4:6,:,:));

            fd_ = @(uc,force,tau,vel,dt) struct('units','frac', ...
                'bas',uc.bas,'bas2pc',uc.bas2pc,'tau2pc',uc.tau2pc,...
                'symb',{{uc.symb{:}}},'mass',uc.mass,'nspecies',uc.nspecies, ...
                'natoms',uc.natoms,'force',force,'tau',tau,'vel',vel,'species',uc.species, ...
                'dt',dt,'nsteps',size(tau,3));
            md = fd_(uc,fd(4:6,:,:),fd(1:3,:,:),cat(3,zeros(3,uc.natoms),(mod_(diff(fd(1:3,:,:),1,3)+.5)-.5)/dt),dt);
            
            % match to uc for saftey
            md = match_cell(md,uc);
            
            fprintf(' (%.f secs)\n',toc);
        end

        function [en]    = load_band_energies(nbands,fbands)
            %
            % Loads outcar preprocessed with outcar2en.sh (generate_script).
            %
            
            import am_lib.*
            import am_dft.*
            
            % count number of lines in file and check that all runs completed properly
            nlines = count_lines_(fbands); if mod(nlines,nbands)~=0; error('lines appear to be missing.'); end

            % open file and parse
            nsteps=nlines/nbands; fid=fopen(fbands); en=reshape(fscanf(fid,'%f'),nbands,nsteps); fclose(fid);
        end
        
        
        % symmetry

        function [T,H,S,R]    = get_symmetries(pc)
            % T = all possible translations which restore the crystal to iteself
            % H = holohogries (all possible rotations which restore the bravais lattice onto iteself)
            % S = space group symmetries
            % R = point group symmetries
            
            import am_lib.*
            import am_dft.*
            
            % define function to check first two dimensions
            check3_ = @(A) all(all(abs(A)<am_lib.tiny,1),2);
            
            % define function to sort atoms and species into a unique order (reference)
            X_ = @(tau,species) sortc_([species;mod_(tau)]); X = X_(pc.tau(:,:,1),pc.species);

            % get vectors that preserve periodic boundary conditions
            N=1; T=mod_(pc.tau(:,pc.species==pc.species(N))-pc.tau(:,N)); nTs=size(T,2); T_ck=false(1,nTs);
            for j = 1:nTs; T_ck(j) = check3_( X_(pc.tau(1:3,:,1)-T(:,j),pc.species)-X ); end
            T=[T(:,T_ck),eye(3)]; T=T(:,rankc_(normc_(T))); nTs = size(T,2);

            if nargout == 1; return; end
            
            % get arithmetic holodries (symmetries for which R'*g*R = g; g = bas'*bas)
            N=9; Q=[-1:1]; nQs=numel(Q);[Y{N:-1:1}]=ndgrid(1:nQs); L=reshape(Q(reshape(cat(N+1,Y{:}),[],N)).',3,3,[]);
            get_holodries_frac_ = @(M) L(:,:,check3_(matmul_(matmul_(permute(L,[2,1,3]),M'*M),L)-M'*M));
            H = get_holodries_frac_(pc.bas); nHs = size(H,3);
            id = member_(flatten_(eye(3)),reshape(H,3^2,[])); H(:,:,[1,id])=H(:,:,[id,1]);

            if nargout == 2; return; end
            
            % get seitz operators which leave the atomic basis invariant
            S = zeros(4,4,nHs*nTs); S(4,4,:)=1; nSs=0;
            for i = 1:nHs; for j = 1:nTs
                if check3_( X_(H(:,:,i)*pc.tau+T(:,j),pc.species) - X ); nSs=nSs+1; S(1:3,1:4,nSs)=[ H(:,:,i), T(:,j) ]; end
            end; end; S = S(:,:,1:nSs); 
        
            % set identity first
            id = member_(flatten_(eye(4)),reshape(S,4^2,[])); S(:,:,[1,id])=S(:,:,[id,1]);
            
            if nargout == 3; return; end

            % get point symmetries
            R  = reshape(uniquecol_( reshape(S(1:3,1:3,:),[9,nSs]) ),3,3,[]);
            
            % set identity first
            id = member_(flatten_(eye(3)),reshape(R,3^2,[])); R(:,:,[1,id])=R(:,:,[id,1]);
        end

        function [MT,E,I]     = get_multiplication_table(S)
            % get multiplication table: S(:,:,i)*S(:,:,j) = S(:,:,MT(i,j)

            import am_lib.*
            import am_dft.*
            
            if     size(S,1) == 4
                % seitz operator (applies mod to translational components)
                md_ = @(X) [X(1:12,:);mod_(X(13:15,:));X(16:end,:)];
                rs_ = @(X) md_(reshape(X,4^2,[]));
                nSs = size(S,3);

                MT = reshape( member_( rs_(matmulp_(S,permute(S,[1,2,4,3]))) , rs_(S) ) , nSs,nSs);
                
            elseif size(S,1) == 3 
                % point operator
                rs_ = @(X) reshape(X,3^2,[]);
                nSs = size(S,3);

                MT = reshape( member_( rs_(matmulp_(S,permute(S,[1,2,4,3]))) , rs_(S) ) , nSs,nSs);
                
            elseif size(S,1) == 1
                % seitz operator combined with permutation (represented as a two-part cell)
                if     size(S{1},1) == 4
                    % seitz operator (applies mod to translational components)
                    md_ = @(X) [X(1:12,:);mod_(X(13:15,:));X(16:end,:)];
                    rs_ = @(X) md_(reshape(X,4^2,[]));
                    nSs = size(S{1},3);
                    
                    ref = [rs_(S{1});reshape(S{2},size(S{2},1),[])];
                    opr = [rs_(matmulp_(S{1},permute(S{1},[1,2,4,3]))); reshape(operm_(S{2},S{2}),size(S{2},1),[])];

                    MT = reshape( member_( opr , ref ) , nSs,nSs);
                elseif size(S{1},1) == 3
                    % point operator
                    rs_ = @(X) reshape(X,3^2,[]);
                    nSs = size(S{1},3);
                    
                    ref = [rs_(S{1});reshape(S{2},size(S{2},1),[])];
                    opr = [rs_(matmulp_(S{1},permute(S{1},[1,2,4,3]))); reshape(operm_(S{2},S{2}),size(S{2},1),[])];

                    MT = reshape( member_( opr , ref ) , nSs,nSs);
                end
            end
            
            if any(MT(:)==0)
                error('MT is incorrect. Check for mistakes in the symmetry and ensure that the symmetry is in the primitive basis');
            end
            
            % find identity 
            if nargout>1; E = find(all(MT==[1:nSs].',1)); end

            % get inverse indicies
            if nargout>2; I = [MT==E]*[1:nSs].'; end
        end

        function [CT,c_id,irr]= get_irreps(S)
            % s2c = identifies the class to which symmetries belong
            % CT = character table
            % irreps
            
            import am_lib.*
            import am_dft.*
            
            % get regular rep G by putting identity along diagonal of multiplciation table
            [MT,~,I] = get_multiplication_table(S); nGs = size(MT,2);
            G = double(accessc_(MT,I)==permute([1:nGs],[1,3,2])); 

            % initialize decomposition loop
            U = eye(nGs); inds = ones(nGs,1); ninds = 1;

            % loop until irreps are fully decomposed
            while true
                % loop over cycle structures
                for j = 1:max(inds)
                    ex_ = inds==j;

                    H = dixon_decomposition_( G(ex_,ex_,:) );
                    [Vp,E] = eig(H,'vector'); [Vp] = orth_(Vp,E);
                    for ig = 1:nGs
                        G(ex_,ex_,ig) = Vp\G(ex_,ex_,ig)*Vp;
                    end
                    U(:,ex_) = (U(:,ex_)*Vp);
                end

                inds = [1:nGs]*merge_(double(sum(abs(G),3)>am_lib.eps)); 
                if ninds == max(inds); break; else
                    ninds = max(inds);
                end
            end

            % get character table
            CT = zeros(nGs,max(inds));
            for i = 1:max(inds); for j = 1:nGs
                CT(j,i) = trace(G(inds==i,inds==i,j));
            end; end

            % get irreducible irreps
            [CT,ir] = unique(round(CT).','rows'); CT=CT.';

            % get irreducible classes
            [CT,~,c_id] = unique(round(CT),'rows','stable'); CT=CT.'; c_id=c_id(:).';

            % get irreducible representations
            nirreps = numel(ir); irr = cell(1,nirreps);
            for i = 1:nirreps
                irr{i} = G(inds==ir(i),inds==ir(i),1:nGs); 
                % this makes them look nice but multiplcation table is not preserved
                % for j = 1:nGs
                %     irrep{i}(:,:,j) = diag(sort(eig(irrep{i}(:,:,j))));
                % end
            end

            function H = dixon_decomposition_(rr)
                nbases=size(rr,1);
                nsyms =size(rr,3);

                for r = 1:nbases
                for s = 1:nbases
                    Hrs = zeros(nbases);
                    if     r==s
                        Hrs(r,s) = 1;
                    elseif r>s
                        Hrs(r,s) = 1;
                        Hrs(s,r) = 1;
                    elseif r<s
                        Hrs(r,s) = sqrt(-1);
                        Hrs(s,r) =-sqrt(-1);
                    end

                    H(1:nbases,1:nbases) = 0;
                    for q = 1:nsyms
                        H = H + rr(:,:,q)' * Hrs * rr(:,:,q);
                    end
                    H = H / nsyms;

                    if any(abs( H(1,1)*eye(nbases)-H ) >am_lib.eps); return; end
                end
                end
                H = eye(nbases);
            end
        end
        
        function c_id         = identify_classes(MT)
            %
            % for AX = XB, if elements A and B are conjugate pairs for some other element X in the group,  they are in the same class
            %
            import am_lib.*
            import am_dft.*

            nsyms = size(MT,2);
            % allocate space for conjugacy class
            c_id = zeros(nsyms,1);
            % allocate space for conjugate elements
            conjugates = zeros(nsyms,1);
            % get inverse indicies
            I = (MT==1) * [1:size(MT,1)].';
            % determine conjugacy classes
            k = 0;
            for i = 1:nsyms
            if c_id(i)==0
                k=k+1;
                % conjugate each element with all other group elements
                % A = X(j) * B * X(j)^-1
                for j = 1:nsyms
                    conjugates(j) = MT(j,MT(i,I(j)));
                end
                % for each subgroup element created by conjugation find the corresponding index of the element in the group
                % in order to save the class class_id number
                for j = 1:nsyms
                    c_id( conjugates(j) ) = k;
                end
            end
            end
            % relabel classes based on how many elements each class has
            c_id = reindex_using_occurances(c_id);
        end
        
        function ps_id        = identify_point_symmetries(R)
            nsyms=size(R,3); ps_id = zeros(1,nsyms);
            for i = 1:nsyms
                % get trace and determinant (fractional)
                tr = trace(R(1:3,1:3,i)); dt = det(R(1:3,1:3,i));
                if     and(tr==+3,dt==+1); ps_id(i) = 1;  % 'e'
                elseif and(tr==-1,dt==+1); ps_id(i) = 2;  % 'c_2'
                elseif and(tr==+0,dt==+1); ps_id(i) = 3;  % 'c_3'
                elseif and(tr==+1,dt==+1); ps_id(i) = 4;  % 'c_4'
                elseif and(tr==+2,dt==+1); ps_id(i) = 5;  % 'c_6'
                elseif and(tr==-3,dt==-1); ps_id(i) = 6;  % 'i'
                elseif and(tr==+1,dt==-1); ps_id(i) = 7;  % 's_2'
                elseif and(tr==+0,dt==-1); ps_id(i) = 8;  % 's_6'
                elseif and(tr==-1,dt==-1); ps_id(i) = 9;  % 's_4'
                elseif and(tr==-2,dt==-1); ps_id(i) = 10; % 's_3'
                else;                      ps_id(i) = 0;  % unknown
                end
            end
        end
        
        function pg_code      = identify_pointgroup(R)
            % 
            % Point symmetries in fractional coordinates so that they are nice integers which can be easily classified.
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
            import am_lib.*
            import am_dft.*
            %
            nsyms = size(R,3);
            % identify point symmetries
            ps_id = identify_point_symmetries(R);
            % count each type of symmetry
            nc2 = sum(ps_id==2);  % 'c_2'
            nc3 = sum(ps_id==3);  % 'c_3'
            nc4 = sum(ps_id==4);  % 'c_4'
            nc6 = sum(ps_id==5);  % 'c_6'
            ni  = sum(ps_id==6);  % 'i' 
            ns2 = sum(ps_id==7);  % 's_2'
            ns4 = sum(ps_id==9);  % 's_4'
            % identify point group by comparing number and types of symmetries
            if         nsyms==1         ; pg_code=1 ;
            elseif     nsyms==48        ; pg_code=32;
            elseif     nsyms==16        ; pg_code=20;
            elseif     nsyms==3         ; pg_code=9 ;
            elseif and(nsyms==2 , ni==1); pg_code=2 ;
            elseif and(nsyms==2 ,nc2==1); pg_code=3 ;
            elseif and(nsyms==2 ,ns2==1); pg_code=4 ;
            elseif and(nsyms==4 , ni==1); pg_code=5 ;
            elseif and(nsyms==4 ,nc2==3); pg_code=6 ;
            elseif and(nsyms==4 ,ns2==2); pg_code=7 ;
            elseif and(nsyms==4 ,nc4==2); pg_code=14; % nc4 == 2 is correct, rather than nc4 == 1.
            elseif and(nsyms==4 ,ns4==2); pg_code=15;
            elseif and(nsyms==6 , ni==1); pg_code=10;
            elseif and(nsyms==6 ,nc2==3); pg_code=11;
            elseif and(nsyms==6 ,ns2==3); pg_code=12;
            elseif and(nsyms==6 ,nc2==1); pg_code=21;
            elseif and(nsyms==6 ,ns2==1); pg_code=22;
            elseif and(nsyms==8 ,ns2==3); pg_code=8 ;
            elseif and(nsyms==8 ,ns2==1); pg_code=16;
            elseif and(nsyms==8 ,ns2==0); pg_code=17;
            elseif and(nsyms==8 ,ns2==4); pg_code=18;
            elseif and(nsyms==8 ,ns2==2); pg_code=19;
            elseif and(nsyms==12,ns2==3); pg_code=13;
            elseif and(nsyms==12,ns2==1); pg_code=23;
            elseif and(nsyms==12,nc2==7); pg_code=24;
            elseif and(nsyms==12,ns2==6); pg_code=25;
            elseif and(nsyms==12,ns2==4); pg_code=26;
            elseif and(nsyms==12,nc3==8); pg_code=28;
            elseif and(nsyms==24,nc6==2); pg_code=27;
            elseif and(nsyms==24, ni==1); pg_code=29;
            elseif and(nsyms==24,nc4==6); pg_code=30;
            elseif and(nsyms==24,ns4==6); pg_code=31;
            else;                         pg_code=0 ;
            end
        end
        
        function sg_code      = identify_spacegroup(pg_code)
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
            %     [S{i},MT{i}] = generate_sg(i); 
            %     % determine pg corresponding to each sgi %4i , ...\n',cell2mat(pg))
            %     pg{i} = identify_pointgroup( reshape(uniquecol_( reshape(S{i}(1:3,1:3,:),9,[]) ),3,3,[]) );
            %     % count number of symmetries
            %     nsyms{i} = size(MT{i},1);
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
              96   96    3    6    6    6    6   12   12 ]; 
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
              20   14    3    6    3    3    3    6    6 ]; 
          
            sg_code = find( (pg_database==pg_code) );
            
        end
        
        function brv_name     = decode_holohodry(pg_code)
            % point group dataset
            brav={'triclinic','','','','monoclinic','','','orthorhombic', ...
                  '','','','','trigonal','','','','','','','tetragonal',...
                  '','','','','','','hexagonal','','','','','cubic'};
            % print point group name
            brv_name = brav{pg_code};
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
            sg_name = sg_name_database(sg_code);
        end

        function order        = get_order(S)

            import am_lib.*
            import am_dft.*
            
            check_ = @(x) any(abs(x(:))>am_lib.eps);

            if size(S,1)==3
                order = 1; X=S(1:3,1:3);
                while check_(X - eye(3))
                    order = order + 1;
                    X = X*S;
                end
            elseif size(S,1) == 4
                order = 1; X=S(1:4,1:4);
                while check_(X - eye(4))
                    order = order + 1;
                    X = X*S; X(1:3,4)=mod_(X(1:3,4));
                end
            end
        end
        
        function [S,MT]       = generate_sg(sg_code)

            import am_lib.*
            import am_dft.*

            % initialize
            S = zeros(4,4,192);

            % load recipe
            recipe = get_recipe(sg_code);

            % mix ingredients 
            nsyms = 0; k = 0;
            nsyms = nsyms+1; S(1:4,1:4,nsyms) = eye(4);
            k = k+1; if recipe(1) == '1'
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

            % bake the cake (add elements until no new elements are generated)
            nsyms_last=0; MT = zeros(192,192);
            while nsyms ~= nsyms_last
                nsyms_last=nsyms;
                for i = 1:nsyms_last
                for j = 1:nsyms_last
                    if MT(i,j)==0
                        A = S(:,:,i)*S(:,:,j); A(1:3,4)=mod_(A(1:3,4));
                        A_id = member_(A(:),reshape(S(:,:,1:nsyms),16,[]));
                        if A_id == 0
                            nsyms = nsyms+1; S(:,:,nsyms) = A; MT(i,j) = nsyms;
                        else
                            MT(i,j) = A_id;
                        end
                    end
                end
                end
            end

            MT = MT(1:nsyms,1:nsyms);
            S  = S(:,:,1:nsyms);

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

        
        % unit cells

        function [uc,pc,ic]   = get_cells(fposcar)
            % wrapper routine
            % fname = 'infile.supercell' (poscar)
            
            import am_lib.*
            import am_dft.*

            % time
            fprintf(' ... solving for cells and symmetries'); tic
            
            % load poscar
            uc = load_poscar(fposcar);

            % get primitive cell
            [pc,p2u,u2p] = get_primitive_cell(uc);  % write_poscar(get_supercell(pc,eye(3)*2),'POSCAR.2x2')

            % get irreducible cell
            [ic,i2p,p2i] = get_irreducible_cell(pc); 

            % complete mapping
            u2i = round(p2i(u2p)); i2u = round(p2u(i2p));

            % save mapping to cells
            pc.p2i = p2i; pc.i2p = i2p;
            pc.p2u = p2u; pc.u2p = u2p;

            ic.i2p = i2p; ic.p2i = p2i;
            ic.i2u = i2u; ic.u2i = u2i;

            uc.u2p = u2p; uc.p2u = p2u;
            uc.u2i = u2i; uc.i2u = i2u;

            % save bas2pc and tau2pc to convert [uc-frac] to [pc-frac]
            uc.bas2pc = pc.bas/uc.bas; uc.tau2pc = pc.bas\uc.bas;
            
            % print basic symmetry info
            [~,H,~,R] = get_symmetries(pc); 
            hg_code = identify_pointgroup(H); fprintf(', %s',decode_holohodry(hg_code));
            pg_code = identify_pointgroup(R); fprintf(', %s',decode_pg(pg_code));
            
            % beta (POSSIBLE SPACE GROUPS)
            sg_code = identify_spacegroup(pg_code);
            fprintf(' (%s)',strrep(cell2mat(join(decode_sg(sg_code),',')),' ',''));
            
            fprintf(' (%.f secs)\n',toc);
        end

        function [uc]         = get_supercell(pc,B)
            
            import am_lib.*
            import am_dft.*

            % basic check
            if mod(det(B),1)~=0; error('determinant of B must be an integer'); end 

            % generate primitive lattice vectors
            n=round(sum(abs(B),1)); [Y{3:-1:1}]=ndgrid(1:n(1),1:n(2),1:n(3)); nLs=prod(n); L=reshape(cat(3+1,Y{:})-1,[],3).'; 

            % expand atoms, coordinates supercell fractional, and reduce to primitive supercell
            X=uniquecol_([ reshape(repmat([1:pc.natoms],nLs,1),1,[]); mod_(inv(B)*osum_(L,pc.tau,2)) ]);

            % create mapping
            u2p = X(1,:); [~,p2u]=unique(u2p); p2u=p2u(:).';

            % define irreducible cell creation function and make structure
            uc_ = @(uc,tau,B,s2u) struct('units','frac','bas',uc.bas*B,'bas2pc',inv(B),'tau2pc',B,...
                'symb',{{uc.symb{unique(uc.species(s2u))}}},'mass',uc.mass,'nspecies',sum(unique(uc.species(s2u)).'==uc.species(s2u),2).', ...
                'natoms',numel(s2u),'tau',tau,'species',uc.species(s2u));
            uc = uc_(pc,X(2:4,:),B,u2p);
            
            % add maps
            uc.u2p = u2p; 
            uc.p2u = p2u; 
            uc.u2i = pc.p2i(uc.u2p); 
            uc.i2u = uc.p2u(pc.i2p);
            
        end
        
        function [dc,idc]     = get_displaced_cell(pc,bvk,n,kpt,amp,mode,nsteps)
            % Note: can set mode=[] for interactive selection
            % n=[4;4;4]; kpt=[0;0;1/4]; amp=10; mode=6; nsteps=51;
            % [~,md] = get_displaced_cell(pc,bvk,n,kpt,amp,mode,nsteps); 
            % clf; [F]=plot_md_cell(md,'view',[0;1;0]); movie(F,3); % save_poscar(md,'POSCAR_test')
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

            import am_lib.*
            import am_dft.*
            
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
                [cc,c2d,d2c] = get_primitive_cell(dc); cc.bas = double(cc.bas);

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
        
        function [h]          = plot_cell(pc)
            
            import am_lib.*
            import am_dft.*
            
            % initialize figure
            set(gcf,'color','w'); hold on;
            
            % plot atoms
            h = scatter3_(pc.bas*pc.tau,50*sqrt(pc.mass(pc.species)),pc.species,'filled');
            
            % plot pc boundaries
            plothull_(pc.bas*[0,1,0,1,0,1,0,1;0,0,1,1,0,0,1,1;0,0,0,0,1,1,1,1]);

            hold off; daspect([1 1 1]); box on;
            
        end
        
        function [F]          = plot_md_cell(md,varargin)
            % n=[4;4;4]; kpt=[0;0;1/4]; amp=10; mode=6; nsteps=51;
            % [~,md] = get_displaced_cell(pc,bvk,n,kpt,amp,mode,nsteps); 
            % clf; [F]=plot_md_cell(md,'view',[0;1;0]); movie(F,3); % save_poscar(md,'POSCAR_test')
            
            import am_lib.*
            import am_dft.*
            
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
        
        
        % brillouin zones

        function [fbz,ibz]    = get_zones(pc,n)
            
            import am_lib.*
            import am_dft.*
            
            % get full brillouin zone
            [fbz] = get_fbz(pc,n);

            % get irreducible zone
            [ibz,i2f,f2i] = get_ibz(fbz,pc);

            % save mapping to zones
            fbz.f2i = f2i; fbz.i2f = i2f;
            ibz.i2f = i2f; ibz.f2i = f2i;
        end
        
        function [bzp]        = get_bz_path(pc,n,brav)

            import am_lib.*
            import am_dft.*
            
            % define kpoint path
            if     contains( lower(brav), 'fcc-short' )
                G=[0;0;0];  X1=[0;1;1]/2; X2=[2;1;1]/2; 
                L=[1;1;1]/2; K=[6;3;3]/8;
                % short path
                ql={'G','X','K','G','L'}; 
                qs=[G,X2,K,G]; 
                qe=[X1,K,G,L];
            elseif contains( lower(brav), 'fcc' )
                G=[0;0;0];  X1=[0;1;1]/2; W=[1;3;2]/4;
                U=[2;5;5]/8; L=[1;1;1]/2; K=[3;6;3]/8;
                % long path
                ql={'G','X','W','K','G','L','U','W','L','K'}; 
                qs=[G,X1,W,K,G,L,U,W,L]; 
                qe=[X1,W,K,G,L,U,W,L,K];
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
            bzp_ = @(recbas,ql,qt,nks,x,k) struct('units','frac-recp', ...
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
        
        function [bzs]        = get_bz_surf(pc,n,vy,vx)
            % bzs=get_bz_surf(pc,[101,101],[1;0;0],[0;1;0]);
            % bzs=get_bvk_dispersion(bvk,bzs);
            % plot_bz_surf(bzs,1)
            
            import am_lib.*
            import am_dft.*
            
            % get number of kpoints
            nks=prod(n); recbas = inv(pc.bas).';

            % get surface in [cart-recp] then convert back to [frac-recp]
            Q_ = @(n) [0:(n-1)]./(n-1); [Y{1:2}]=meshgrid(Q_(n(2)),Q_(n(1)));
            k = Y{1}(:).'.*vx + Y{2}(:).'.*vy; k = recbas\k;

            % create path object
            bzs_ = @(recbas,n,nks,k) struct('units','frac-recp', ...
                'recbas',recbas,'nks',nks,'n',n,'k',k);
            bzs = bzs_(recbas,n,nks,k);
        end

        function [bzl]        = get_bz_line(pc,n,vs,ve)
            % surf: vs and ve are the start and end vectors in cart
            
            import am_lib.*
            import am_dft.*
            
            % get number of kpoints
            nks=n; recbas = inv(pc.bas).';
            
            % define path (includes both boundaries)
            path_ = @(k,q,N) cumsum([zeros(3,1),repmat((k-q)/(N-1),1,N-1)],2)+repmat(q,1,N);
            x_    = @(k,q,N) [0, repmat(norm((k-q)/(N-1)),1,N-1) ];
            
            % build path in recp.-frac
            k(1:3,:) = path_(ve-vs,[0;0;0],n); k = recbas\k;
            x(    :) =    x_(ve-vs,[0;0;0],n); x = cumsum(x);

            % create path object
            bzl_ = @(recbas,n,nks,x,k) struct('units','frac-recp',...
                'recbas',recbas,'nks',nks,'n',n,'x',x,'k',k);
            bzl = bzl_(recbas,n,nks,x,k);
        end

        function [bz]         = get_dispersion(model,bz,flag)
            % Get phonon and electron dispersions on bz.
            % model is either bvk or tb
            %
            % Q: Are the eigenvalues and eigenvectors the same at
            %    symmetry-equivalent positions in reciprocal space?
            % A: Eigenvalues, yes; eigenvectors are (to within a sign)
            %    rotated by the same point symmetry relating the different
            %    k-points. Check it out with:
            %
            %     % generate point symmetries [recp-frac]
            %     [~,~,~,R] = get_symmetries(pc); nRs = size(R,3); R = permute(R,[2,1,3]);
            %     % get orbit of a rank k point
            %     k = matmul_(R,rand(3,1)); nks = size(k,2); 
            %     % get eigenvectors and eigenvalues at each orbit point
            %     for i = 1:nks
            %     input = num2cell([bvk.fc{:},(pc.bas.'\k(:,i)).',bvk.mass]); % [cart]
            %     [U(:,:,i),hw(:,i)] = eig( force_hermiticity_(bvk.D(input{:})) ,'vector');
            %     end
            % 
            %     % covert symmetries [recp-frac] -> [recp-cart]
            %     sym_rebase_ = @(B,R) matmul_(matmul_(B,R),inv(B));
            %     R_cart = sym_rebase_(inv(pc.bas).',R);
            % 
            %     % rotate the first eigenvector and compare with other kpoints
            %     X = matmul_( blkdiag_(R_cart,R_cart) ,  U(:,1,1) );
            %     abs(X) - abs(squeeze(U(:,1,:))) < am_lib.eps
            % 
            %     % compare eigenvalues
            %     hw-hw(:,1) < am_lib.eps
            %

            import am_lib.*
            import am_dft.* 
            
            fprintf(' ... computing dispersion '); tic; 
            switch lower(strtrim(flag))
                case 'electron'

                    % get eigenvalues
                    bz.nbands = model.nbands; 
                    bz.E = zeros(bz.nbands,bz.nks); 
                    bz.V = zeros(bz.nbands,bz.nbands,bz.nks);
                    for i = 1:bz.nks
                        % define input ...
                        input = num2cell([model.vsk{:},[bz.recbas*bz.k(:,i)].']);
                        % ... and evaluate (V are column vectors)
                        [bz.V(:,:,i),bz.E(:,i)] = eig(  force_hermiticity_(model.H(input{:})) ,'vector');
                        % sort energies
                        [bz.E(:,i),inds]=sort(bz.E(:,i)); bz.V(:,:,i)=bz.V(:,inds,i);
                    end
                    
                case 'phonon'
            
                    % get eigenvalues
                    bz.nbranches = model.nbranches; 
                    bz.hw = zeros(bz.nbranches,bz.nks); 
                    bz.U  = zeros(bz.nbranches,bz.nbranches,bz.nks);
                    for i = 1:bz.nks
                        % define input ...
                        % input = num2cell([bvk.fc{:},bz.k(:,i).',bvk.mass]); % [pc-frac]
                        input = num2cell([model.fc{:},(bz.recbas*bz.k(:,i)).',model.mass]); % [cart]
                        % ... and evaluate (U are column vectors)
                        [bz.U(:,:,i),bz.hw(:,i)] = eig( force_hermiticity_(model.D(input{:})) ,'vector');
                        % correct units
                        bz.hw(:,i) = sqrt(real(bz.hw(:,i))) * am_lib.units_eV;
                        % sort energies
                        [bz.hw(:,i),inds]=sort(bz.hw(:,i)); bz.U(:,:,i)=bz.U(:,inds,i);
                    end
            end
            fprintf('(%.f secs)\n',toc);
        end
        
        function [nesting]    = get_nesting(fbz,ibz,degauss,Ep)
            % get nesting function on ibz at probing energies Ep using smearing degauss
            % degauss, degauss = 0.04 61x61x61 kpoint mesh
            
            import am_lib.*
            import am_dft.*
            
            % number of probing energies
            nEps = numel(Ep);
            
            % get momentum conserving q-point triplets
            qqq = get_qqq(fbz,ibz);

            % copy irreducible energies onto fbz mesh
            E = ibz2fbz(fbz,ibz,ibz.E);
            
            % compute spectral function A on the full mesh
            A = reshape( sum(lorentz_(( reshape(E,[1,size(E)]) - Ep(:) )/degauss)/degauss,2) , [nEps,fbz.n]);

            % compute nesting on ibz and transfer to ibz
            nesting = zeros(nEps,ibz.nks); access_ = @(x,i) reshape(x(i),[],1); m = 2;
            for j = 1:nEps
                nesting(j,:) = accumarray( fbz.f2i( qqq(1,:,m)).' , ...
                           access_(A(j,:,:,:),qqq(2,:,m))     ...
                        .* access_(A(j,:,:,:),qqq(3,:,m))   , [], @sum  )./prod(ibz.n);
            end
        end
        
        function plot_interpolated(fbz,bzp,x,varargin)
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
            plot(bzp.x,fftinterp_(x,bzp.k,fbz.n),varargin{:});
            axs_(gca,bzp.qt,bzp.ql); axis tight; xlabel('Wavevector k');
        end
              
        function plot_dispersion(model,bzp,flag,varargin)
            % model is either bvk or tb
            
            import am_lib.*
            import am_dft.*
            
            % define figure properties
            fig_ = @(h)       set(h,'color','white');
            axs_ = @(h,qt,ql) set(h,'Box','on','XTick',qt,'Xticklabel',ql);
            fig_(gcf);
            
            switch lower(strtrim(flag))
                case 'electron'
                    % get electron band structure along path
                    bzp = get_dispersion(model,bzp,flag);
                    % plot results
                    plot(bzp.x,sort(bzp.E),varargin{:});
                    axs_(gca,bzp.qt,bzp.ql); axis tight; ylabel('Energy [eV]'); xlabel('Wavevector k');
                case 'phonon'
                    % get phonon band structure along path
                    bzp = get_dispersion(model,bzp,flag);
                    % plot results
                    plot(bzp.x,sort(real(bzp.hw)*1E3),'-k',bzp.x,-sort(abs(imag(bzp.hw))),varargin{:});
                    axs_(gca,bzp.qt,bzp.ql); axis tight; ylabel('Energy [meV]'); xlabel('Wavevector k');
            end
            

        end
        
        function plot_nesting(ibz,fbz,bzp,degauss,Ep,varargin)
            
            import am_lib.*
            import am_dft.*

            % plot results
            plot_interpolated(fbz,bzp, ibz2fbz(fbz,ibz,get_nesting(fbz,ibz,degauss,Ep)) ,varargin{:})
        end
        
        function plot_bz(fbz)
            
            import am_lib.*
            import am_dft.*
            
            % initialize figure
            set(gcf,'color','w'); hold on;
            
            % plot points
            % h = scatter3_(uc2ws(fbz.recbas*fbz.k,fbz.recbas),'.');
            
            % generate all possible points halfway between reciprocal lattice vectors
            P = [-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1;...
                 -1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1;...
                 -1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1]/2;
            P = uc2ws(fbz.recbas*P,fbz.recbas); P = uniquecol_(P);

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
            X = uc2ws(X*0.999,fbz.recbas)/0.999; X = uniquecol_(X);

            % plot convex hull 
            plothull_(X);

            hold off; daspect([1 1 1]); box on;
        end
        
        function plot_bz_path(bzp)
            
            import am_lib.*
            import am_dft.*

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
        
        function plot_bz_surf(bzs,band)
            % bzs=get_bz_surf(pc,[101,101],[1;0;0],[0;1;0]);
            % bzs=get_bvk_dispersion(bvk,bzs);
            % plot_bz_surf(bzs,1)
            
            import am_lib.*
            import am_dft.*

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
        

        % phonons (harmonic)

        function [bvk,pp]     = get_bvk(pc,uc,md,cutoff2,flags)
            
            import am_lib.*
            import am_dft.*

            if nargin<5; flags='-identify -model -fit -enforce'; end

            % get irreducible shells
            fprintf(' ... identifying pairs ');
            if contains(flags,'-identify'); tic;
                [bvk,pp] = get_pairs(pc,uc,cutoff2);
                fprintf('(%.f secs)\n',toc);
            else; fprintf('(skipped)\n'); end
            
            % force constant model
            fprintf(' ... determining harmonic force constant interdependancies and dynamical matrix ');
            if contains(flags,'-model'); tic;
                bvk = get_bvk_model(bvk,pp,uc);
                fprintf('(%.f secs)\n',toc);
            else; fprintf('(skipped)\n'); end
            
            % get force constants
            fprintf(' ... solving for harmonic force constants ');
            if contains(flags,'-fit'); tic;
                bvk = get_bvk_force_constants(uc,md,bvk,pp);
                fprintf('(%.f secs)\n',toc);
            else; fprintf('(skipped)\n'); end

            % enforce asr
            fprintf(' ... enforcing acoustic sum rules '); 
            if contains(flags,'-enforce'); tic;
                bvk = set_bvk_acoustic_sum_rules(bvk,pp);
                fprintf('(%.f secs)\n',toc);
            else; fprintf('(skipped)\n'); end
        end
 
        function [bvk]        = get_bvk_model(ip,pp,uc)
            
            import am_lib.*
            import am_dft.*

            % set sym digits
            digits(10);

            % get form of force constants for irreducible prototypical bonds
            for i = 1:ip.nshells
                % use stabilzer group to determine crystallographic symmetry relations; A*B*C' equals kron(C,A)*B(:)
                W = sum(kron_( pp.Q{1}(1:3,1:3,ip.s_ck(:,i)) , pp.Q{1}(1:3,1:3,ip.s_ck(:,i)))-eye(9),3);

                % enforce intrinsic symmetry (immaterial order of differentiation: c == c.')
                F = zeros(9,9); F(sub2ind([9,9],[1:9],[1,4,7,2,5,8,3,6,9])) = 1; W = W + F-eye(9);

                % get linearly-independent nullspace and normalize to first nonzero element
                W=real(null(W)); W=frref_(W.').'; W(abs(W)<am_lib.eps)=0; W(abs(W-1)<am_lib.eps)=1; W=W./accessc_(W,findrow_(W.').');

                % define parameters
                c = sym(sprintf('c%02i_%%d%%d',i),[3,3],'real'); c = c(findrow_(double(W).'));

                % get symmetry adapted force constants
                phi = reshape( sym(W,'d')*c(:), [3,3]);

                % save important stuff (sort W to be in line with c, matlabFunction sorts D variables)
                [sav.c{i},n] = sort(c(:).'); sav.W{i} = W(:,n); sav.phi{i} = phi;
            end

            % create bvk structure
            bvk_ = @(pp,ip,sav) struct('units','cart','bas',pp.bas2pc*pp.bas, ...
                'symb',{pp.symb},'mass',pp.mass,'species',pp.species(pp.p2u),'cutoff',pp.cutoff,'natoms',pp.pc_natoms,...
                'nbranches',3*pp.pc_natoms,'nshells',size(sav.W,2),'W',{sav.W},'phi',{sav.phi},'d',ip.d,'v',ip.v,'xy',ip.xy);
            bvk = bvk_(pp,ip,sav);

            % [cart] define function to get bond vector
            vec_ = @(xy) uc2ws(uc.bas*(uc.tau(:,xy(2,:))-uc.tau(:,xy(1,:))),pp.bas);
            
            % construct symbolic dynamical matrix
            D=sym(zeros(bvk.nbranches)); kvec=sym('k%d',[3,1],'real'); mass=sym('m%d',[1,numel(pp.i2u)],'positive');
            for p = 1:pp.pc_natoms
            for j = 1:pp.npairs(p)
                % get indicies and already permute xy,mn,mp,np if necessary by iq
                i = pp.i{p}(j);iq = pp.iq{p}(j);
                x = pp.c{p}(1); y = pp.o{p}(j,1); xy = [x;y];
                m = pp.u2p(x); mp = [1:3]+3*(m-1); n = pp.u2p(y); np = [1:3]+3*(n-1);

                % rotate force constants and bond vector (NOTE #1)
                rij = vec_(xy); rij(abs(rij)<am_lib.eps) = 0;
                phi = sym(pp.Q{1}(1:3,1:3,iq)) * permute(bvk.phi{i},pp.Q{2}(:,iq)) * sym(pp.Q{1}(1:3,1:3,iq)).';
                
                % build dynamical matrix
                D(mp,np) = D(mp,np) + phi .* exp(sym(2i*pi * rij(:).','d') * kvec(:) );
            end
            end

            % simplify (speeds evaluation up significantly later)
            D = simplify(D,'steps',500);
            
            % multiply by 1/sqrt(mass)
            mass = repelem(mass(pp.species(pp.p2u)),1,3); mass = 1./sqrt(mass.' * mass); D = D .* mass;
            
            % attach symbolic dynamical matrix to bvk
            bvk.D = matlabFunction(D);
        end

        function [bvk]        = get_bvk_force_constants(uc,md,bvk,pp)
            % Extracts symmetry adapted force constants.
            %
            % Q: What is the difference of these methods? 
            % A: Force constants extracted, dispersions look different:
            %    Method 1:    8.9932    0.2321   -0.4806   -0.1065   -0.0337   -4.9099   14.6953   -0.5217    0.6199   -0.2590
            %    Method 2:   10.0033    0.4073   -0.5706    0.0775   -0.0314   -5.3528   16.0851   -0.5830    0.2604   -0.4623
            % Q: How do the methods comapre in terms of memory, speed, and
            %    ability to reproduce AIMD forces? 
            % A:        METHOD     1           2      
            %              R^2     0.668       0.633  
            %      MATRIX SIZE     3*m*nFCs    3*n*m  
            %          TIME [s]    0.322       0.031  
            %      TIMES SLOWER    10x         1x     

            import am_lib.*
            import am_dft.*

            % [cart] get displacements and forces
            u = matmul_( md.bas, mod_( md.tau-uc.tau +.5 )-.5 );
            f = matmul_( md.bas, md.force );

            % select a method (method 1 is the most accurate)
            algo = 1;
            switch algo
                case 1
                    % Method incorporating full intrinsic and crystalographic symmetries 
                    % Note that ASR are automatically enforced when this method is used!
                    % F [3m * 1] = - U [3m * nFcs] * FC [ nFCs * 1]: n pairs, m atoms ==> FC = - U \ F
                    % 
                    % Q: How are roto-inversions incorporated?
                    % A: Check with:
                    %
                    %         % define input prep
                    %         prep_ = @(x) [x(1),x(2),x(3),0,0,0,0,0,0;0,0,0,x(1),x(2),x(3),0,0,0;0,0,0,0,0,0,x(1),x(2),x(3)];
                    %         % select irreducible pair j and symmetry i
                    %         j=2;i=9;
                    %         % define symbolic displacement u and force constants c
                    %         u=sym('u_%d',[3,1]); c=sym('c_%d',[size(bvk.W{j},2),1]); 
                    %         % get rotations and inverse
                    %         R=pp.Q{1}(1:3,1:3,pp.q{1}(i)); iR=pp.Q{1}(1:3,1:3,pp.iq{1}(i));
                    %         % check without rotation
                    %         prep_(u)*bvk.W{j} - equationsToMatrix(reshape(bvk.W{j}*c,3,3)*u,c)
                    %         % check with rotation
                    %         R*prep_(iR*u)*bvk.W{j} - equationsToMatrix(R*reshape(bvk.W{j}*c,3,3)*iR*u,c)
                    %
                    
                    % get U matrix and indexing I
                    [U,I] = get_bvk_U_matrix(bvk,pp,u);
                    
                    % solve for force constants
                    fc = - U \ f(I);
                    
                    % save force constants
                    s_id = repelem([1:bvk.nshells],cellfun(@(x)size(x,2),bvk.W));
                    for s = 1:bvk.nshells; bvk.fc{s} = double(fc(s_id==s).'); end
                case 2
                    % basic method using full 3x3 second-order tensor (ignores intrinsic force constant symmetry)
                    % F [3 * m] = - FC [3 * 3n] * U [3n * m]: n pairs, m atoms ==> FC = - F / U
                    phi=[];
                    for m = 1:pp.pc_natoms
                        % get forces : f = [ (x,y,z), (1:natoms)*nsteps ] 
                        % get displacements : u [ (x,y,z)*orbits, (1:natoms)*nsteps ]
                        %  ... and solve for the generalized force constants: FC = - f / u 
                        fc = - reshape( f(:,pp.c{m},:) ,3,pp.ncenters(m)*md.nsteps) /...
                               reshape( u(:,pp.o{m},:) ,3*pp.npairs(m),pp.ncenters(m)*md.nsteps);
                        % reshape into 3x3 matrices
                        phi = double(cat(3,phi,reshape(fc,3,3,[])));
                    end
                    
                    % transform fc from orbit to irrep
                    q = cat(1,pp.q{:}); phi = matmul_( matmul_( pp.Q{1}(1:3,1:3,q) ,phi), permute(pp.Q{1}(1:3,1:3,q),[2,1,3]) );
                    for j = 1:size(phi,3); phi(:,:,j) = permute( phi(:,:,j), pp.Q{2}(:,q(j)) ); end

                    % solve for symmetry adapted force constants : A x = B
                    for i = 1:bvk.nshells
                        ex_ = cat(1,pp.i{:})==i;
                        if any(ex_)
                            A = repmat(double(bvk.W{i}),sum(ex_),1);
                            B = reshape(phi(:,:,ex_),[],1);
                            % get force constants as row vectors
                            bvk.fc{i} = reshape( A \ B , 1, []);
                        end
                    end
            end
        end

        function [md]         = run_bvk_md(bvk,pp,uc,dt,nsteps,Q,T)
            % set time step [ps ~ 0.1], number of MDs steps, Nose-Hoover "mass" Q, and temperature T [K]
            % dt = 0.1; nsteps = 10000; Q = 1; T = 300; [md] = run_bvk_md(bvk,pp,uc,dt,nsteps,Q,T)

            import am_lib.*
            import am_dft.*

            % build force constants
            for m = 1:pp.pc_natoms
                phi{m} = zeros(3,3*pp.npairs(m));
            for j = 1:pp.npairs(m)
                % get indicies
                i = pp.i{m}(j); iq = pp.iq{m}(j);
                % get irrep force constant indicies
                iphi = reshape(bvk.W{i}*bvk.fc{i}(:),3,3);
                % rotate force constants from irrep to orbit
                phi{m}(1:3,[1:3]+3*(j-1)) = pp.Q{1}(1:3,1:3,iq) * permute(iphi,pp.Q{2}(:,iq)) * pp.Q{1}(1:3,1:3,iq).';
            end
            end

            % allocate arrays and constants [cart]
            k_boltz = 8.6173303E-5; % [eV/K]
            u = single(zeros(3,uc.natoms,nsteps)); KE=zeros(1,nsteps);
            v = single(zeros(3,uc.natoms,nsteps)); PE=zeros(1,nsteps); 
            f = single(zeros(3,uc.natoms,nsteps));

            % set initial value conditions: small displacement just to get atoms moving
            u(:,:,1) = (.5-rand(3,uc.natoms))*0.00001; 
            v(:,:,1) =    zeros(3,uc.natoms); 

            % % set velocity gaussianly distributed sqrt( k_B * K / amu) = 0.000911856 [Ang/fs]
            % r1 = rand(3,uc.natoms); r2 = rand(3,uc.natoms); u(:,:,1) = zeros(3,uc.natoms);
            % v(:,:,1) = sqrt(-2*log(r1)).*cos(2*pi*r2).*sqrt(T./uc.mass(uc.species)/2) * 0.000911856;
            
            % run md using verlet algorithm
            fprintf('%10s   %10s     %10s   %10s   %10s \n','step','temp','PE','KE','PE+KE');
            for j = 1:nsteps
                % compute force [eV/Ang]
                for m = 1:pp.pc_natoms; f(:,pp.c{m},j) = - phi{m} * reshape(u(:,pp.o{m},j), size(pp.o{m}).*[3,1]); end            

                % compute potential energy [eV/Ang * Ang -> eV]
                PE(j) = - flatten_(u(:,:,j)).'*flatten_(f(:,:,j));

                % 4) compute kinetic energy : amu * (Ang/fs)^2 -> 103.6382 eV
                KE(j) = uc.mass(uc.species)*sum((uc.bas*v(:,:,j)).^2,1).'/2 * 103.6382;

                % 6) current temperature
                Tj = 2/3*KE(j)/uc.natoms/k_boltz;
                
                % 5) compute Nose-Hoover drag: p_eta = KE - TE
                nosehoover = v(:,:,j)/Q * ( Tj - T ) / uc.natoms;

                % 6) get acceleration
                acc = f(:,:,j) ./ uc.mass(uc.species);

                % ***) update md [frac]: x' = x + v * dt; v' = v + a * dt; Nose-Hoover dv/dt becomes a - p_eta / Q * v;
                if j ~= nsteps
                    u(:,:,j+1) = u(:,:,j) + dt * v(:,:,j);
                    v(:,:,j+1) = v(:,:,j) + dt * (acc - nosehoover);
                end

                % print
                if mod(j,50)==1
                    fprintf('%10i   %10f K   %10f   %10f   %10f \n',j,Tj,PE(j),KE(j),PE(j)+KE(j));
                    %
                    figure(1); set(gcf,'color','white'); 
                    plot([1:nsteps],[KE;PE;KE+PE].');legend('KE','PE','KE+PE')
                    xlabel('step'); ylabel('energy [eV]'); xlim([0,nsteps]); drawnow;
                end
            end

            % convert [cart] to [frac] and u to tau
            tau = matmul_(inv(uc.bas),u)+uc.tau;
            v   = matmul_(inv(uc.bas),v);
            f   = matmul_(inv(uc.bas),f);
            
            % define md creation function [frac]
            md_ = @(uc,force,tau,vel,dt) struct('units','frac',...
                'bas',uc.bas,'symb',{{uc.symb{:}}},'mass',uc.mass,'nspecies',uc.nspecies, ...
                'natoms',uc.natoms,'force',force,'tau',tau,'vel',vel,'species',uc.species, ...
                'dt',dt,'nsteps',size(tau,3));
            md = md_(uc,f,tau,v,dt);
        end

        function [bvk]        = interpolate_bvk(bvk_1,bvk_2,n)
            % interpolates force constants and masses from bvk_1 and bvk_2 on n points (includes end points)

            import am_lib.*
            import am_dft.*
            
            bvk_ = @(bvk,mass,fc) struct('units','cart','bas',bvk.bas,'recbas',bvk.recbas,'natoms',bvk.natoms,'mass',mass, ...
                'nshells',bvk.nshells,'W',{bvk.W},'shell',{bvk.shell},'nbranches',bvk.nbranches,'D',bvk.D,'fc',{fc});

            fc_interp = nlinspace( [bvk_1.fc{:}] , [bvk_2.fc{:}] , n );
            mu_interp = nlinspace( [bvk_1.mass]  , [bvk_2.mass]  , n );

            % get dimensions
            for i = 1:numel(bvk_1.fc); m(i) = numel(bvk_1.fc{i}); end; E = cumsum(m); S = E - m + 1;

            % create bvks cells
            for i = 1:n
                for j = 1:numel(E); fc{j} = fc_interp(S(j):E(j),i).'; end
                bvk(i) = bvk_(bvk_1,mu_interp(:,i).',fc);
            end
        end

        function plot_bvk_vs_aimd(uc,md,bvk,pp,bvt,pt)
            
            import am_lib.*
            import am_dft.*
            
            % select algo 2 which is faster
            algo = 2;
            
            % [cart] get displacements
            u = matmul_( md.bas, mod_( md.tau-uc.tau +.5 )-.5 );
            % [cart] get aimd forces
            f_aimd= matmul_( md.bas, md.force );
            % [cart]  get harmonic forces
            f_har = get_bvk_forces(bvk,pp,u,algo); 
            % [cart] get anharmnoic forces
            f_anh = 0; 
            if nargin == 6; if and(~isempty(bvt),~isempty(pt))
            f_anh = get_bvt_forces(bvt,pt,u,algo);
            end; end
            
            % plot correlation
            plotcorr_(f_anh(:)+f_har(:),f_aimd(:)); xlabel('aimd force [eV/Ang]'); ylabel('bvk force [eV/Ang]');
        end

        function [T,KE,PE,q_sk] = plot_md_stats(uc,md,fbz,bvk)
            
            import am_lib.*
            import am_dft.*
            
            % match uc to md
            uc = match_cell(uc,md);
            
            % [cart] get displacements and forces; [amu * (Ang/fs)^2 ] = 103.6382 [eV]
            u = matmul_( md.bas, mod_( md.tau-uc.tau +.5 )-.5 );
            f = matmul_( md.bas, md.force );
            v = matmul_( md.bas, md.vel ) * sqrt(103.6382);
            
            % get potentoal energy : [eV/Ang * Ang] = [eV]
            PE(:,1) = -dot(reshape(u,[],md.nsteps),reshape(f,[],md.nsteps),1)/2; 
            % get kinetic energy : [amu * (Ang/fs)^2 ] = 103.6382 [eV] (incorporate above in to v)
            KE(:,1) =  reshape(sum(uc.mass(uc.species).*dot(v,v,1),2),1,[])/2;
            
            % get temperature : amu * (Ang/fs)^2/ k_B = 1.20267E6 K
            k_boltz = 8.6173303E-5; T = 2/3*KE/uc.natoms/k_boltz;   

            % perform normal-mode analysis
            if nargin == 4
                % convert phonon energies back to wierd units
                hw = real(fbz.hw)./am_lib.units_eV;
                
                % get normal transformations from uc eigenvector (Wallace p 115 eq 10.48)
                U   = expand_bvk_eigenvectors(bvk,uc,fbz);
                u2q = U'  .* repelem(sqrt(uc.mass(uc.species)).',3,1).';
                v2p = U.' .* repelem(sqrt(uc.mass(uc.species)).',3,1).';

                % get normal modes (Wallace p 115 eq 10.49)
                q_sk = (u2q*reshape(u,[],md.nsteps));
                p_sk = (v2p*reshape(v,[],md.nsteps));

                % get energies (Wallace p 115 eq 10.53)
                PE(:,2) = dot(abs(q_sk).*hw(:),abs(q_sk).*hw(:),1)/2;
                KE(:,2) = dot(abs(p_sk)       ,abs(p_sk)       ,1)/2;
            end
                
            % if no output is requested, plot and print stuff
            if nargout == 0
                % print a few time steps
                Z = [[1:md.nsteps].',T(:,1),PE(:,1),KE(:,1),KE(:,1)+PE(:,1)];
                fprintf('%10s   %10s   %10s   %10s   %10s \n','step [#]','T [K]','PE [eV]','KE [eV]','PE+KE [eV]');
                if md.nsteps<100; fprintf('%10i   %10f   %10f   %10f   %10f \n',Z(1:end,:).'); else
                                  fprintf('%10i   %10f   %10f   %10f   %10f \n',Z(1:99,:).');
                                  fprintf('%10i   %10f   %10f   %10f   %10f \n',Z(100:50:end,:).');
                end    
                              
                % plot energies
                set(gcf,'color','w');
                if nargin == 4
                subplot(3,4,1:3); plot(1:md.nsteps,KE(:,1),'-',1:md.nsteps,PE(:,1),'-',1:md.nsteps,KE(:,1)+PE(:,1),'-',...
                                       1:md.nsteps,KE(:,2),'.',1:md.nsteps,PE(:,2),'.',1:md.nsteps,KE(:,2)+PE(:,2),'.');
                                legend('KE','PE','KE+PE','nKE','nPE','nKE+nPE');
                                xlabel('time step'); ylabel('energy [eV]'); axis tight;
                subplot(3,4,4); plot3(repelem([1:md.nsteps].',1,fbz.nks*bvk.nbranches),real(q_sk.*hw(:)).',real(p_sk).'); 
                                view([1 0 0]); box on; ax=axis; maxax=max(abs(ax(3:6))); axis([ax(1:2),[-1 1 -1 1].*maxax]); daspect([md.nsteps/maxax 1 1]);
                else
                subplot(3,1,1); plot(1:md.nsteps,KE(:,1),'-',1:md.nsteps,KE(:,2),'o',1:md.nsteps,KE(:,1)+PE(:,1),'^');
                                legend('KE','PE','KE+PE');
                                xlabel('time step'); ylabel('energy [eV]'); axis tight;
                end

                % plot position and velocity histograms
                nbins  = 101;
                dist_v = reshape(normc_(v),uc.natoms,md.nsteps); 
                dist_u = reshape(normc_(u),uc.natoms,md.nsteps); 
                bin_v  = linspace(0,max(dist_v(:)),nbins);
                bin_u  = linspace(0,max(dist_u(:)),nbins);
                hist_v = zeros(nbins-1,md.nsteps);
                hist_u = zeros(nbins-1,md.nsteps);
                for i = 1:md.nsteps
                    [hist_v(:,i)] = histcounts(dist_v(:,i),bin_v);
                    [hist_u(:,i)] = histcounts(dist_u(:,i),bin_u);
                end
                [Y{1:2}]=meshgrid(1:md.nsteps,bin_v(1:(nbins-1))); hist_v(:,1)=[]; Y{1}(:,1)=[]; Y{2}(:,1)=[];
                subplot(3,1,2); surf(Y{1},Y{2},hist_v,'Facecolor','interp','edgecolor','none');
                                view(2); axis tight; xlabel('time step'); ylabel('velocity [Ang/fs]');
                [Y{1:2}]=meshgrid(1:md.nsteps,bin_u(1:(nbins-1))); hist_u(:,1)=[]; Y{1}(:,1)=[]; Y{2}(:,1)=[];
                subplot(3,1,3); surf(Y{1},Y{2},hist_u,'Facecolor','interp','edgecolor','none');
                                view(2); axis tight; xlabel('time step'); ylabel('displacement [Ang]');
                
            end
        end

        
        % phonons (anharmonic)
        
        function [bvt,pt,bvk] = get_bvt(pc,uc,md,cutoff3,bvk,pp,flags)
            
            import am_lib.*
            import am_dft.*
            
            if nargin<7; flags='-identify -model -fit'; end

            % get irreducible shells
            fprintf(' ... identifying triplets '); 
            if contains(flags,'-identify'); tic;
                [bvt,pt] = get_triplets(pc,uc,cutoff3);
                fprintf('(%.f secs)\n',toc);
            else; fprintf('(skipped)\n'); end
            
            % get force constant model
            fprintf(' ... determining anharmonic force constant interdependancies ');
            if contains(flags,'-model'); tic;
                bvt = get_bvt_model(bvt,pt);
                fprintf('(%.f secs)\n',toc);
            else; fprintf('(skipped)\n'); end
            
            % get force constants
            fprintf(' ... solving for anharmonic (3rd-order) force constants '); 
            if contains(flags,'-fit'); tic;
                [bvt,bvk] = get_bvt_force_constants(uc,md,bvk,pp,bvt,pt);
                fprintf('(%.f secs)\n',toc);
            else; fprintf('(skipped)\n'); end
        end

        function [bvt]        = get_bvt_model(it,pt)
            % 
            % Q: What is the intrinsict symmetry of third-order force
            %    constants due to permutation of coordinate axes x,y,z?
            % A: Each page of the tensor is symmetry.
            %
            %     phi(:,:,1) =				   phi(1,:,:) =                
            %         [ c_111, c_211, c_311]   		[ c_111, c_211, c_311] 
            %         [ c_211, c_221, c_321]   		[ c_211, c_221, c_321] 
            %         [ c_311, c_321, c_331]   		[ c_311, c_321, c_331] 
            %     phi(:,:,2) =				   phi(2,:,:) =                
            %         [ c_211, c_221, c_321]   		[ c_211, c_221, c_321] 
            %         [ c_221, c_222, c_322]   		[ c_221, c_222, c_322] 
            %         [ c_321, c_322, c_332]   		[ c_321, c_322, c_332] 
            %     phi(:,:,3) =				   phi(3,:,:) =                
            %         [ c_311, c_321, c_331]   		[ c_311, c_321, c_331] 
            %         [ c_321, c_322, c_332]   		[ c_321, c_322, c_332] 
            %         [ c_331, c_332, c_333]   		[ c_331, c_332, c_333] 
            %
            %    And the same thing for the middle phi(:,1:3,:). Check if out with the code:
            %
            %         % get [x,y,z] flip/permutation operators
            %         a = reshape([1:27],3,3,3); p = perms(1:3).'; nps = size(p,2); F = zeros(27,27,nps); 
            %         for i = 1:nps; F(:,:,i) = sparse( flatten_(a), flatten_(permute(a,p(:,i))), ones(27,1), 27, 27); end
            %         % enforce intrinsic symmetry (immaterial order of differentiation: c == c.')
            %         W = sum(F-eye(27),3);
            %         % get linearly-independent nullspace and normalize to first nonzero element
            %         W=real(null(W)); W=frref_(W.').'; W(abs(W)<am_lib.eps)=0; W(abs(W-1)<am_lib.eps)=1; W=W./accessc_(W,findrow_(W.').');
            %         % define parameters
            %         c = sym(sprintf('c_%%d%%d%%d',i),[3,3,3],'real'); c = c(findrow_(double(W).'));
            %         % get symmetry adapted force constants
            %         phi = reshape( sym(W)*c(:), [3,3,3]);
            %
            % Q: How do second and third rank tensors rotate?
            % A:
            %         %% rotate 2nd rank tensor A
            %         clear;clc;rng(1); R=rand(3,3); A = rand(3,3);
            %         AR = zeros(3,3);
            %         for m = 1:3; for n = 1:3; 
            %         for i = 1:3; for j = 1:3
            %             AR(m,n) = AR(m,n) + A(i,j)*R(m,i)*R(n,j);
            %         end; end
            %         end; end
            %         AR - R*A*R.'
            %         AR - reshape(kron(R,R)*A(:),3,3)
            % 
            %         %% rotate 3rd rank tensor A
            %         clear;clc;rng(1); R=rand(3,3); A = rand(3,3,3);
            %         AR = zeros(3,3,3);
            %         for m = 1:3; for n = 1:3; for o = 1:3
            %         for i = 1:3; for j = 1:3; for k = 1:3
            %             AR(m,n,o) = AR(m,n,o) + A(i,j,k)*R(m,i)*R(n,j)*R(o,k);
            %         end; end; end
            %         end; end; end
            %         AR - reshape(kron(kron(R,R),R)*A(:),[3,3,3])
            %
            % Q: How do triple dot products work? For example, when
            %    fourier-transforming third-order force constants.
            % A: The triple dot product PHI ... u1 u2 u3, as described by Ziman,
            %    produces as a scalar and is given in Einstein summation as 
            %    PHI(ijk) u1(i) u2(j) u3(k). 
            %
            %         % apply double dot product
            %         F = zeros(3,1);
            %         for m = 1:3
            %         for i = 1:3; for j = 1:3
            %             F(m) = F(m) + AR(i,j,m)*u(i)*w(j);
            %         end; end
            %         end
            %
            % Q: What about a double dot product? For example, when
            %    third-order force constants multiply displacements to
            %    obtain forces.  
            % A: Check it out with the code below:
            %
            %         % double dot product of rotated 3rd rank tensor
            %         clear;clc;rng(1); R=rand(3,3); u = rand(3,1); w = rand(3,1); A = rand(3,3,3);
            %         % rotate A
            %         AR = zeros(3,3,3);
            %         for m = 1:3; for n = 1:3; for o = 1:3
            %         for i = 1:3; for j = 1:3; for k = 1:3
            %             AR(m,n,o) = AR(m,n,o) + A(i,j,k)*R(m,i)*R(n,j)*R(o,k);
            %         end; end; end
            %         end; end; end
            %         % apply double dot product
            %         F = zeros(3,1);
            %         for m = 1:3
            %         for i = 1:3; for j = 1:3
            %             F(m) = F(m) + AR(i,j,m)*u(i)*w(j);
            %         end; end
            %         end
            %         % equivalent formulations incorproating rotation:
            %         F - reshape(kron(kron(R,R),R)*A(:),9,3).'*reshape(u(:)*w(:).',[],1)
            %         % equivalent formulations without incorporating rotation
            %         F - reshape(          AR          ,9,3).'*reshape(u(:)*w(:).',[],1)
            %         F - [ u.'*AR(:,:,1)*w;  u.'*AR(:,:,2)*w ; u.'*AR(:,:,3)*w ]
            %         F - matmul_(AR,w).'*u
            %         F - reshape(AR,9,3).'*flatten_(u*w.')
            %         F - squeeze(sum(sum( AR.*(u*w.') ,1),2))
            %

            import am_lib.*
            import am_dft.*
        
            % set sym digits
            digits(10);

            % get flip operators
            a = reshape([1:27],3,3,3); p = perms(1:3).'; nps = size(p,2); F = zeros(27,27,nps); 
            for i = 1:nps; F(:,:,i) = sparse( flatten_(a), flatten_(permute(a,p(:,i))), ones(27,1), 27, 27); end

            % get form of force constants for irreducible prototypical bonds
            for i = 1:it.nshells
                % use stabilzer group to determine crystallographic symmetry relations; A*B*C' equals kron(C,A)*B(:)
                W = sum( kronpow_( pt.Q{1}(1:3,1:3,it.s_ck(:,i)) , 3) - eye(27), 3);

                % enforce intrinsic symmetry (immaterial order of differentiation: c == c.')
                W = W + sum(F-eye(27),3);

                % get linearly-independent nullspace and normalize to first nonzero element
                W=real(null(W)); W=frref_(W.').'; W(abs(W)<am_lib.eps)=0; W(abs(W-1)<am_lib.eps)=1; W=W./accessc_(W,findrow_(W.').');

                % define parameters
                c = sym(sprintf('c%02i_%%d%%d',i),[3,3,3],'real'); c = c(findrow_(double(W).'));

                % get symmetry adapted force constants
                phi = reshape( sym(W,'d')*c(:), [3,3,3]);

                % save important stuff (sort W to be in line with c, matlabFunction sorts D variables)
                [sav.c{i},n] = sort(c(:).'); sav.W{i} = W(:,n); sav.phi{i} = phi;
            end

            % augment bvk structure with triplet info
            bvt_ = @(pt,it,sav) struct('units','cart','bas',pt.bas2pc*pt.bas, ...
                'symb',{pt.symb},'mass',pt.mass,'species',pt.species(pt.p2u),'cutoff',pt.cutoff,'natoms',pt.pc_natoms,...
                'nshells',size(sav.W,2),'W',{sav.W},'phi',{sav.phi},'xyz',it.xyz);
            bvt = bvt_(pt,it,sav);

        end
        
        function [bvt,bvk]    = get_bvt_force_constants(uc,md,bvk,pp,bvt,pt)
            % Extracts symmetry adapted third-order force constant and
            % second order force constants 
            %
            % Q: How to ensure that the extraction process is working
            %    correctly? Is it working correctly?
            % A: Yes, it works. Checked it by doing the reverse operation and
            %    extracting force constants from a known displacement-force
            %    set using: 
            % 
            %         u = matmul_( md.bas, mod_( md.tau-uc.tau +.5 )-.5 );
            %         md.force = matmul_( inv(md.bas), get_bvt_forces(bvt,pt,u,1) + get_bvk_forces(bvk,pp,u,1) );
            %         [cmp] = get_bvt_force_constants(uc,md,bvk,pp,bvt,pt);
            %         plot([bvt.fc{:}],[bvt2.fc{:}],'.')
            %
            % Q: How to factor out displacements from third-rank tensors?
            % A: Need to create the matrix 
            % 
            %       [ u_1*w_1, u_2*w_1, u_3*w_1, u_1*w_2, u_2*w_2, u_3*w_2, u_1*w_3, u_2*w_3, u_3*w_3,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0]
            %       [       0,       0,       0,       0,       0,       0,       0,       0,       0, u_1*w_1, u_2*w_1, u_3*w_1, u_1*w_2, u_2*w_2, u_3*w_2, u_1*w_3, u_2*w_3, u_3*w_3,       0,       0,       0,       0,       0,       0,       0,       0,       0]
            %       [       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0, u_1*w_1, u_2*w_1, u_3*w_1, u_1*w_2, u_2*w_2, u_3*w_2, u_1*w_3, u_2*w_3, u_3*w_3]
            %
            %     using the Z matrix below:
            %
            %       Z = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            %            0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            %            0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            %            0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            %            0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
            %            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;
            %            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;
            %            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0;
            %            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1].';
            %
            %    That is: u = sym('u_%d',[3,1]); w = sym('w_%d',[3,1]); reshape(Z* reshape(u*w.',[],1),3,27)
            %

            import am_lib.*
            import am_dft.*

            % [cart] get displacements and forces
            u = matmul_( md.bas, mod_( md.tau-uc.tau +.5 )-.5 );
            f = matmul_( md.bas, md.force );

            % choose the algorthim to extract third order anharmonic force
            % constants; algorithm 2 is more accurate
            switch 2
                case 1
                    % subtract harmonic forces first than get anharmonic
                    % force constants from the residuals
                    
                    % [cart] subtract harmonic forces
                    if and(~isempty(bvk),~isempty(pp))
                        f = f - get_bvk_forces(bvk,pp,u);
                    end

                    % get U matrix and indexing I
                    [U,I] = get_bvt_U_matrix(bvt,pt,u);

                    % [cart] solve for force constants
                    fc = - U \ f(I);

                    % save force constants
                    s_id = repelem([1:bvt.nshells],cellfun(@(x)size(x,2),bvt.W));
                    for s = 1:bvt.nshells; bvt.fc{s} = double(fc(s_id==s)).'; end
                    
                case 2
                    % get the second order harmonic and third order anharmonic
                    % force constants all in one go without first subtrating the
                    % harmonic forces.
                    
                    % get U matrix and indexing I
                    [Uh,Ih] = get_bvk_U_matrix(bvk,pp,u);
                    % get U matrix and indexing I
                    [Ua,Ia] = get_bvt_U_matrix(bvt,pt,u);
                    % check that indicies match
                    if maxabs_(Ih-Ia)>am_lib.eps
                        error('indicies do not match');
                        % this case has not yet been coded...
                    else
                        I = Ih;
                    end

                    % solve for force constants
                    fc = - [ Uh, Ua ] \ f(I);

                    % save harmonic and anharmnoic force constants
                    h_id = repelem([1:bvk.nshells],cellfun(@(x)size(x,2),bvk.W));
                    a_id = repelem([1:bvt.nshells],cellfun(@(x)size(x,2),bvt.W));
                    s_id = [h_id,bvk.nshells+a_id];
                    for s = 1:bvk.nshells; bvk.fc{s} = double(fc(s_id==s            ).'); end
                    for s = 1:bvt.nshells; bvt.fc{s} = double(fc(s_id==s+bvk.nshells).'); end
            end
        end
        
        
        % electrons

        function [tb,pp]      = get_tb(pc,uc,dft,cutoff,spdf,nskips)
            
            
            import am_lib.*
            import am_dft.*
            
            % get irreducible shells
            fprintf(' ... identifying pairs'); tic;
            [tb,pp] = get_pairs(pc,uc,cutoff);
            fprintf(' (%.f secs)\n',toc);
            
            % [cart] print shell results
            print_pairs(uc,pp)
            
            % tight binding model
            fprintf(' ... solving for symbolic matrix elements and hamiltonian'); tic;
            tb = get_tb_model(tb,pp,uc,spdf);
            fprintf(' (%.f secs)\n',toc);
            
            % get tight binding matrix elements
            fprintf(' ... solving for tight binding matrix elements '); tic;
            tb = get_tb_matrix_elements(tb,dft,nskips);
            fprintf('(%.f secs)\n',toc);

        end

        function [tb]         = get_tb_model(ip,pp,uc,spdf)
            % set oribtals per irreducible atom: spdf = {'d','p'};
            % may wish to do this to set it per species: x={'p','d'}; spdf={x{ic.species}};
            import am_lib.*
            import am_dft.*

            % set sym digits
            digits(10);

            % initialize irreducible atom properties: for each irreducible atom,
            % set azimuthal quantum numbers J{:}, symmetries D{:}, and parity-transpose F{:} 
            [J,D,F] = get_tb_model_initialize_atoms(spdf,pp.Q{1}(1:3,1:3,:));

            % primitive cell atoms define hamiltonian blocks dimensions and start/end sections
            p2i=uc.u2i(uc.p2u); for p=[1:pp.pc_natoms]; d(p)=sum(J{p2i(p)}*2+1); end; E=cumsum(d); S=E-d+1; nbands=E(end);

            % get form of force constants for irreducible prototypical bonds
            for p = 1:ip.nshells
                % get indicies
                x = ip.xy(1,p); i = uc.u2i(x); m = uc.u2p(x); dm = d(m);
                y = ip.xy(2,p); j = uc.u2i(y); n = pp.u2p(y); dn = d(n);

                % use stabilzer group to determine crystallographic symmetry relations; A*B*C' equals kron(C,A)*B(:)
                W = sum(kron_( D{j}(:,:,ip.s_ck(:,p)) , D{i}(:,:,ip.s_ck(:,p)) ) - eye(dm*dn),3);

                % partity transpose 
                if (i==j); W = W + F{i}-eye(dm*dn); end
                
                % get linearly-independent nullspace and normalize to first nonzero element
                W=real(null(W)); W=frref_(W.').'; W(abs(W)<am_lib.eps)=0; W(abs(W-1)<am_lib.eps)=1; W=W./accessc_(W,findrow_(W.').');

                % define parameters
                c = sym(sprintf('c%02i_%%d%%d',p),[dm,dn],'real'); c = c(findrow_(W.'));

                % get symmetry adapted force constants
                vsk = reshape( sym(W,'d')*c(:), [dm,dn]);

                % save important stuff (sort W to be in line with c, matlabFunction sorts D variables)
                [sav.c{p},n] = sort(c(:).'); sav.W{p} = W(:,n); sav.vsk{p} = vsk;
            end

            % create bvk structure
            tb_ = @(pp,ip,sav,nbands) struct('units','cart','bas',pp.bas2pc*pp.bas, ...
                'symb',{pp.symb},'mass',pp.mass,'species',pp.species(pp.p2u),'cutoff',pp.cutoff,'natoms',pp.pc_natoms,...
                'nbands',nbands,'nshells',size(sav.W,2),'W',{sav.W},'vsk',{sav.vsk},'d',ip.d,'v',ip.v,'xy',ip.xy);
            tb = tb_(pp,ip,sav,nbands);

            % define function to get bond vector
            vec_ = @(xy) uc2ws(uc.bas*(uc.tau(:,xy(2,:))-uc.tau(:,xy(1,:))),pp.bas);
            
            % construct symbolic hamiltonian matrix
            H=sym(zeros(tb.nbands)); kvec=sym('k%d',[3,1],'real');
            for p = 1:pp.pc_natoms
            for u = 1:pp.npairs(p)
                % get indicies:
                %    ir = irreducible shell index
                %    iq = symmetry which takes ir -> orbit
                %    x,y=        unit cell atomic indicies
                %    i,j= irreducible cell atomic indicies
                %    m,n=   primitive cell atomic indicies
                %  mp,np= vector spanning the part of the hamiltonian
                %         corresponding to primitive atoms m and n
                ir= pp.i{p}(u); iq = pp.iq{p}(u);
                x = pp.c{p}(1); y = pp.o{p}(u,1); xy = [x;y];
                i = uc.u2i(x); m = uc.u2p(x); mp = S(m):E(m);
                j = uc.u2i(y); n = pp.u2p(y); np = S(n):E(n);

                % rotate force constants and bond vector from irrep to bond
                rij = vec_(xy); rij(abs(rij)<am_lib.eps) = 0;
                vsk =  sym(D{i}(:,:,iq)) * permute(tb.vsk{ir},pp.Q{2}(:,iq)) * sym(D{j}(:,:,iq))';
                
                % build hamiltonian matrix
                H(mp,np) = H(mp,np) + vsk .* exp(sym(2i*pi * rij(:).','d') * kvec(:) );
            end
            end
            
            % simplify (speeds evaluation up significantly later)
            H = simplify(rewrite(H,'cos'),'steps',20);
            
            % attach symbolic dynamical matrix to bvk
            tb.H = matlabFunction(H);
        end

        function [tb]         = get_tb_matrix_elements(tb,dft,nskips)
            % nskips : number of dft bands to skip (e.g. 5)

            import am_lib.*
            import am_dft.*
            
            % copy number of bands to skip
            tb.nskips = nskips;

            % fit neighbor parameter at high symmetry points using poor man's simulated anneal
            d4fc = repelem(tb.d,cellfun(@(x)size(x,2),tb.W)); nfcs=numel(d4fc); x=zeros(1,nfcs); 
            d=unique(rnd_(d4fc)); d=conv([d,Inf],[1 1]/2,'valid'); nds = numel(d); r_best = Inf;
            
            % set simulated annealing temeprature and optimization options
            kT = 20; kT_decay_ = @(kT,i) kT .* exp(-i/10); rand_ = @(x) (0.5-rand(size(x))).*abs(x./max(x));
            opts = optimoptions('lsqnonlin','Display','None','MaxIter',7);
            
            % select bands
            bnd_id = [1:tb.nbands]+tb.nskips;
            
            % define cost function
            kpt_id = 1:max(round(dft.nks/20),1):dft.nks;
            % kpt_id = [1:dft.nks];
            cost_ = @(x) dft.E(bnd_id,kpt_id) - eval_energies_(tb,x,dft.k(:,kpt_id));

            % poor man's simulated annealing: loop over distances, incorporating each shell at a time
            % it appears that ignoring the loop over distance is better, at least for cases with small pair cutoffs
            for j = nds % 2:nds
                for i = 1:30
                    % simulated annealing
                    if i ~= 1; x = x_best + rand_(x_best) * kT_decay_(kT,i); end

                    % optimize
                    [x,r] = lsqnonlin_(cost_,x,[d4fc>d(j)],[],[],opts);

                    % save r_best parameter
                    if r < r_best; r_best = r; x_best = x; 
                        % plot band structure (quick and dirty)
                        plot([1:dft.nks], eval_energies_(tb,x,dft.k),'-k',...
                             [1:dft.nks], dft.E(bnd_id,:),':r');
                        set(gca,'XTick',[]); axis tight; grid on; 
                        ylabel('Energy E'); xlabel('Wavevector k'); drawnow;
                    end
                end
            end

            % redefine cost function on all kpoints
            kpt_id = [1:dft.nks];
            cost_ = @(x) dft.E(bnd_id,kpt_id) - eval_energies_(tb,x,dft.k(:,kpt_id));
            
            % final pass with all parameters and all kpoints
            [x,~] = lsqnonlin_(cost_,x,false(1,nfcs),[],[],opts);

            % save refined matrix elements and conform to bvk
            for i = [1:tb.nshells]; d(i)=size(tb.W{i},2); end; Evsk=cumsum(d); Svsk=Evsk-d+1;
            for i = [1:tb.nshells]; tb.vsk{i} = x(Svsk(i):Evsk(i)); end

        end

        function plot_tb_vs_dft(tb,dft)
            
            import am_lib.*
            import am_dft.*
            
            % plot correlation
            plotcorr_( flatten_( dft.E([1:tb.nbands]+tb.nskips,:)        ) , ...
                       flatten_( eval_energies_(tb,[tb.vsk{:}],dft.k) ) );
            xlabel('dft energies [eV]'); ylabel('tb energies [eV]');  
        end

        
        % pairs and triplets

        function [ip,pp]      = get_pairs(pc,uc,cutoff)
            %
            % IMPORTANT: 
            %    clc; X=PM(A(3,:),:)
            %    XRef=repmat(X(1,:),size(X,1),1); 
            %    Qi=findrow_(XRef==X(:,E));
            %    % These last two are equivalent:
            %    accessc_(X.',MT(:,I(Qi))).'-XRef %   Qi  takes orbits to prototype
            %    accessc_(XRef.',MT(:,Qi)).'-X    % I(Qi) takes prototype to orbit
            %
            % Q: How is pp organized?
            % A: Check it out with the code:
            %         
            %         % select a primitive atom
            %         m  = 1; 
            %         % construct array containing positions of pairs for each equivalent atom
            %         X  = reshape(uc.tau(:,pp.o{m}(:,:)),[3,size(pp.o{m})]);
            %         % get reference atom in pair
            %         X0 = permute(uc.tau(:,pp.c{m}),[1,3,2]);
            %         % subtract reference and convert to [cart]
            %         U  = matmul_(uc.bas , X-X0 );
            %         % get rotations Rq which take bond -> irrep
            %         Rq = pp.Q{1}(1:3,1:3,pp.q{m});
            %         % apply rotations
            %         UR = permute(matmul_(Rq,permute(U,[1,3,2])),[1,3,2]);
            %         % translate to wigner-seitz cell
            %         UR = reshape(uc2ws_mex(UR,uc.bas,am_lib.tiny),size(UR));

            import am_lib.*
            import am_dft.*

            % readjust cutoff based on unitcell
            cutoff = min([normc_(uc.bas)/2,cutoff]);
            
            % step 1: get pair symmetries symmetries [pc-frac]
            
                % get space symmetries
                [~,~,S] = get_symmetries(pc); nSs = size(S,3); 

                % save space symmetry combined with permutation of atomic positions as Q
                M = perms([2:-1:1]).'; Q{1} = repmat(S,1,1,size(M,2)); Q{2} = repelem(M,1,nSs); 

                % get multiplication table, list of inverse elements, and identity
                [MT,E,I]= get_multiplication_table(Q); nQs = size(MT,1);
            
            % step 2: [PM, V, ip2pp, and pp2ip]

                % get all possible pairs which have bond lengths below the cutoff 
                d_cart_ = @(dX) normc_(uc2ws(uc.bas*dX,uc.bas));
                [Y{1:2}]=ndgrid(1:uc.natoms,uc.p2u); x=[Y{2}(:),Y{1}(:)].';
                ex_ = d_cart_(uc.tau(:,x(2,:))-uc.tau(:,x(1,:)))<cutoff;

                % [pc-frac] compute action of space symmetries on pair positions
                seitz_apply_ = @(S,tau) reshape(matmul_(S(1:3,1:3,:),tau),3,[],size(S,3)) + S(1:3,4,:);
                pc_tau = uc.tau2pc*mod_(uc.tau);
                tau(:,:,:,1) = seitz_apply_(Q{1},pc_tau(:,x(1,ex_)));
                tau(:,:,:,2) = seitz_apply_(Q{1},pc_tau(:,x(2,ex_)));
                for iq = 1:nQs; tau(:,:,iq,Q{2}(:,iq)) = tau(1:3,:,iq,:); end

                % [uc-frac] shift reference atom to primitive cell and record uc index
                %     relax matching criteria here by a factor of 10;
                %     solves a problem for systems with atoms are at 1/3
                %     position where one of the coordinates may be -0.6666
                %     and the other 0.3334. Applying mod takes -0.6666 to
                %     0.3334 causing a difference of 0.001.        
                G_ = @(tau) tau - mod_(tau); tau = mod_(matmul_(inv(uc.tau2pc),tau-G_(tau(:,:,:,1))));
                P1 = member_(tau(:,:,:,1)/10,uc.tau/10);
                P2 = member_(tau(:,:,:,2)/10,uc.tau/10);

                % create a unique pair label
                [V,~,V_p2i]=unique([P1(:),P2(:)],'rows'); V=V.';

                % get permutation representation (entries are unique pair indicies)
                PM = reshape(V_p2i,size(P1)); A = get_connectivity(PM); 

                % get map
                ip2pp = findrow_(A); pp2ip = [1:size(A,1)]*A;
            
            % step 3: [xy, qi, iqi]
                
                % get symmetry which takes irrep to orbit
                qi = findrow_(PM==PM(ip2pp(pp2ip),E)); iqi = I(qi);

                % get uc indicies, vectors, and stabilizers
                xy = V(:,PM(:,E)); v = uc2ws(uc.bas*(uc.tau(:,xy(2,:))-uc.tau(:,xy(1,:))),uc.bas); s_ck = [PM==PM(:,E)].';
                
                % create "irreducible" structure
                ip_ = @(uc,s_ck,xy,d,v) struct('units','cart','bas',uc.bas2pc*uc.bas, ...
                    'symb',{uc.symb},'mass',uc.mass,'natoms',numel(uc.p2u),'species',uc.species(uc.p2u),...
                    'cutoff',cutoff,'nshells',size(xy,2),'s_ck',s_ck,'xy',xy,'d',d,'v',v);
                ip = ip_(uc,s_ck(:,ip2pp),xy(:,ip2pp),normc_(v(:,ip2pp)),v(:,ip2pp));
            
            % step 4: [c_id, o_id, i_id, q_id]
            
                mn = uc.u2p(xy); pc_natoms = numel(uc.p2u);
                for m = 1:pc_natoms
                    % record unit cell atoms of primitive type m
                    c_id{m} = find(uc.u2p==m); ncenters = numel(c_id{m});
                    % count number of orbits involving primitive cell atom m
                    npairs = sum(mn(1,:)==m);
                    % allocate space
                    o_id{m} = zeros(npairs,ncenters);
                    i_id{m} = zeros(npairs,1);
                    q_id{m} = zeros(npairs,1);
                   iq_id{m} = zeros(npairs,1);
                    % loop over centers
                    for n = 1:ncenters
                        % [uc-frac] find the closest primitive lattice vector to atom n
                        G = uc.tau(:,c_id{m}(n))-uc.tau(:,uc.p2u(m));
                        % [uc-frac] shift atom n to the primitive cell
                        tau = mod_(uc.tau - G); 
                        % find orbits around c_id{m}(n)
                        ex_ = member_(uc.tau(:,xy(1,:)),tau).'==c_id{m}(n);
                        % record uc id for the pairing atom for each orbit
                        o_id{m}(:,n) = member_(uc.tau(:,xy(2,ex_)),tau).';
                    end
                    % irreducible pair index (independent of n)
                    i_id{m}(:) = pp2ip(ex_);
                    % symmetry which takes bond to irrep (independent of n)
                    q_id{m}(:) =  qi(ex_);
                   iq_id{m}(:) = iqi(ex_);
                end

            % define primitive pair saving function
            pp_ = @(uc,c_id,o_id,i_id,q_id,iq_id,Q) struct(...
                'units','cart',...
                'bas',uc.bas,'bas2pc',uc.bas2pc,'tau2pc',uc.tau2pc,...
                'symb',{uc.symb},'mass',uc.mass,'natoms',uc.natoms,'tau',uc.bas*uc.tau,'species',uc.species,...
                'u2p',uc.u2p,'u2i',uc.u2i,'p2u',uc.p2u,'i2u',uc.i2u, ...
                'cutoff',cutoff,'pc_natoms',numel(uc.p2u),...
                'npairs',cellfun(@(x)size(x,1),o_id),...
                'ncenters',cellfun(@(x)size(x,2),o_id), ...
                'c',{c_id},'o',{o_id},'i',{i_id},'q',{q_id},'iq',{iq_id},...
                'nQs',size(Q{1},3),'Q',{Q});
            
            % covert symmetries [pc-frac] -> [cart]
            sym_rebase_ = @(B,S) [[ matmul_(matmul_(B,S(1:3,1:3,:)),inv(B)), ...
                reshape(matmul_(B,S(1:3,4,:)),3,[],size(S,3))]; S(4,1:4,:)];
            Q{1} = sym_rebase_(uc.bas2pc*uc.bas,Q{1});
            
            % correct rounding errors in cart
            for i = 1:numel(Q{1}); for wdv = [0,1,.5,sqrt(3)/2]
                if abs(abs(Q{1}(i))-wdv)<am_lib.tiny; Q{1}(i)=wdv*sign(Q{1}(i)); end
            end;end
        
            % save "primitive" pairs
            pp = pp_(uc,c_id,o_id,i_id,q_id,iq_id,Q);
        end

        function [it,pt]      = get_triplets(pc,uc,cutoff)

            import am_lib.*
            import am_dft.*

            % readjust cutoff based on unitcell
            cutoff = min([normc_(uc.bas)/2,cutoff]); 

            % step 1: get pair symmetries symmetries [pc-frac]

                % get space symmetries
                [~,~,S] = get_symmetries(pc); nSs = size(S,3); 

                % save space symmetry combined with permutation of atomic positions as Q
                M = perms([1:3]).'; Q{1} = repmat(S,1,1,size(M,2)); Q{2} = repelem(M,1,nSs);

                % get multiplication table, list of inverse elements, and identity
                [MT,E,I]= get_multiplication_table(Q); nQs = size(MT,1);

            % step 2: [PM, V, ip2pp, and pt2it]

                % get all possible triplets for which every bond length is below the cutoff 
                [Y{1:3}]=ndgrid(1:uc.natoms,1:uc.natoms,uc.p2u); x=[Y{3}(:),Y{2}(:),Y{1}(:)].'; ex_=true(1,size(x,2));
                ex_(ex_) = normc_(uc2ws(uc.bas*(uc.tau(:,x(2,ex_))-uc.tau(:,x(1,ex_))),uc.bas))<cutoff;
                ex_(ex_) = normc_(uc2ws(uc.bas*(uc.tau(:,x(3,ex_))-uc.tau(:,x(2,ex_))),uc.bas))<cutoff;
                ex_(ex_) = normc_(uc2ws(uc.bas*(uc.tau(:,x(1,ex_))-uc.tau(:,x(3,ex_))),uc.bas))<cutoff;

                % [pc-frac] compute action of space symmetries on pair positions
                % NOTE: Q{1} is active and Q{2} is passive!
                seitz_apply_ = @(S,tau) reshape(matmul_(S(1:3,1:3,:),tau),3,[],size(S,3)) + S(1:3,4,:);
                pc_tau = uc.tau2pc*mod_(uc.tau);
                tau(:,:,:,1) = seitz_apply_(Q{1},pc_tau(:,x(1,ex_)));
                tau(:,:,:,2) = seitz_apply_(Q{1},pc_tau(:,x(2,ex_)));
                tau(:,:,:,3) = seitz_apply_(Q{1},pc_tau(:,x(3,ex_)));
                for iq = 1:nQs; tau(:,:,iq,Q{2}(:,iq)) = tau(1:3,:,iq,:); end

                % [uc-frac] shift reference atom to primitive cell and record uc index
                G_ = @(tau) tau - mod_(tau); tau = mod_(matmul_(inv(uc.tau2pc),tau-G_(tau(:,:,:,1))));
                P1 = member_(tau(:,:,:,1)/10,uc.tau/10);
                P2 = member_(tau(:,:,:,2)/10,uc.tau/10);
                P3 = member_(tau(:,:,:,3)/10,uc.tau/10);

                % create a unique pair label
                [V,~,V_p2i]=unique([P1(:),P2(:),P3(:)],'rows'); V=V.';

                % get permutation representation (entries are unique pair indicies)
                PM = reshape(V_p2i,size(P1)); A = get_connectivity(PM); 

                % get map
                it2pt = findrow_(A); pt2it = [1:size(A,1)]*A;

            % step 3: [xy, qi, iqi]

                % get symmetry which takes irrep to orbit
                qi = findrow_(PM==PM(it2pt(pt2it),E)); iqi = I(qi); % i=2; X=accessc_(PM(pt2it==i,:).',MT(:,qi(pt2it==i))).'

                % get uc indicies, vectors, and stabilizers
                xyz = V(:,PM(:,E)); s_ck = [PM==PM(:,E)].';

                % create "irreducible" structure
                it_ = @(uc,s_ck,xyz) struct('units','cart','bas',uc.bas2pc*uc.bas, ...
                    'cutoff',cutoff,'symb',{uc.symb},'mass',uc.mass,'natoms',numel(uc.p2u),'species',uc.species(uc.p2u),...
                    'nshells',size(xyz,2),'s_ck',s_ck,'xyz',xyz);
                it = it_(uc,s_ck(:,it2pt),xyz(:,it2pt));

            % step 4: [c_id, o_id, i_id, q_id]

                mno = uc.u2p(xyz); pc_natoms = numel(uc.p2u);
                for m = 1:pc_natoms
                    % record unit cell atoms of primitive type m
                    c_id{m} = find(uc.u2p==m); ncenters = numel(c_id{m});
                    % count number of orbits involving primitive cell atom m
                    npairs = sum(mno(1,:)==m);
                    % allocate space
                    o_id{m} = zeros(npairs,ncenters,2);
                    i_id{m} = zeros(npairs,1);
                    q_id{m} = zeros(npairs,1);
                   iq_id{m} = zeros(npairs,1);
                    % loop over centers
                    for n = 1:ncenters
                        % [uc-frac] find the closest primitive lattice vector to atom n
                        G = uc.tau(:,c_id{m}(n))-uc.tau(:,uc.p2u(m));
                        % [uc-frac] shift atom n to the primitive cell
                        tau = mod_(uc.tau - G); 
                        % find orbits around c_id{m}(n)
                        ex_ = member_(uc.tau(:,xyz(1,:)),tau).'==c_id{m}(n);
                        % record uc id for the pairing atom for each orbit
                        o_id{m}(:,n,:) = reshape(member_(uc.tau(:,xyz(2:3,ex_)),tau).',2,[]).';
                    end
                    % irreducible pair index (independent of n)
                    i_id{m}(:) = pt2it(ex_);
                    % symmetry which takes bond to irrep (independent of n)
                    q_id{m}(:) =  qi(ex_);
                   iq_id{m}(:) = iqi(ex_);
                end

            % define primitive pair saving function
            pt_ = @(uc,c_id,o_id,i_id,q_id,iq_id,Q) struct(...
                'units','cart',...
                'bas',uc.bas,'bas2pc',uc.bas2pc,'tau2pc',uc.tau2pc,...
                'symb',{uc.symb},'mass',uc.mass,'natoms',uc.natoms,'tau',uc.bas*uc.tau,'species',uc.species,...
                'u2p',uc.u2p,'u2i',uc.u2i,'p2u',uc.p2u,'i2u',uc.i2u, ...
                'cutoff',cutoff,'pc_natoms',numel(uc.p2u),...
                'npairs',cellfun(@(x)size(x,1),o_id),...
                'ncenters',cellfun(@(x)size(x,2),o_id), ...
                'c',{c_id},'o',{o_id},'i',{i_id},'q',{q_id},'iq',{iq_id},...
                'nQs',size(Q{1},3),'Q',{Q});

            % covert symmetries [pc-frac] -> [cart]
            sym_rebase_ = @(B,S) [[ matmul_(matmul_(B,S(1:3,1:3,:)),inv(B)), ...
                reshape(matmul_(B,S(1:3,4,:)),3,[],size(S,3))]; S(4,1:4,:)];
            Q{1} = sym_rebase_(uc.bas2pc*uc.bas,Q{1});
            % correct rounding errors in cart
            for i = 1:numel(Q{1}); for wdv = [0,1,.5,sqrt(3)/2]
                if abs(abs(Q{1}(i))-wdv)<am_lib.tiny; Q{1}(i)=wdv*sign(Q{1}(i)); end
            end;end

            % save "primitive" pairs
            pt = pt_(uc,c_id,o_id,i_id,q_id,iq_id,Q);

        end
        
        function print_pairs(uc,pp)
            
            import am_lib.*
            import am_dft.*
            
            vec_ = @(xy) uc2ws(uc.bas*(uc.tau(:,xy(2,:))-uc.tau(:,xy(1,:))),uc.bas); Z=[]; 
            bar_ = @(x) repmat('-',[1,x]); fprintf('%s primitive shells %s\n', bar_(30), bar_(30) );
            for m = 1:pp.pc_natoms
                Y=[]; ex_ = uniquemask_(pp.i{m});
                fprintf('atom %i: %i shells\n', m, sum(ex_));
                fprintf('%-10s    %-30s   %-4s   %-7s   %-7s   %-4s\n', 'd [cart]','bond [cart]','#','ic(i,j)','pc(m,n)','irr.'); 
                fprintf('%-10s    %-30s   %-4s   %-7s   %-7s   %-4s\n', bar_(10),bar_(30),bar_(4),bar_(7),bar_(7),bar_(4));
            for i = 1:pp.npairs(m)
                if ex_(i)
                    % record basic info
                    xyp = [pp.c{m}(1);pp.o{m}(i,1)]; % uc indicies
                    mn  = uc.u2p(xyp).'; % pc indicies
                    ij  = uc.u2i(xyp).'; % ic indicies
                    v   = vec_(xyp); d = normc_(v);
                    ir  = pp.i{m}(i); % irreducible index 
                    w   = sum(pp.i{m}==ir); % number of points in orbit
                    % save stuff [ d(1), r(2,3,4), w(5), ij(6,7), mn(8,9), irres(10)]
                    Y = [Y,[d;v;w;ij;mn;ir]];
                end
            end
                fprintf('%10.5f  %10.5f %10.5f %10.5f   %4i   %-3i-%3i   %-3i-%3i   %4i\n', Y(:,rankc_(Y(1,:))) ); fprintf('\n');
                Z=[Z,Y];
            end
            w = accumarray(Z(end,:).',Z(5,:).',[],@sum); Z = Z(:,uniquemask_(Z(end,:).')); Z(5,:) = w; 
            fprintf('%s irreducible shells %s\n', bar_(29), bar_(29) );
            fprintf('%-10s    %-30s   %-4s   %-7s   %-7s   %-4s\n', 'd [cart]','bond [cart]','#','ic(i,j)','pc(m,n)','irr.'); 
            fprintf('%-10s    %-30s   %-4s   %-7s   %-7s   %-4s\n', bar_(10),bar_(30),bar_(4),bar_(7),bar_(7),bar_(4));
            fprintf('%10.5f  %10.5f %10.5f %10.5f   %4i   %-3i-%3i   %-3i-%3i   %4i\n', Z(:,rankc_(Z(1,:))) );
        end

        function print_triplets(uc,pt)
            
            import am_lib.*
            import am_dft.*
            
            vec_ = @(xy) uc2ws(pt.tau(:,xy(2,:)) - pt.tau(:,xy(1,:)),pt.bas); Z=[]; 
            bar_ = @(x) repmat('-',[1,x]); fprintf('%s primitive shells %s\n', bar_(30), bar_(30) );
            for m = 1:pt.pc_natoms
                Y=[]; ex_ = uniquemask_(pt.i{m});
                fprintf('atom %i: %i shells\n', m, sum(ex_));
                fprintf('  %-30s    %-30s   %-4s   %-11s   %-11s   %-4s\n','tau_1 [cart]','tau_2 [cart]','#','ic(i,j,k)','pc(m,n,o)','irr.'); 
                fprintf('  %-30s    %-30s   %-4s   %-11s   %-11s   %-4s\n',      bar_(30),      bar_(30),bar_(4),bar_(11),bar_(11),bar_(4));
            for i = 1:pt.npairs(m)
                if ex_(i)
                    % record basic info
                    xyzp = [pt.c{m}(1);pt.o{m}(i,1,1);pt.o{m}(i,1,2)]; % uc indicies
                    mno  = uc.u2p(xyzp).'; % pc indicies
                    ijk  = uc.u2i(xyzp).'; % ic indicies
                    ir   = pt.i{m}(i); % irreducible index 
                    v1   = vec_(xyzp([1,2]));
                    v2   = vec_(xyzp([1,3]));
                    d    = normc_(v1)+normc_(v2);
                    w    = sum(pt.i{m}==ir); % number of points in orbit
                    % save stuff [ d(1), r(2,3,4), r(5,6,7), w(8), ij(9,10), mn(11,12), irres(13)]
                    Y = [Y,[d;v1;v2;w;ijk;mno;ir]];
                end
            end
                fprintf('%10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f   %4i   %3i-%3i-%3i   %3i-%3i-%3i   %4i\n', Y(2:end,rankc_(Y(1,:)))); fprintf('\n');
                Z=[Z,Y];
            end

            w = accumarray(Z(end,:).',Z(8,:).',[],@sum); Z = Z(:,uniquemask_(Z(end,:).')); Z(8,:) = w; 
            fprintf('%s irreducible shells %s\n', bar_(29), bar_(29) );
            fprintf('  %-30s    %-30s   %-4s   %-11s   %-11s   %-4s\n','tau_1 [cart]','tau_2 [cart]','#','ic(i,j,k)','pc(m,n,o)','irr.'); 
            fprintf('  %-30s    %-30s   %-4s   %-11s   %-11s   %-4s\n',      bar_(30),      bar_(30),bar_(4),bar_(11),bar_(11),bar_(4));
            fprintf('%10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f   %4i   %3i-%3i-%3i   %3i-%3i-%3i   %4i\n', Z(2:end,rankc_(Z(1,:)))); 
        end
    end
        
    
    % aux library
    
    methods (Static)
        
        % aux unit cells
        
        function [pc,p2u,u2p] = get_primitive_cell(uc)
            % NOTE: saves p2u entries which share a common closest
            % primitive lattice vector, not just the first primitive atoms
            % produced by the matrix A. When building shells, this property
            % is exploited.

            import am_lib.*
            import am_dft.*
            
            % build permutation matrix for atoms related by translations
            T = get_symmetries(uc); nTs=size(T,2); PM=zeros(uc.natoms,nTs);
            for i = [1:nTs]; PM(:,i)=rankc_( [mod_(uc.tau(:,:,1)+T(1:3,i));uc.species] ); end

            % construct a sparse binary representation 
            A=zeros(uc.natoms); A(sub2ind([1,1]*uc.natoms,repmat([1:uc.natoms].',nTs,1),PM(:)))=1; A=frref_(A); A=A(~all(A==0,2),:);

            % set basis (the three smallest vectors which preserve periodic boundary conditions)
            inds=[0,0,0];
            for j = 1:nTs; if any(abs(T(:,j))>am_lib.eps); inds(1)=j; break; end; end
            for j = 1:nTs; if any(abs(cross(T(:,2),T(:,j)))>am_lib.eps); inds(2)=j; break; end; end
            for j = 1:nTs; inds(3)=j; if abs(det(T(:,inds))+eye(3)*eps) > am_lib.eps; break; end; end
            B=T(:,inds); if det(B)<0; B=fliplr(B); end
            
            % set identifiers (see NOTE: cannot simply using p2u = findrow_(A)!)
            p2u = member_(mod_(B*mod_(B\uc.tau(:,findrow_(A),1))),mod_(uc.tau(:,:,1))).'; u2p = ([1:size(A,1)]*A);

            % define primitive cell creation function and make structure
            pc_ = @(uc,B,p2u) struct('units','frac','bas',uc.bas*B, ...
                'symb',{uc.symb},'mass',uc.mass,'nspecies',sum(unique(uc.species(p2u)).'==uc.species(p2u),2).', ...
                'natoms',numel(p2u),'tau',mod_(B\uc.tau(:,p2u)),'species',uc.species(p2u) );
            pc = pc_(uc,B,p2u);
        end

        function [ic,i2p,p2i] = get_irreducible_cell(pc)
            % idenitifes irreducible atoms
            
            import am_lib.*
            import am_dft.*
            
            % get seitz matrices
            [~,~,S] = get_symmetries(pc);

            % define function to apply symmetries to position vectors
            seitz_apply_ = @(S,tau) mod_(reshape(matmul_(S(1:3,1:3,:),tau),3,[],size(S,3)) + S(1:3,4,:));

            % get permutation matrix and construct a sparse binary representation 
            PM = member_(seitz_apply_(S,pc.tau),pc.tau); A = get_connectivity(PM);

            % set identifiers
            i2p = round(findrow_(A)).'; p2i = round(([1:size(A,1)]*A));
            
            % define irreducible cell creation function and make structure
            ic_ = @(uc,i2p) struct('units','frac','bas',uc.bas, ...
                'symb',{uc.symb},'mass',uc.mass,'nspecies',sum(unique(i2p).'==i2p,2).', ...
                'natoms',numel(i2p),'tau',uc.tau(1:3,i2p),'species',uc.species(i2p));
            ic = ic_(pc,i2p);
        end

        function [uc,inds]    = match_cell(uc,uc_ref)
            
            import am_lib.*
            import am_dft.*
            
            % find closest atom
            [~,inds] = min(reshape(normc_( ...
                mod_(osum_(-uc.tau(:,:,1),uc_ref.tau(:,:,1),2)+.5)-.5 ),...
                                                    uc.natoms,uc.natoms) );

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
        

        % aux brillouin zones
        
        function [fbz]         = get_fbz(pc,n)

            import am_lib.*
            import am_dft.*
            
            % check
            if any(mod(n,1)~=0); error('n must be integers'); end
            if numel(n)~=3; error('n must be three integers'); end
                
            % generate primitive lattice vectors
            Q_ = @(i) [0:(n(i)-1)]./n(i); [Y{1:3}]=ndgrid(Q_(1),Q_(2),Q_(3)); k=reshape(cat(3+1,Y{:}),[],3).';

            % define irreducible cell creation function and make structure
            fbz_ = @(uc,n,k) struct('units','frac-recp','recbas',inv(uc.bas).',...
                'n',n,'nks',size(k,2),'k',k,'w',ones([1,size(k,2)]));
            fbz = fbz_(pc,n,k);

        end

        function [ibz,i2f,f2i] = get_ibz(fbz,pc)
            
            import am_lib.*
            import am_dft.*
            
            % get point symmetries [real-frac --> rec-frac] by transposing R
            [~,~,~,R] = get_symmetries(pc); R = permute(R,[2,1,3]); 

            % build permutation matrix for kpoints related by point symmetries
            PM = member_(mod_(matmul_(R,fbz.k)),fbz.k); A = get_connectivity(PM);

            % set identifiers
            i2f = round(findrow_(A)).'; f2i = round(([1:size(A,1)]*A)); w=sum(A,2).';
            if abs(sum(w)-prod(fbz.n))>am_lib.eps; error('mismatch: kpoint mesh and point group symmetry'); end

            % get irreducible tetrahedra
            [tet,~,tet_f2i] = unique(sort(f2i(get_tetrahedra(fbz.recbas,fbz.n))).','rows'); tet=tet.'; tetw = hist(tet_f2i,[1:size(tet,2)].'-.5);

            % define irreducible cell creation function and make structure
            ibz_ = @(fbz,i2f,w,tet,tetw) struct('units','frac-recp','recbas',fbz.recbas,...
                'n',fbz.n,'nks',numel(i2f),'k',fbz.k(:,i2f),'w',w,'ntets',size(tet,2),'tet',tet,'tetw',tetw);
            ibz = ibz_(fbz,i2f,w,tet,tetw);
        end

        function tet           = get_tetrahedra(recbas,n)
            % divide mesh into boxes
            box = am_lib.grid2box(n); nboxes = size(box,2);
            % divide a single box into six tetrahedron
            tetrahedron = am_lib.box2tetrahedron(recbas);
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
        end

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
        
        function [y]           = ibz2fbz(fbz,ibz,x)
            % copy ibz values on fbz grid
            %
            %       x [ n , ibz.nks ] ==> y [ n , fbz.nks ]
            % 
            % n can be any number of rows
            %
            y = zeros(size(x,1),fbz.nks);
            for i = [1:ibz.nks]; y(:,fbz.f2i==i) = repmat( x(:,i) , [1,sum(fbz.f2i==i)] ); end

        end
        
        function [qqq]         = get_qqq(fbz,ibz)
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
        
        
        % aux phonons (harmonic)

        function [f]   = get_bvk_forces(bvk,pp,u,algo)
            %
            % Get forces f from the displacement u.
            %
            % Selecting a METHOD is optional (defaults to METHOD 2):
            %
            % METHOD 1: Build U matrix and compute forces using:
            %           F [3m * 1] = - U [3m * nFcs] * FC [ nFCs * 1] for n pairs, m atoms
            %
            % METHOD 2: Compute F directly from tensor double dot products using:
            %           F [3 * m] = - FC [3 * 3n] * U [3n * m]: 
            %           for n pairs, m atoms (3 spans displacements)
            %
            % Get [cart] displacement using, for example:
            %     u = matmul_( md.bas, mod_( md.tau-uc.tau +.5 )-.5 );
            %
            % Compare forces computed using the two methos:
            %   u = matmul_( md.bas, mod_( md.tau-uc.tau +.5 )-.5 );
            %   tic; f1 = get_bvk_forces(bvk,pp,u,1); toc
            %   tic; f2 = get_bvk_forces(bvk,pp,u,2); toc
            %   plot(f1(:),f2(:),'.')
            %

            import am_lib.*
            import am_dft.*            

            if nargin ~= 4; algo=2; end
            
            switch algo
                case 1
                    % get U matrix and indexing
                    [U,I] = get_bvk_U_matrix(bvk,pp,u);
                    % solve for forces f [3m * 1] = - U [3m * nFcs] * FC [ nFCs * 1]: n pairs, m atoms
                    f(I) = - U * [bvk.fc{:}].';
                    % rearrange forces 
                    f = reshape(f,size(u));
                case 2
                    % get sizes
                    [~,natoms,nsteps] = size(u);
                    % build force constants
                    for m = 1:bvk.natoms
                        phi{m} = zeros(3,3*pp.npairs(m));
                    for j = 1:pp.npairs(m)
                        % get indicies
                        i = pp.i{m}(j); iq = pp.iq{m}(j);
                        % get irrep force constant indicies
                        iphi = reshape(bvk.W{i}*bvk.fc{i}(:),3,3);
                        % rotate force constants from irrep to orbit
                        phi{m}(1:3,[1:3]+3*(j-1)) = pp.Q{1}(1:3,1:3,iq) * permute(iphi,pp.Q{2}(:,iq)) * pp.Q{1}(1:3,1:3,iq).';
                    end
                    end
                    % compute forces on every atom at every step
                    f = zeros(3,natoms,nsteps);
                    for j = 1:nsteps
                        for m = 1:bvk.natoms
                            f(1:3,pp.c{m},j) = - phi{m} * reshape(u(:,pp.o{m},j), size(pp.o{m}).*[3,1]);
                        end
                    end
            end
        end
        
        function [U,I] = get_bvk_U_matrix(bvk,pp,u)
            
            import am_lib.*
            import am_dft.*

            % the Z matrix factors out force constants leaving displacements
            Z = [1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;
                 0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0;
                 0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1].';

            % initialize arrays 
            nsteps = size(u,3);
            nFCs = sum(cellfun(@(x)size(x,2),bvk.W));
            U = zeros(3*sum(nsteps*pp.ncenters),nFCs);
            I = zeros(3*sum(nsteps*pp.ncenters),1);
            X = reshape(1:numel(u),size(u));

            % record which shell FCs belong to
            s_id = repelem([1:bvk.nshells],cellfun(@(x)size(x,2),bvk.W));
            m_id = repelem([1:bvk.natoms],3*nsteps*pp.ncenters);

            % F [3m * 1] = - U [3m * nFcs] * FC [ nFCs * 1]: n pairs, m atoms ==> FC = - U \ F
            for m = 1:bvk.natoms
                % get inds to properly reordering of force constants
                I(m_id==m) = reshape(X(:,pp.c{m},:),[],1);
            for s = 1:bvk.nshells
                ex_    = [pp.i{m}==s]; 
                npairs = sum(ex_);
                nFCs   = size(bvk.W{s},2);
                natoms = pp.ncenters(m)*nsteps; 

                % define array reshape functions
               ushp_ = @(X) reshape(X,3,  npairs,natoms);
                shp_ = @(X) reshape(X,3,9,npairs,natoms);
               fshp_ = @(X) reshape(X,3,    nFCs,natoms);

               % get rotation matrices
                R  = pp.Q{1}(1:3,1:3, pp.q{m}(ex_)); 
                iR = pp.Q{1}(1:3,1:3,pp.iq{m}(ex_));

                % get working matrix UW
                UW = ushp_(u(:,pp.o{m}(ex_,:),:));
                UW = matmul_(R,permute(UW,[1,4,2,3]));
                UW = shp_(matmul_(Z,UW));
                UW = shp_(matmul_(iR,UW));
                UW = fshp_(matmul_(sum(UW,3),bvk.W{s}));

                % construct U
                U(m_id==m,s_id==s) = reshape(permute(UW,[1,3,2]), 3*natoms, nFCs );
            end
            end
        end
        
        function [q2u] = expand_bvk_eigenvectors(bvk,uc,bz)
            % Expand primitive cell eigenvectors onto the unit cell; bz
            % must be the full mesh with the same dimensions as the
            % super unitcell (relative to that of the primitive cell).
            %  
            % Use the expanded eigenvectors to convert normal phonon coordinates to
            % displacements and velocities (vectorized Wallace Eq. 10.41, p. 113): 
            % U = q2u [ 1:3 * uc.natoms , bvk.nbranches * ibz.nks ]  * q_ks
            %
            % Note: Hermitian matrices eigendecomposed as A = U E U* for 
            % a unitary matrix U and eigenvalues E -- meaning: 
            %       q2u'*q2u = q2u*q2u' = Identity 
            %
            % Q: Are the eigenvectors of a hermitian matrix unitary? 
            % A: It appears to be so. Check with:
            %
            %        for i = 1:bz.nks; 
            %               uni(i) = max(max(abs(inv(bz.U(:,:,i))-bz.U(:,:,i)'))); 
            %        end; semilogy(uni)
            %
            % Q: Does U diagonalize D? 
            % A: Yes, U does diagonalize D when D is constructed properly including
            %    1/sqrt(M) factors! Otherwise it will leave it in block diagonal form
            %    with each block spanning the dimensions of the primitive cell basis.
            %    Check with:
            %
            %         n=[5;5;5];
            %         [bz,~] = get_zones(pc,n,''); bz = get_bvk_dispersion(bvk,bz); 
            %         U   = expand_bvk_eigenvectors(bvk,uc,bz); N = normc_(real(U));
            %         q2u = real(U)  ./ N   ./ repelem(sqrt(uc.mass(uc.species)).',3,1);
            %         u2q = real(U)' ./ N.' .* repelem(sqrt(uc.mass(uc.species)).',3,1).';
            % 
            %         M = repelem(uc.mass(uc.species),1,3);
            %         D = zeros(3*uc.natoms,3*uc.natoms);
            %         for m = 1:pp.pc_natoms
            %             for i = 1:pp.ncenters(m)
            %             for j = 1:pp.npairs(m)
            %                 x  = pp.c{m}(i);   xp = [1:3]+3*(x-1); ir = pp.i{m}(j);
            %                 y  = pp.o{m}(j,i); yp = [1:3]+3*(y-1); iq = pp.iq{m}(j);
            % 
            %                 iphi = reshape(bvk.W{ir}*bvk.fc{ir}(:),3,3);
            % 
            %                 D(xp,yp) = D(xp,yp) + ...
            %                     pp.Q{1}(1:3,1:3,iq) * permute(iphi,pp.Q{2}(:,iq)) * pp.Q{1}(1:3,1:3,iq)';
            %             end
            %             end
            %         end
            %         D = diag(1./sqrt(M)) * D * diag(1./sqrt(M));
            %         spy(abs(U'*D*U)>1E-8); view(2);  axis tight;
            %         spy(abs(u2q*D*q2u)>1E-4); view(2);  axis tight;
            %
            % Q: Are the eigenvalues obtained from the diagonalization the same?
            % A: The values are the same as that obtained by diagonalizing the
            %    primitive cell dynamical matrix. Eigenvenvalues obtained from the unit
            %    cell dynamical matrix are the same but sorted in a different order.
            %    Check with:
            % 
            %       plot(bz.hw(:)./am_lib.units_eV,real(sqrt(real(diag(U'*D*U)))),'.');
            %
            % Q: Does real(U) span a complete basis? 
            % A: No. While U does span a complete basis, the rank of real(U) is reduced and
            %    the null space of real(U) is increased by the same amount relative to
            %    that of U, indicating that real(U) vectors overlap on a smaller manifold.
            %
            
            import am_lib.*
            import am_dft.*

            % orthonormalize U
            % U = bz.U;
            % for i = 1:bz.nks
            %     [~,~,a]=unique(rnd_(bz.hw(:,i)));
            %     for j = 1:max(a)
            %         U(:,a==j,i)=orth(U(:,a==j,i));
            %     end
            % end
            
            % find closest primitive lattice vector
            % G_ = @(tau) tau - inv(uc.tau2pc)*mod_(uc.tau2pc*tau); 
            G_ = @(tau) tau;

            % expand primitive eigenvectors to unit cell
            % NOTE: exp in must ABSTOLUTELY be positive, i.e. +2i*pi and NOT -2i*pi
            u2p = flatten_([1:3].'+3*(uc.u2p-1));
            U = reshape(bz.U(u2p,:,:),3*uc.natoms,bvk.nbranches*bz.nks);
            E = repelem(exp(2i*pi*(uc.tau2pc*G_(uc.tau)).'*bz.k),3,bvk.nbranches);
            q2u = U .* E * sqrt(bvk.natoms/uc.natoms);
            
            % % explicit
            % q2u = zeros(uc.natoms,bvk.nbranches*bz.nks);
            % for i = 1:uc.natoms
            %     y = 0;  x = i; xp = [1:3] + 3*(x-1); 
            %     m = uc.u2p(x); mp = [1:3] + 3*(m-1);
            % for j = 1:bz.nks; for k = 1:bvk.nbranches; y = y+1;
            %     q2u(xp,y) = bz.U(mp,k,j) * exp(-2i*pi*dot(uc.tau2pc*G_(uc.tau(:,x)),bz.k(:,j))) * sqrt(bvk.natoms/uc.natoms);
            % end;end
            % end
        end

        function [bvk] = set_bvk_acoustic_sum_rules(bvk,pp)

            import am_lib.*
            import am_dft.*

            % build force constants
            phi = zeros(3,3,bvk.nshells);
            for i = 1:bvk.nshells; phi(:,:,i) = reshape(bvk.W{i}*bvk.fc{i}.',3,3); end

            % enforce acoustic sum rule
            for i = 1:bvk.nshells
                % check if it is a 0-th neighbor shell
                if bvk.xy(1,i)==bvk.xy(2,i)
                    % get index of primitive cell atom corresponding to this shell
                    m = pp.u2p(bvk.xy(1,i));
                    % get self forces
                    asr = zeros(3,3,pp.npairs(m));
                    for j = 1:pp.npairs(m)
                        % get irrep->orbit symmetry
                        iq = pp.iq{m}(j,1); q = pp.q{m}(j,1);
                        % rotate force constants from irrep to orbit
                        asr(:,:,j) = permute( pp.Q{1}(1:3,1:3,iq) * phi(:,:,pp.i{m}(j)) * pp.Q{1}(1:3,1:3,q), pp.Q{2}(:,iq) );
                    end
                    % impose asr on self-forces
                    asr = -sum(asr(:,:,pp.o{m}(:,1)~=pp.c{m}(1)),3);
                    % solve for symmetry-adapted force constants
                    A = double(bvk.W{i}); B = reshape(asr,[],1);
                    % get force constants as row vectors
                    bvk.fc{i} = reshape( A \ B , 1, []);
                end
            end
        end

        
        % aux phonons (anharmonic)
        
        function [f]   = get_bvt_forces(bvt,pt,u,algo)
            %
            % Get forces f [cart] from the displacement u [cart].
            %
            % Selecting a METHOD is optional (defaults to METHOD 2, which is ~ 10x faster):
            %
            % METHOD 1: Build U matrix and compute forces using:
            %           F [3m * 1] = - U [3m * nFcs] * FC [ nFCs * 1] for n pairs, m atoms
            %
            % METHOD 2: Compute F directly from tensor double dot products using:
            %           F [3 * m] = - FC [3 * 9n] * U [9n * m]: 
            %           for n pairs, m atoms (9 ispans  outer displacement product)
            %
            % Compare forces obtained from the two methods:
            %     u = matmul_( md.bas, mod_( md.tau-uc.tau +.5 )-.5 );
            %     tic; f1 = get_bvt_forces(bvt,pt,u,1); toc
            %     tic; f2 = get_bvt_forces(bvt,pt,u,2); toc
            %     plot(f1(:),f2(:),'.')
            %
            import am_lib.*
            import am_dft.*            

            if nargin ~= 4; algo=2; end
            
            switch algo
                case 1
                    % get U matrix and indexing
                    [U,I] = get_bvt_U_matrix(bvt,pt,u);
                    % solve for forces f [3m * 1] = - U [3m * nFcs] * FC [ nFCs * 1]: n pairs, m atoms
                    f(I) = - U * [bvt.fc{:}].';
                    % rearrange forces 
                    f = reshape(f,size(u));
                case 2
                    % get sizes
                    [~,natoms,nsteps] = size(u);
                    % build anharmonic force constant matrices for bonds
                    for m = 1:bvt.natoms
                        phi{m} = zeros(3,9*pt.npairs(m));
                    for j = 1:pt.npairs(m)
                        % get indicies
                        i = pt.i{m}(j); iq = pt.iq{m}(j);
                        % get irrep force constant indicies
                        iphi = reshape(bvt.W{i}*bvt.fc{i}(:),3,3,3);
                        % rotate force constants from irrep to orbit (NOTE: permuting iphi here
                        % does nothing because iphi already incorporates intrinsic symmetries.)
                        phi{m}(1:3,[1:9]+9*(j-1)) = reshape( kronpow_(pt.Q{1}(1:3,1:3,iq),3) * iphi(:) , [9,3] ).';
                    end
                    end
                    % compute forces on every atom at every step
                    f = zeros(3,natoms,nsteps);
                    for m = 1:bvt.natoms
                        % construct displacement outer products
                        ux = outerc_( reshape(u(:,flatten_(pt.o{m}(:,:,1)),:),3,[]) , ...
                                      reshape(u(:,flatten_(pt.o{m}(:,:,2)),:),3,[]) );
                        % reshape ux [3 x 3 x (npairs*natoms)] -> [9 x npairs x ncenters x nsteps]
                        ux = reshape(ux,9*pt.npairs(m),pt.ncenters(m),nsteps);
                        % compute forces
                        f(1:3,pt.c{m},:) = - matmul_(phi{m} , ux);
                    end
            end
        end
        
        function [U,I] = get_bvt_U_matrix(bvt,pt,u)
            
            import am_lib.*
            import am_dft.*
            
            % the Z matrix factors out force constants leaving displacements
            Z = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                 0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                 0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                 0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0;
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1].';
            % initialize arrays
            nsteps = size(u,3);
            nFCs = sum(cellfun(@(x)size(x,2),bvt.W));
            U = zeros(3*sum(nsteps*pt.ncenters),nFCs);
            I = zeros(3*sum(nsteps*pt.ncenters),1);
            X = reshape(1:numel(u),size(u));

            % record which shell FCs belong to
            s_id = repelem([1:bvt.nshells],cellfun(@(x)size(x,2),bvt.W));
            m_id = repelem([1:bvt.natoms],3*nsteps*pt.ncenters);

            % F [3m * 1] = - U [3m * nFcs] * FC [ nFCs * 1]: n pairs, m atoms ==> FC = - U \ F
            for m = 1:bvt.natoms
                % get inds to properly reordering of force constants
                I(m_id==m) = reshape(X(:,pt.c{m},:),[],1);
            for s = 1:bvt.nshells
                ex_ = [pt.i{m}==s]; 
                if any(ex_)
                    % get number of things
                    npairs = sum(ex_);
                    nFCs   = size(bvt.W{s},2);
                    natoms = pt.ncenters(m)*nsteps; 
                    % define array reshape functions
                   ushp_ = @(X) reshape(X,9,   npairs,natoms);
                    shp_ = @(X) reshape(X,3,27,npairs,natoms);
                   fshp_ = @(X) reshape(X,3,     nFCs,natoms);
                    % get rotation matrices
                    R  = pt.Q{1}(1:3,1:3, pt.q{m}(ex_)); 
                    iR = pt.Q{1}(1:3,1:3,pt.iq{m}(ex_));
                    % construct displacement outer products
                    ux = outerc_( reshape(u(:,flatten_(pt.o{m}(ex_,:,1)),:),3,[]) , ...
                                  reshape(u(:,flatten_(pt.o{m}(ex_,:,2)),:),3,[]) );
                    % reshape ux [3 x 3 x (npairs*natoms)] -> [9 x npairs x natoms]
                    ux = ushp_(ux);

                    % construct U matrix in a four-step process using a working matrix W:
                    % 1) transform outer dispacement product to irreducible orientation 
                    %    R^^2 [9 x 9] * UW [9 x 1  x npairs x natoms] = UW [9  x 1  x npairs x natoms]
                    W = matmul_(kron_(R,R),permute(ux,[1,4,2,3])); 
                    % 2) factor out force constants, leaving displacements
                    %    Z [81 x 9] * UW [9 x 1  x npairs x natoms] = UW [81 x 1  x npairs x natoms] reshaped to UW [3  x 27 x npairs x natoms]
                    W = shp_(matmul_(Z,W));
                    % 3) return displacement to bond orientation
                    %    iR [3 x 3] * UW [3 x 27 x npairs x natoms] = UW [3  x 27 x npairs x natoms]
                    W = shp_(matmul_(iR,W));
                    % 4) take into account intrinsic and crystallographic symmetries
                    %    UW [3 x 27 x npairs x natoms] * bvt.W [ 27 * nfcs ] = UW [3 x ncfs x npairs x natoms] sum over pairs -> UW [3 x ncfs x 1 x natoms] 
                    W = fshp_(matmul_(sum(W,3),bvt.W{s}));

                    % construct U
                    U(m_id==m,s_id==s) = reshape(permute(W,[1,3,2]), 3*natoms, nFCs );
                end
            end
            end
        end
        
        
        % aux electrons

        function [J,D,F] = get_tb_model_initialize_atoms(spdf,R)
            % set symmetries D{:}, and parity-transpose F{:} for each
            % irreducible atom given a list of orbitals for each
            % irreducible atom, spdf = {'sp','d'}

            import am_lib.*
            import am_dft.*

            % get symmetries
            nRs=size(R,3);

            % transform symmetries to the tight binding representation (wiger functions)
            W=cell(1,3); for j=[1:3]; W{j} = get_wigner(j,R); end

            % set orbitals J{:}, symmetries D{:}, and parity-transpose T{:} for each irreducible atom
            natoms=numel(spdf); F=cell(1,natoms);  D=cell(1,natoms);
            for i = 1:natoms
                % set orbitals
                J{i} = findrow_('spdf'==spdf{i}(:)).'-1;

                % set start and end points for J
                E=cumsum(J{i}*2+1); S=E-(J{i}*2+1)+1;

                % construct D matrix and lay the ground work construction of parity super-operator
                d = max(E); P = zeros(1,d); D{i} = zeros(d,d,nRs);
                for j = 1:length(J{i})
                    if J{i}(j)==0 % s orbital
                        D{i}(S(j):E(j),S(j):E(j),:) = 1;
                    else % p,d,f orbitals
                        D{i}(S(j):E(j),S(j):E(j),:) = W{J{i}(j)};
                    end
                    P(S(j):E(j)) = (-1).^j;
                end

                % construct parity super-operator    
                f_ = @(x) x(:); A=(P.'*P).*reshape([1:d^2],[d,d]); 
                F{i}=zeros(d^2,d^2); F{i}(sub2ind([d^2,d^2],abs(f_(A')),abs(f_(A))))=sign(f_(A'));
            end
            
            % correct rounding errors in sym (non-exauhstive)
            for i = 1:numel(D); for j = 1:numel(D{i}); for wdv = [0,1,.5,sqrt(3)/2]
                if abs(abs(D{i}(j))-wdv)<am_lib.eps; D{i}(j)=wdv*sign(real(D{i}(j))); end
            end;end;end
        end
        
        function E       = eval_energies_(tb,x,k)
            % get hamiltonians
            nks = size(k,2); E = zeros(tb.nbands,nks); recbas = inv(tb.bas).';
            for m = 1:nks
                % build input
                input = num2cell([x,(recbas*k(:,m)).']);
                % evaluate H
                E(:,m) = real(eig(tb.H(input{:})));
            end
            % sort values
            E = sort(E);
        end
        
        
        % aux pairs and triplets
        
        function [C]     = get_connectivity(PM)

            import am_lib.*
            import am_dft.*

            % binary 
            [natoms,nRs] = size(PM);

            % exclude all rows containing all zeros
            PM = PM(~all(PM==0,2),:);

            % construct sparse vectors
            m = size(PM,1); t = zeros(1,m); for i = [1:m]; t(i) = PM(i,find(PM(i,:),1)); end
            v = [ repmat(t(:),nRs,1), PM(:) ]; 

            % exlcude zeros, make symmetric, ensure diagonals, and remove repeat
            v = v(~any(v==0,2),:); v=[v;[v(:,2),v(:,1)]]; v=[v;[v(:,1),v(:,1)]]; v = unique(v,'rows');

            % construct a sparse binary representation 
            C = sparse(v(:,1),v(:,2),ones(size(v,1),1),natoms,natoms); % A = double((A'*A)~=0);

            % merge and reduce binary rep
            C = merge_(C); C(abs(C)<am_lib.tiny)=0; C(abs(C)>am_lib.tiny)=1; C=full(C(any(C~=0,2),:)); 

            % convert to logical
            C = logical(C);
        end

    end
    
    
    % unix functions and scripts
    
    methods (Static)
        
        % stand-alone
        
        function generate_scripts()
            
            import am_lib.*
            import am_dft.*
            
            % initialize counter
            i=0; fprintf('Generating scripts:\n')
            
            % scripts functions
            fname='outcar2fp.sh'; fid=fopen([fname],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); fprintf(' ... %s (succeeded)\n',fname);
            %{
            #!/bin/bash
            # Parse outcar to extract forces and displacements at each timestep.
            # Antonio Mei Nov/2014
            # Antonio Mei Jan/2017
            usage_ () {
                echo "Creates infile.force_position based on supplied outcar files."
                echo ""
                echo "Usage: $0 [-h] [-t] [-f] -o <outcar_list> -n <natoms> [-c <compress_name>]"
                echo ""
                echo "Example: $0 -f -t -o \"\$(find . -name "OUTCAR*" | grep 4.00-300)\" -n 250"
                echo ""
                echo "-h : prints this message"
                echo "-n : [REQUIRED] number of atoms in the simulation cell"
                echo "-o : [REQUIRED] list of outcar files to parse"
                echo "-t : trims the last md run (useful for removing runs which have not completed)"
                echo "-f : overwrites existing infile.force_position"
                echo "-c : compresses infile.force_position to a tar.gz file"
                echo ""
                echo "infile.force_position file contents:"
                echo "   x position   y position   z position     x force      y force      z force"
                exit 1
            }
            
            main_ () {
                # trim the last md run which may not have completed
                trim_ () { tac $1 | awk '!found && /POSITION/{found=1;next}1' | tac ; }
                # get position and forces
                get_  () { cat $2 | grep -h -A $(($1+1)) POSITION  ; }
                # cut header lines
                cut_  () { cat $1 | sed '/^--$/d' | sed '/--/d' | sed '/POSITION/d' ; }
                # compress produced infile.force_position
                compress_ () { tar -zcvf infile.force_position.tar.gz infile.force_position ; }
                #
                if ${ISFORCE}; then
                    if [ -f "./infile.force_position" ]; then
                        rm ./infile.force_position
                        printf " ... ./infile.force_position overwritten\n"
                    fi
                fi
                #
                if ${ISTRIM}; then
                    printf " ... trim:\n"
                    for F in "${FLIST}"; do
                        printf " ...     %-100s\n" "${F}"
                        trim_ ${F} | get_ ${NATOMS} | cut_ >> infile.force_position
                    done
                else
                    printf " ... batch parsing without trim\n"
                    get_ ${NATOMS} "${FLIST}" | cut_ >> infile.force_position
                fi
                #
                printf " ... infile.force_position created\n"
                #
                if ${ISCOMPRESS}; then
                    printf " ... infile.force_position.tar.gz compressed\n"
                    compress_ 
                fi
            }
            
            ISCOMPRESS=false; ISTRIM=false; ISFORCE=false;
            if (($# == 0)); then usage_; exit 1; fi
            while getopts "n:o:htfc" o; do
                case "${o}" in
                    o)  FLIST=${OPTARG} ;;
                    n)  NATOMS=${OPTARG} ;;
                    c)  ISCOMPRESS=true ;;
                    t)  ISTRIM=true ;;
                    f)  ISFORCE=true ;;
                    h)  usage_; exit 0 ;;
                    *)  usage_; exit 1 ;;
                esac
            done
            main_
            %}
            
            fname='outcar2en.sh'; fid=fopen([fname],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); fprintf(' ... %s (succeeded)\n',fname);
            %{
            #!/bin/bash
            # Preprocess outcar to remove the last run, which may not have finished
            # Antonio Mei Nov/2014
            # Antonio Mei Jan/2017
            usage_ () {
                echo "Creates infile.electron_energies based on supplied outcar files."
                echo ""
                echo "Usage: $0 [-h] [-t] [-f] -o <outcar_list> -n <nbands> [-c <compress_name>]"
                echo ""
                echo "Example: $0 -f -t -o \"\$(find . -name "OUTCAR*" | grep 4.00-300)\" -n 751"
                echo ""
                echo "-h : prints this message"
                echo "-n : [REQUIRED] number of bands"
                echo "-o : [REQUIRED] list of outcar files to parse"
                echo "-t : trims the last md run (useful for removing runs which have not completed)"
                echo "-f : overwrites existing infile.electron_energies"
                echo "-c : compresses infile.electron_energies to a tar.gz file"
                echo ""
                echo "infile.electron_energies file contents:"
                echo "   n index    En energy    fn occupation"
                exit 1
            }
            main_ () {
                # trim the last md run which may not have completed
                trim_ () { tac $1 | awk '!found && /POSITION/{found=1;next}1' | tac ; }
                # get energies
                get_  () { cat $2 | grep -h -A ${1} occupation  ; }
                # cut header lines
                cut_  () { cat $1 | sed '/^--$/d' | sed '/--/d' | sed '/occupation/d' ; }
                # compress produced infile.electron_energies
                compress_ () { tar -zcvf infile.electron_energies.tar.gz infile.electron_energies ; }
                #
                if ${ISFORCE}; then
                    if [ -f "./infile.electron_energies" ]; then
                        rm ./infile.electron_energies
                        printf " ... ./infile.electron_energies overwritten\n"
                    fi
                fi
                # 
                if ${ISTRIM}; then
                    printf " ... trim:\n"
                    for F in "${FLIST}"; do
                        printf " ...     %-100s\n" "${F}"
                        trim_ ${F} | get_ ${NBANDS} | cut_ >> infile.electron_energies
                    done
                else
                    printf " ... batch parsing without trim\n"
                    get_ ${NBANDS} "${FLIST}" | cut_ >> infile.electron_energies
                fi
                #
                awk '{ print $2 }' infile.electron_energies > infile.electron_energies.tmp && mv infile.electron_energies.tmp infile.electron_energies
                #
                printf " ... infile.electron_energies created\n"
                #
                if ${ISCOMPRESS}; then
                    printf " ... infile.electron_energies.tar.gz compressed\n"
                    compress_ 
                fi
            }
            ISCOMPRESS=false; ISTRIM=false; ISFORCE=false;
            if (($# == 0)); then usage_; exit 1; fi
            while getopts "n:o:htfc" o; do
                case "${o}" in
                    o)  FLIST=${OPTARG} ;;
                    n)  NBANDS=${OPTARG} ;;
                    c)  ISCOMPRESS=true ;;
                    t)  ISTRIM=true ;;
                    f)  ISFORCE=true ;;
                    h)  usage_; exit 0 ;;
                    *)  usage_; exit 1 ;;
                esac
            done
            main_
            %}
                
        end
        
    end    
    

    % mex functions
    
    methods (Static)
        
        function compile_mex()
            % compile mex functions
            %
            % to compile the library manually
            %
            %       fort -O3 -parallel -fpp -c am_mex_lib.f90
            %
            % to compile individual interfaces manually
            %
            %       ifort -O3 -parallel -fpp -fPIC -L/Applications/MATLAB_R2016b.app/bin/maci64
            %       -I/Applications/MATLAB_R2016b.app/extern/include -lmx -lmex -lmat -nofor_main 
            %       -bundle ./uc2ws_mex.f90 -o uc2ws_mex.mexmaci64    
            %
            
            import am_lib.*
            import am_dft.*
            
            % initialize counter
            i=0;
            
            % mex library
            flib = 'am_mex_lib'; fid=fopen([flib,'.f90'],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); 
            
% START OF FORTRAN CODE 
%{
module am_mex_lib

    implicit none

    public

    !  $$$$$$\   $$$$$$\  $$\   $$\  $$$$$$\ $$$$$$$$\  $$$$$$\  $$\   $$\ $$$$$$$$\  $$$$$$\
    ! $$  __$$\ $$  __$$\ $$$\  $$ |$$  __$$\\__$$  __|$$  __$$\ $$$\  $$ |\__$$  __|$$  __$$\
    ! $$ /  \__|$$ /  $$ |$$$$\ $$ |$$ /  \__|  $$ |   $$ /  $$ |$$$$\ $$ |   $$ |   $$ /  \__|
    ! $$ |      $$ |  $$ |$$ $$\$$ |\$$$$$$\    $$ |   $$$$$$$$ |$$ $$\$$ |   $$ |   \$$$$$$\
    ! $$ |      $$ |  $$ |$$ \$$$$ | \____$$\   $$ |   $$  __$$ |$$ \$$$$ |   $$ |    \____$$\
    ! $$ |  $$\ $$ |  $$ |$$ |\$$$ |$$\   $$ |  $$ |   $$ |  $$ |$$ |\$$$ |   $$ |   $$\   $$ |
    ! \$$$$$$  | $$$$$$  |$$ | \$$ |\$$$$$$  |  $$ |   $$ |  $$ |$$ | \$$ |   $$ |   \$$$$$$  |
    !  \______/  \______/ \__|  \__| \______/   \__|   \__|  \__|\__|  \__|   \__|    \______/

    integer    , parameter :: sp      = kind(0.0E0) !> single precision
    integer    , parameter :: dp      = kind(0.0D0) !> double precision
    real(dp)   , parameter :: tiny    =  1.0D-5
    real(dp)   , parameter :: halfpi  =  1.570796326794897_dp
    real(dp)   , parameter :: sqrtpi  =  1.772453850905516_dp
    real(dp)   , parameter :: pi      =  3.141592653589793_dp
    real(dp)   , parameter :: twopi   =  6.283185307179586_dp
    real(dp)   , parameter :: fourpi  = 12.566370614359172_dp
    complex(dp), parameter :: cmplx_i = cmplx(0,1,dp)
    complex(dp), parameter :: itwopi  = cmplx_i*twopi

    ! $$$$$$\ $$\   $$\ $$$$$$$$\ $$$$$$$$\ $$$$$$$\  $$$$$$$$\  $$$$$$\   $$$$$$\  $$$$$$$$\  $$$$$$\
    ! \_$$  _|$$$\  $$ |\__$$  __|$$  _____|$$  __$$\ $$  _____|$$  __$$\ $$  __$$\ $$  _____|$$  __$$\
    !   $$ |  $$$$\ $$ |   $$ |   $$ |      $$ |  $$ |$$ |      $$ /  $$ |$$ /  \__|$$ |      $$ /  \__|
    !   $$ |  $$ $$\$$ |   $$ |   $$$$$\    $$$$$$$  |$$$$$\    $$$$$$$$ |$$ |      $$$$$\    \$$$$$$\
    !   $$ |  $$ \$$$$ |   $$ |   $$  __|   $$  __$$< $$  __|   $$  __$$ |$$ |      $$  __|    \____$$\
    !   $$ |  $$ |\$$$ |   $$ |   $$ |      $$ |  $$ |$$ |      $$ |  $$ |$$ |  $$\ $$ |      $$\   $$ |
    ! $$$$$$\ $$ | \$$ |   $$ |   $$$$$$$$\ $$ |  $$ |$$ |      $$ |  $$ |\$$$$$$  |$$$$$$$$\ \$$$$$$  |
    ! \______|\__|  \__|   \__|   \________|\__|  \__|\__|      \__|  \__| \______/ \________| \______/

    interface foursort
        module procedure d_foursort, i_foursort
    end interface ! foursort

    interface fourrank
        module procedure d_fourrank, i_fourrank
    end interface ! fourrank

    contains

    !  $$$$$$\ $$$$$$$$\  $$$$$$\ $$$$$$$$\ $$$$$$\  $$$$$$\ $$$$$$$$\ $$$$$$\  $$$$$$\   $$$$$$\
    ! $$  __$$\\__$$  __|$$  __$$\\__$$  __|\_$$  _|$$  __$$\\__$$  __|\_$$  _|$$  __$$\ $$  __$$\
    ! $$ /  \__|  $$ |   $$ /  $$ |  $$ |     $$ |  $$ /  \__|  $$ |     $$ |  $$ /  \__|$$ /  \__|
    ! \$$$$$$\    $$ |   $$$$$$$$ |  $$ |     $$ |  \$$$$$$\    $$ |     $$ |  $$ |      \$$$$$$\
    !  \____$$\   $$ |   $$  __$$ |  $$ |     $$ |   \____$$\   $$ |     $$ |  $$ |       \____$$\
    ! $$\   $$ |  $$ |   $$ |  $$ |  $$ |     $$ |  $$\   $$ |  $$ |     $$ |  $$ |  $$\ $$\   $$ |
    ! \$$$$$$  |  $$ |   $$ |  $$ |  $$ |   $$$$$$\ \$$$$$$  |  $$ |   $$$$$$\ \$$$$$$  |\$$$$$$  |
    !  \______/   \__|   \__|  \__|  \__|   \______| \______/   \__|   \______| \______/  \______/

    pure function factorial(n) result(y)
        !
        implicit none
        !
        integer, intent(in) :: n
        real(dp) :: y
        !
        if (n.ge.1) then
            y = product([1:n])
        else
            y = 1
        endif
    end function  factorial

    pure function nchoosek(n,k) result (res)
        !
        implicit none
        !
        integer, intent(in) :: n
        integer, intent(in) :: k
        integer :: res
        !
        res=factorial(n)/(factorial(k)*factorial(n-k))
        !
    end function  nchoosek

    function      perms(n) result(PT)
        !
        implicit none
        !
        integer, intent(in) :: n
        integer, allocatable :: P(:,:)
        integer, allocatable :: PT(:,:)
        integer :: m
        !
        m = product([1:n])
        !
        allocate(P(m,n))
        !
        call permutate([1:n],P)
        !
        allocate(PT,source=transpose(P))
        !
        contains
        recursive subroutine permutate(E, P)
            !
            implicit none
            !
            integer, intent(in)  :: E(:) ! array of objects 
            integer, intent(out) :: P(:,:) ! permutations of E 
            integer :: N, Nfac, i, k, S(size(P,1)/size(E), size(E)-1) 
            N = size(E); Nfac = size(P,1); 
            do i = 1, N
              if( N>1 ) call permutate((/E(:i-1), E(i+1:)/), S) 
              forall(k=1:Nfac/N) P((i-1)*Nfac/N+k,:) = (/E(i), S(k,:)/) 
            enddo 
        end subroutine permutate 
    end function  perms

    ! $$\   $$\ $$\   $$\ $$\      $$\ $$$$$$$\  $$$$$$$$\ $$$$$$$\        $$$$$$$$\ $$\   $$\ $$$$$$$$\  $$$$$$\  $$$$$$$\ $$\     $$\
    ! $$$\  $$ |$$ |  $$ |$$$\    $$$ |$$  __$$\ $$  _____|$$  __$$\       \__$$  __|$$ |  $$ |$$  _____|$$  __$$\ $$  __$$\\$$\   $$  |
    ! $$$$\ $$ |$$ |  $$ |$$$$\  $$$$ |$$ |  $$ |$$ |      $$ |  $$ |         $$ |   $$ |  $$ |$$ |      $$ /  $$ |$$ |  $$ |\$$\ $$  /
    ! $$ $$\$$ |$$ |  $$ |$$\$$\$$ $$ |$$$$$$$\ |$$$$$\    $$$$$$$  |         $$ |   $$$$$$$$ |$$$$$\    $$ |  $$ |$$$$$$$  | \$$$$  /
    ! $$ \$$$$ |$$ |  $$ |$$ \$$$  $$ |$$  __$$\ $$  __|   $$  __$$<          $$ |   $$  __$$ |$$  __|   $$ |  $$ |$$  __$$<   \$$  /
    ! $$ |\$$$ |$$ |  $$ |$$ |\$  /$$ |$$ |  $$ |$$ |      $$ |  $$ |         $$ |   $$ |  $$ |$$ |      $$ |  $$ |$$ |  $$ |   $$ |
    ! $$ | \$$ |\$$$$$$  |$$ | \_/ $$ |$$$$$$$  |$$$$$$$$\ $$ |  $$ |         $$ |   $$ |  $$ |$$$$$$$$\  $$$$$$  |$$ |  $$ |   $$ |
    ! \__|  \__| \______/ \__|     \__|\_______/ \________|\__|  \__|         \__|   \__|  \__|\________| \______/ \__|  \__|   \__|

    pure function primes(nprimes)
        !
        ! naive approach to generating the first n prime numbers
        !
        implicit none
        !
        integer, intent(in)  :: nprimes
        integer, allocatable :: primes(:) ! array that will hold the primes
        integer :: at, found, i
        logical :: is_prime
        !
        allocate (primes(nprimes))
        !
        primes(1) = 2
        at = 2
        found = 1
        do
            is_prime = .true. ! assume prime
            do i = 1, found
                if (modulo(at,primes(i)).eq.0) then ! if divisible by any other element
                    is_prime = .false.               ! in the array, then not prime.
                    at = at + 1
                    continue
                end if
            end do
            found = found + 1
            primes(found) = at
            at = at + 1
            if (found == nprimes) then ! stop when all primes are found
                exit
            endif
        end do
        !
    end function  primes

    !  $$$$$$\   $$$$$$\  $$$$$$$\ $$$$$$$$\ $$$$$$\ $$\   $$\  $$$$$$\
    ! $$  __$$\ $$  __$$\ $$  __$$\\__$$  __|\_$$  _|$$$\  $$ |$$  __$$\
    ! $$ /  \__|$$ /  $$ |$$ |  $$ |  $$ |     $$ |  $$$$\ $$ |$$ /  \__|
    ! \$$$$$$\  $$ |  $$ |$$$$$$$  |  $$ |     $$ |  $$ $$\$$ |$$ |$$$$\
    !  \____$$\ $$ |  $$ |$$  __$$<   $$ |     $$ |  $$ \$$$$ |$$ |\_$$ |
    ! $$\   $$ |$$ |  $$ |$$ |  $$ |  $$ |     $$ |  $$ |\$$$ |$$ |  $$ |
    ! \$$$$$$  | $$$$$$  |$$ |  $$ |  $$ |   $$$$$$\ $$ | \$$ |\$$$$$$  |
    !  \______/  \______/ \__|  \__|  \__|   \______|\__|  \__| \______/


    pure function d_foursort(x) result(y)
        !
        implicit none
        !
        real(dp), intent(in) :: x(4)
        real(dp) :: y(4)
        !
        y = x
        if (y(1).gt.y(2)) then; y([1,2]) = y([2,1]); endif
        if (y(3).gt.y(4)) then; y([3,4]) = y([4,3]); endif
        if (y(1).gt.y(3)) then; y([1,3]) = y([3,1]); endif
        if (y(2).gt.y(4)) then; y([2,4]) = y([4,2]); endif
        if (y(2).gt.y(3)) then; y([2,3]) = y([3,2]); endif
        !
    end function  d_foursort

    pure function i_foursort(x) result(y)
        !
        implicit none
        !
        integer, intent(in) :: x(4)
        integer :: y(4)
        !
        y = x
        if (y(1).gt.y(2)) then; y([1,2]) = y([2,1]); endif
        if (y(3).gt.y(4)) then; y([3,4]) = y([4,3]); endif
        if (y(1).gt.y(3)) then; y([1,3]) = y([3,1]); endif
        if (y(2).gt.y(4)) then; y([2,4]) = y([4,2]); endif
        if (y(2).gt.y(3)) then; y([2,3]) = y([3,2]); endif
        !
    end function  i_foursort

    pure function d_fourrank(n) result(ind)
        ! borrowed from olle, who took it from stack overflow
        implicit none
        !
        real(dp),intent(in) :: n(4)
        integer :: ind(4)
        integer :: low1,high1,low2,high2,highest,lowest,middle1,middle2
        !    
        if ( n(1) <= n(2) ) then
            low1 = 1
            high1 = 2
        else 
            low1 = 2
            high1 = 1
        endif

        if ( n(3) <= n(4) ) then
            low2 = 3
            high2 = 4
        else
            low2 = 4
            high2 = 3
        endif

        if ( n(low1) <= n(low2) ) then
            lowest = low1
            middle1 = low2
        else
            lowest = low2
            middle1 = low1
        endif

        if ( n(high1) >= n(high2) ) then
            highest = high1
            middle2 = high2
        else
            highest = high2
            middle2 = high1
        endif

        if ( n(middle1) < n(middle2) ) then
            ind=(/lowest,middle1,middle2,highest/)
        else
            ind=(/lowest,middle2,middle1,highest/)
        endif
        !
    end function  d_fourrank

    pure function i_fourrank(n) result(ind)
        ! borrowed from olle, who took it from stack overflow
        implicit none
        !
        integer,intent(in) :: n(4)
        integer :: ind(4)
        integer :: low1,high1,low2,high2,highest,lowest,middle1,middle2
        !    
        if ( n(1) <= n(2) ) then
            low1 = 1
            high1 = 2
        else 
            low1 = 2
            high1 = 1
        endif

        if ( n(3) <= n(4) ) then
            low2 = 3
            high2 = 4
        else
            low2 = 4
            high2 = 3
        endif

        if ( n(low1) <= n(low2) ) then
            lowest = low1
            middle1 = low2
        else
            lowest = low2
            middle1 = low1
        endif

        if ( n(high1) >= n(high2) ) then
            highest = high1
            middle2 = high2
        else
            highest = high2
            middle2 = high1
        endif

        if ( n(middle1) < n(middle2) ) then
            ind=(/lowest,middle1,middle2,highest/)
        else
            ind=(/lowest,middle2,middle1,highest/)
        endif
        !
    end function  i_fourrank

    pure function d_rank(v) result(ind)
      !
      ! based on dlasrt2
      !
      implicit none
      !
      real(dp), intent(in) :: v(:)
      real(dp), allocatable :: d(:)
      integer, allocatable :: ind(:)
      integer :: select, endd, i, j, start, stkpnt, tmpind, stack(2,32), n
      real(dp) :: d1, d2, d3, dmnmx, tmp
      !
      select = 20
      !
      n = size(v)
      !
      allocate(d,source=v)
      !
      allocate(ind,source=[1:n])
      !
      if (n.le.1) return
      !
      stkpnt = 1
      stack(1,1) = 1
      stack(2,1) = n
      !
      do
         start = stack(1, stkpnt)
         endd = stack(2, stkpnt)
         stkpnt = stkpnt - 1
         if(endd-start.gt.0) then
            do i = start + 1, endd
               do j = i, start + 1, -1
                  if(d(j).lt.d(j-1)) then
                     dmnmx = d(j)
                     d(j) = d(j-1)
                     d(j-1) = dmnmx
                     tmpind = ind(j)
                     ind(j) = ind(j-1)
                     ind(j-1) = tmpind
                  else
                     exit
                  end if
               enddo
            enddo
         else if(endd-start.gt.select) then
            d1 = d(start)
            d2 = d(endd)
            i = (start+endd) / 2
            d3 = d(i)
            if(d1.lt.d2) then
               if(d3.lt.d1) then
                  dmnmx = d1
               else if(d3.lt.d2) then
                  dmnmx = d3
               else
                  dmnmx = d2
               end if
            else
               if(d3.lt.d2) then
                  dmnmx = d2
               else if(d3.lt.d1) then
                  dmnmx = d3
               else
                  dmnmx = d1
               end if
            end if
            i = start - 1
            j = endd + 1
            do
               do
                  j = j - 1
                  if (.not.(d(j).gt.dmnmx)) exit
               enddo
               do
                  i = i + 1
                  if (.not.(d(i).lt.dmnmx)) exit
               enddo
               if(i.gt.j) then
                  exit   
               else 
                  tmp = d(i)
                  d(i) = d(j)
                  d(j) = tmp
                  tmpind = ind(j)
                  ind(j) = ind(i)
                  ind(i) = tmpind
               end if
            enddo
            if(j-start.gt.endd-j-1) then
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = start
               stack(2, stkpnt) = j
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = j + 1
               stack(2, stkpnt) = endd
            else
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = j + 1
               stack(2, stkpnt) = endd
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = start
               stack(2, stkpnt) = j
            end if
         end if
         if (.not.(stkpnt.gt.0)) exit
      enddo
    end function  d_rank

    pure function i_rank(v) result(ind)
      !
      ! based on dlasrt2
      !
      implicit none
      !
      integer, intent(in) :: v(:)
      integer, allocatable :: d(:)
      integer, allocatable :: ind(:)
      integer :: select, endd, i, j, start, stkpnt, tmpind, stack(2,32), n
      integer :: d1, d2, d3, dmnmx, tmp
      !
      select = 20
      !
      n = size(v)
      !
      allocate(d,source=v)
      !
      allocate(ind,source=[1:n])
      !
      if (n.le.1) return
      !
      stkpnt = 1
      stack(1,1) = 1
      stack(2,1) = n
      !
      do
         start = stack(1, stkpnt)
         endd = stack(2, stkpnt)
         stkpnt = stkpnt - 1
         if(endd-start.gt.0) then
            do i = start + 1, endd
               do j = i, start + 1, -1
                  if(d(j).lt.d(j-1)) then
                     dmnmx = d(j)
                     d(j) = d(j-1)
                     d(j-1) = dmnmx
                     tmpind = ind(j)
                     ind(j) = ind(j-1)
                     ind(j-1) = tmpind
                  else
                     exit
                  end if
               enddo
            enddo
         else if(endd-start.gt.select) then
            d1 = d(start)
            d2 = d(endd)
            i = (start+endd) / 2
            d3 = d(i)
            if(d1.lt.d2) then
               if(d3.lt.d1) then
                  dmnmx = d1
               else if(d3.lt.d2) then
                  dmnmx = d3
               else
                  dmnmx = d2
               end if
            else
               if(d3.lt.d2) then
                  dmnmx = d2
               else if(d3.lt.d1) then
                  dmnmx = d3
               else
                  dmnmx = d1
               end if
            end if
            i = start - 1
            j = endd + 1
            do
               do
                  j = j - 1
                  if (.not.(d(j).gt.dmnmx)) exit
               enddo
               do
                  i = i + 1
                  if (.not.(d(i).lt.dmnmx)) exit
               enddo
               if(i.gt.j) then
                  exit   
               else 
                  tmp = d(i)
                  d(i) = d(j)
                  d(j) = tmp
                  tmpind = ind(j)
                  ind(j) = ind(i)
                  ind(i) = tmpind
               end if
            enddo
            if(j-start.gt.endd-j-1) then
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = start
               stack(2, stkpnt) = j
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = j + 1
               stack(2, stkpnt) = endd
            else
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = j + 1
               stack(2, stkpnt) = endd
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = start
               stack(2, stkpnt) = j
            end if
         end if
         if (.not.(stkpnt.gt.0)) exit
      enddo
    end function  i_rank

    ! $$\      $$\  $$$$$$\ $$$$$$$$\ $$\   $$\       $$$$$$$$\ $$\   $$\ $$\   $$\  $$$$$$\ $$$$$$$$\ $$$$$$\  $$$$$$\  $$\   $$\  $$$$$$\
    ! $$$\    $$$ |$$  __$$\\__$$  __|$$ |  $$ |      $$  _____|$$ |  $$ |$$$\  $$ |$$  __$$\\__$$  __|\_$$  _|$$  __$$\ $$$\  $$ |$$  __$$\
    ! $$$$\  $$$$ |$$ /  $$ |  $$ |   $$ |  $$ |      $$ |      $$ |  $$ |$$$$\ $$ |$$ /  \__|  $$ |     $$ |  $$ /  $$ |$$$$\ $$ |$$ /  \__|
    ! $$\$$\$$ $$ |$$$$$$$$ |  $$ |   $$$$$$$$ |      $$$$$\    $$ |  $$ |$$ $$\$$ |$$ |        $$ |     $$ |  $$ |  $$ |$$ $$\$$ |\$$$$$$\
    ! $$ \$$$  $$ |$$  __$$ |  $$ |   $$  __$$ |      $$  __|   $$ |  $$ |$$ \$$$$ |$$ |        $$ |     $$ |  $$ |  $$ |$$ \$$$$ | \____$$\
    ! $$ |\$  /$$ |$$ |  $$ |  $$ |   $$ |  $$ |      $$ |      $$ |  $$ |$$ |\$$$ |$$ |  $$\   $$ |     $$ |  $$ |  $$ |$$ |\$$$ |$$\   $$ |
    ! $$ | \_/ $$ |$$ |  $$ |  $$ |   $$ |  $$ |      $$ |      \$$$$$$  |$$ | \$$ |\$$$$$$  |  $$ |   $$$$$$\  $$$$$$  |$$ | \$$ |\$$$$$$  |
    ! \__|     \__|\__|  \__|  \__|   \__|  \__|      \__|       \______/ \__|  \__| \______/   \__|   \______| \______/ \__|  \__| \______/

    ! special functions

    pure function legendre(l,m,x) result(y)
        !
        ! Associated legrende polynomial
        !
        ! W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, Numerical Recipes in Fortran 77: The Art
        ! of Scientific Computing, 2 edition (Cambridge University Press, Cambridge England?; New York, 1992), p 246.
        !
        implicit none
        !
        integer , intent(in) :: l,m
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: y(:)
        integer  :: ll
        real(dp) :: pll(size(x))
        real(dp) :: pmm(size(x))
        real(dp) :: pmmp1(size(x))
        real(dp) :: somx2(size(x))
        !
        allocate(y(size(x)))
        !
        if ((m.lt.0).or.(m.gt.l).or.(all(abs(x).gt.1.0_dp))) then
            y = 0
        endif
        !
        pmm=1.0
        if (m > 0) then
            somx2=sqrt((1.0_dp-x)*(1.0_dp+x))
            pmm=product(arth(1.0_dp,2.0_dp,m))*somx2**m
            if (mod(m,2) == 1) pmm=-pmm
        end if
        if (l == m) then
            y=pmm
        else
            pmmp1=x*(2*m+1)*pmm
            if (l == m+1) then
                y=pmmp1
            else
                do ll=m+2,l
                    pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
                    pmm=pmmp1
                    pmmp1=pll
                end do
                y=pll
            end if
        end if
        contains
        pure function arth(first,increment,n)
            !
            implicit none
            !
            real(dp), intent(in) :: first
            real(dp), intent(in) :: increment
            integer , intent(in) :: n
            real(dp) :: arth(n)
            integer  :: k,k2
            real(dp) :: temp
            !
            if (n > 0) arth(1)=first
            if (n <= 16) then
                do k=2,n
                    arth(k)=arth(k-1)+increment
                enddo
            else
                do k=2,8
                    arth(k)=arth(k-1)+increment
                end do
                temp=increment*8
                k=8
                do
                    if (k >= n) exit
                    k2=k+k
                    arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
                    temp=temp+temp
                    k=k2
                enddo
            end if
        end function arth
    end function  legendre

    pure function laguerre(k,p,x) result(y)
        ! associated Laguerre polynomial L_k^p(x)(k,p,x)
        ! Using the expression from Samuel Shaw Ming Wong "Computational Methods in Physics and Engineering", Eq 4-86 p 139
        ! Note, there is a typographical error on http://mathworld.wolfram.com/AssociatedLaguerrePolynomial.html Eq 10
        ! Also see Linus Pauling "Introuction to Quantum Mechancs"
        implicit none 
        !
        integer , intent(in) :: k, p
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: y(:)
        integer :: j
        !
        allocate(y,mold=x)
        y=0.0_dp
        !
        do j = 0, k
            y = y + nchoosek(k+p,k-j) * (-x)**j / factorial(j)
        enddo
        !
    end function  laguerre

    pure function Ylm(l,m,theta,phi)
        !
        ! computes spherical harmonics. Theta and phi in radians. size(theta) must equal size(phi); together they define the radial coordinates.
        !
        ! W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, Numerical Recipes in Fortran 77: The Art
        ! of Scientific Computing, 2 edition (Cambridge University Press, Cambridge England?; New York, 1992), p 246.
        ! 
        implicit none
        !
        integer , intent(in) :: l, m
        real(dp), intent(in) :: theta(:), phi(:)
        complex(dp), allocatable :: Ylm(:)
        integer :: n
        ! 
        n = size(theta)
        !
        allocate(Ylm(n))
        !
        Ylm = sqrt( real(2*l+1,dp)/fourpi * factorial(l-m)/real(factorial(l+m),dp) ) * legendre(l,m,cos(theta)) * exp(cmplx_i*m*phi)
        !
    end function  Ylm

    pure function heavi(m)
        !
        ! A. V. Podolskiy and P. Vogl, Phys. Rev. B 69, 233101 (2004). Eq 15
        !
        implicit none
        !
        integer, intent(in) :: m
        real(dp) :: heavi
        !
        if (m.ge.0) then
            heavi = 1.0_dp
        else
            heavi = 0.0_dp
        endif
        !
    end function  heavi

    ! heaviside theta functions

    pure function fermi_dirac(x) result(y)
        !
        implicit none
        !
        real(dp), intent(in) :: x
        real(dp) :: y
        real(dp) :: maxarg
        maxarg = 200.0_dp
        ! Fermi-Dirac smearing
        if (x.lt.-maxarg) then
            y = 0.0_dp
        elseif (x.gt.maxarg) then
            y = 1.0_dp
        else
            y = 1.0_dp/(1.0_dp + exp(-x))
        endif
    end function  fermi_dirac
    
    pure function methfessel_paxton(x,n) result(y)
        ! theta function : PRB 40, 3616 (1989).
        implicit none
        !
        real(dp), intent(in) :: x
        integer , intent(in) :: n
        real(dp) :: y
        
        real(dp) :: a, hp, arg, hd
        integer :: i, ni
        real(dp), parameter :: maxarg = 200.0_dp
        
        ! Methfessel-Paxton
        y = gauss_freq(x * sqrt(2.0_dp) )
        if (n.eq.0) return
        hd = 0.0_dp
        arg = min(maxarg, x**2)
        hp = exp(- arg)
        ni = 0
        a = 1.0_dp / sqrt (pi)
        do i=1,n
            hd = 2.0_dp*x*hp-2.0_dp*real(ni,dp)*hd
            ni = ni+1
            a = -a/(real(i,dp)*4.0_dp)
            y = y-a*hd
            hp = 2.0_dp*x*hd-2.0_dp*real(ni,dp)*hp
            ni = ni+1
        enddo
        contains
        pure function gauss_freq(x)
            !     gauss_freq(x) = (1+erf(x/sqrt(2)))/2 = erfc(-x/sqrt(2))/2
            implicit none
            !
            real(dp),intent(in) :: x
            real(dp)            :: gauss_freq
            !
            gauss_freq = 0.5_dp * erfc(-x*0.7071067811865475_dp)
            !
        end function  gauss_freq
    end function  methfessel_paxton

    pure function marzari_vanderbilt(x) result(y)
        ! theta function: PRL 82, 3296 (1999)
        ! 1/2*erf(x-1/sqrt(2)) + 1/sqrt(2*pi)*exp(-(x-1/sqrt(2))**2) + 1/2
        implicit none
        !
        real(dp), intent(in) :: x
        real(dp) :: y
        real(dp) :: arg, xp
        real(dp) :: maxarg
        maxarg = 200.0_dp
        ! Cold smearing
         xp = x - 1.0_dp / sqrt (2.0_dp)
         arg = min (maxarg, xp**2)
         y = 0.5d0 * erf (xp) + 1.0_dp / sqrt (2.0_dp * pi) * exp (- arg) + 0.5d0
    end function  marzari_vanderbilt

    pure function erf(x)
      !---------------------------------------------------------------------
      !
      !     Error function - computed from the rational approximations of
      !     W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
      !
      !     for abs(x) le 0.47 erf is calculated directly
      !     for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
      !
      implicit none  
      real(dp), intent(in) :: x
      real(dp) :: x2, p1(4), q1(4)
      real(dp) :: erf
      !
      p1 = [2.426679552305318E2_dp, 2.197926161829415E1_dp, 6.996383488619136_dp,  -3.560984370181538E-2_dp]
      q1 = [2.150588758698612E2_dp, 9.116490540451490E1_dp, 1.508279763040779E1_dp, 1.000000000000000_dp]
      !
      if (abs(x).gt.6.0_dp) then  
         !  erf(6)=1-10^(-17) cannot be distinguished from 1
         erf = sign (1.0_dp, x)  
      else  
         if (abs(x).le.0.47_dp) then  
            x2 = x**2  
            erf = x*(p1(1)+x2*(p1(2)+x2*(p1(3)+x2*p1(4))))/(q1(1)+x2*(q1(2)+x2*(q1(3)+x2*q1(4))))
         else  
            erf = 1.0_dp - erfc(x)  
         endif
      endif
    end function  erf

    pure function erfc(x)
      !     erfc(x) = 1-erf(x)  - See comments in erf
      implicit none  
      real(dp),intent(in) :: x
      real(dp)            :: erfc
      real(dp) :: ax, xm2, p2(8), q2(8), p3(5), q3(5), pim1
      !
      p2 = [ 3.004592610201616E2_dp,  4.519189537118719E2_dp,  3.393208167343437E2_dp,  1.529892850469404E2_dp,  4.316222722205674E1_dp,  7.211758250883094_dp,    5.641955174789740E-1_dp,-1.368648573827167E-7_dp]
      q2 = [ 3.004592609569833E2_dp,  7.909509253278980E2_dp,  9.313540948506096E2_dp,  6.389802644656312E2_dp,  2.775854447439876E2_dp,  7.700015293522947E1_dp,  1.278272731962942E1_dp,  1.000000000000000_dp]
      p3 = [-2.996107077035422E-3_dp,-4.947309106232507E-2_dp,  -2.269565935396869E-1_dp,-2.786613086096478E-1_dp,  -2.231924597341847E-2_dp]
      q3 = [ 1.062092305284679E-2_dp, 1.913089261078298E-1_dp,  1.051675107067932_dp,    1.987332018171353_dp,     1.000000000000000_dp]
      pim1 = 0.56418958354775629_dp  ! sqrt(1/pi)
      ax = abs(x)  
      if (ax > 26.0_dp) then  
         !  erfc(26.0)=10^(-296); erfc(9.0)=10^(-37);
         erfc = 0.0_dp  
      elseif (ax.gt.4.0_dp) then
         xm2=(1.0_dp/ax)**2
         erfc=(1.0_dp/ax)*exp(-2.0_dp)*(pim1+xm2*(p3(1)+xm2*(p3(2)+xm2*(p3(3)+xm2*(p3(4)+xm2*p3(5)))))/(q3(1)+xm2*(q3(2)+xm2*(q3(3)+xm2*(q3(4)+xm2*q3(5))))))
      elseif(ax.gt.0.47_dp)then
         erfc=exp(-x**2)*(p2(1)+ax*(p2(2)+ax*(p2(3)+ax*(p2(4)+ax*(p2(5)+ax*(p2(6)+ax*(p2(7)+ax*p2(8))))))))/(q2(1)+ax*(q2(2)+ax*(q2(3)+ax*(q2(4)+ax*(q2(5)+ax*(q2(6)+ax*(q2(7)+ax*q2(8))))))))
      else
         erfc=1.0_dp-erf(ax)
      endif
      ! erf(-x)=-erf(x)  =>  erfc(-x) = 2-erfc(x)
      if (x < 0.0_dp) erfc = 2.0_dp - erfc 
    end function  erfc

    ! delta functions and theta function derivatives

    pure function fermi_dirac_dydx(x) result(y)
        ! derivative of Fermi-Dirac function: 0.5/(1.0+cosh(x))
        implicit none
        real(dp), intent(in) :: x
        real(dp) :: y
        !
        if (abs(x).le.36.0_dp) then
            y = 1.0_dp/(2.0_dp+exp(-x)+exp(+x))
        else
            y = 0.0_dp
        endif
    end function  fermi_dirac_dydx

    pure function marzari_vanderbilt_dydx(x) result(y)
        ! 1/sqrt(pi)*exp(-(x-1/sqrt(2))**2)*(2-sqrt(2)*x)
        implicit none
        real(dp), intent(in) :: x
        real(dp) :: y
        real(dp) :: arg
        real(dp) :: sqrtpm1
        !
        sqrtpm1 = 0.564189583547756_dp ! 1/sqrt(pi)
        arg = min(200.0_dp,(x-1.0_dp/sqrt(2.0_dp))**2)
        y = sqrtpm1*exp(-arg)*(2.0_dp-sqrt(2.0_dp)*x)
    end function  marzari_vanderbilt_dydx

    pure function methfessel_paxton_dydx(x,n) result(y)
        ! derivative of the corresponding Methfessel-Paxton wgauss
        implicit none
        real(dp), intent(in) :: x
        integer , intent(in) :: n
        real(dp) :: y
        real(dp) :: a, arg, hp, hd
        integer :: i, ni
        real(dp) :: sqrtpm1
        !
        sqrtpm1 = 0.564189583547756_dp ! 1/sqrt(pi)
        arg = min(200.0_dp, x**2)
        y = exp(-arg)*sqrtpm1
        if (n.eq.0) return
        hd = 0.0_dp
        hp = exp(-arg)
        ni = 0
        a  = sqrtpm1
        do i = 1, n
            hd = 2.0_dp*x*hp-2.0_dp*real(ni,dp)*hd
            ni = ni+1
            a  = -a/(real(i,dp)*4.0_dp)
            hp = 2.0_dp*x*hd-2.0_dp*real(ni,dp)*hp
            ni = ni+1
            y  = y+a*hp
        enddo
    end function  methfessel_paxton_dydx

    pure function gauss(x) result(y)
        ! derivative of Fermi-Dirac function: exp( - (x-x_o)^2 / (2*sigma^2) ) / sigma * sqrt( 2 * pi )
        ! degauss = (sqrt(2)*sigma) => exp(-(x-xo)**2/degauss**2)/(degauss*sqrt(pi))
        ! set degauss = 1 and center xo = 0
        implicit none
        real(dp), intent(in) :: x
        real(dp) :: y
        !
        if (abs(x).le.7.0_dp) then
            y = exp(-x**2.0_dp)/sqrtpi
            ! in order to avoid problems for large values of x in the e
        else
            y = 0.0_dp
        endif
    end function  gauss

    pure function lorentz(x) result(y)
        ! x = ( E(j,k)-Ep(i) )/degauss 
        ! derivative of Fermi-Dirac function: 1/( ((x-xo)/degauss)**2 + 1) * 1/(pi * degauss), 
        ! notice that the 1/degauss factor is taken outside of the equation 
        implicit none
        real(dp), intent(in) :: x
        real(dp) :: y
        !
        y = 1.0_dp/( pi * ( x**2 + 1 ) )
        !
    end function  lorentz

    ! $$$$$$\ $$\   $$\ $$$$$$$$\ $$$$$$$$\  $$$$$$\  $$$$$$$\   $$$$$$\ $$$$$$$$\ $$$$$$\  $$$$$$\  $$\   $$\
    !  \_$$  _|$$$\  $$ |\__$$  __|$$  _____|$$  __$$\ $$  __$$\ $$  __$$\\__$$  __|\_$$  _|$$  __$$\ $$$\  $$ |
    !    $$ |  $$$$\ $$ |   $$ |   $$ |      $$ /  \__|$$ |  $$ |$$ /  $$ |  $$ |     $$ |  $$ /  $$ |$$$$\ $$ |
    !    $$ |  $$ $$\$$ |   $$ |   $$$$$\    $$ |$$$$\ $$$$$$$  |$$$$$$$$ |  $$ |     $$ |  $$ |  $$ |$$ $$\$$ |
    !    $$ |  $$ \$$$$ |   $$ |   $$  __|   $$ |\_$$ |$$  __$$< $$  __$$ |  $$ |     $$ |  $$ |  $$ |$$ \$$$$ |
    !    $$ |  $$ |\$$$ |   $$ |   $$ |      $$ |  $$ |$$ |  $$ |$$ |  $$ |  $$ |     $$ |  $$ |  $$ |$$ |\$$$ |
    !  $$$$$$\ $$ | \$$ |   $$ |   $$$$$$$$\ \$$$$$$  |$$ |  $$ |$$ |  $$ |  $$ |   $$$$$$\  $$$$$$  |$$ | \$$ |
    !  \______|\__|  \__|   \__|   \________| \______/ \__|  \__|\__|  \__|  \__|   \______| \______/ \__|  \__|

    ! gaussian method

    function       get_dos_quick(Ep,E,kptw,degauss,flags) result(D)
        ! flags = fermi, mp, mv, gauss, lorentz
        implicit none
        !
        real(dp), intent(in) :: Ep(:)       ! probing energies
        real(dp), intent(in) :: E(:,:)      ! E(nbands,nkpts) band energies
        real(dp), intent(in) :: kptw(:)     ! kptw(nkpts) normalized kpoint weights
        real(dp), intent(in) :: degauss
        character(*), intent(in) :: flags
        real(dp), allocatable :: D(:)
        integer  :: nEs
        integer  :: nbands
        integer  :: nkpts
        integer  :: i,j,k
        real(dp) :: xp
        ! get number of probing energies
        nEs = size(Ep)
        ! get n of band
        nbands = size(E,1)
        ! get n of kpoints
        nkpts = size(kptw)
        ! allocate space for dos
        allocate(D(nEs))
        !$OMP PARALLEL PRIVATE(i,j,k) SHARED(D,E,Ep,nEs,kptw)
        !$OMP DO
        do i = 1, nEs
            ! initialize D(i)
            D(i) = 0.0_dp
            do j = 1, nbands
            do k = 1, nkpts
                xp = ( E(j,k)-Ep(i) )/degauss
                ! get contribution
                if     (index(flags,'fermi').ne.0) then
                    D(i) = D(i) + kptw(k) * fermi_dirac_dydx(x=xp)
                elseif (index(flags,'mp').ne.0) then
                    D(i) = D(i) + kptw(k) * methfessel_paxton_dydx(x=xp,n=1)
                elseif (index(flags,'mv').ne.0) then
                    D(i) = D(i) + kptw(k) * marzari_vanderbilt_dydx(x=xp)
                elseif (index(flags,'gauss').ne.0) then
                    D(i) = D(i) + kptw(k) * gauss(x=xp)
                elseif (index(flags,'lorentz').ne.0) then
                    D(i) = D(i) + kptw(k) * lorentz(x=xp)
                else
                    stop 'ERROR [get_dos_quick]: flag not recognized'
                endif
            enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        D = D / degauss
    end function   get_dos_quick

    function       get_pdos_quick(Ep,E,kptw,degauss,W,flags) result(pD)
        ! flags = fermi, mp, mv, gauss, lorentz
        implicit none
        !
        real(dp), intent(in) :: Ep(:)    ! probing energies
        real(dp), intent(in) :: E(:,:)   ! E(nbands,nkpts) band energies
        real(dp), intent(in) :: W(:,:,:) ! W(nprojections,nbands,nkpts) band weights
        real(dp), intent(in) :: kptw(:)  ! kptw(nkpts) normalized kpoint weights
        real(dp), intent(in) :: degauss
        character(*), intent(in) :: flags
        real(dp), allocatable :: pD(:,:)
        integer  :: nEs
        integer  :: nbands
        integer  :: nkpts
        integer  :: nprojections
        integer  :: i,j,k,l
        real(dp) :: xp,yp
        ! get number of probing energies
        nEs = size(Ep)
        ! get n of band
        nbands = size(E,1)
        ! get n of kpoints
        nkpts = size(kptw)
        ! numbe of projections
        nprojections = size(W,1)
        ! allocate space for dos
        allocate(pD(nprojections,nEs))
        ! initialize
        pD = 0.0_dp
        !$OMP PARALLEL PRIVATE(i,j,k,l) SHARED(pD,E,Ep,nEs,kptw)
        !$OMP DO
        do i = 1, nEs
            do j = 1, nbands
            do k = 1, nkpts
                xp = ( E(j,k)-Ep(i) )/degauss
                ! get contribution from this band
                if     (index(flags,'fermi').ne.0) then
                    yp = fermi_dirac_dydx(x=xp)
                elseif (index(flags,'mp').ne.0) then
                    yp = methfessel_paxton_dydx(x=xp,n=1)
                elseif (index(flags,'mv').ne.0) then
                    yp = marzari_vanderbilt_dydx(x=xp)
                elseif (index(flags,'gauss').ne.0) then
                    yp = gauss(x=xp)
                elseif (index(flags,'lorentz').ne.0) then
                    yp = lorentz(x=xp)
                else
                    stop 'ERROR [get_dos_quick]: flag not recognized'
                endif
                ! multiply by kpoint weight
                yp = yp * kptw(k)
                ! multiply by band weights
                do l = 1, nprojections
                    pD(l,i) = pD(l,i) + yp * W(l,j,k)
                enddo
            enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        pD = pD / degauss
    end function   get_pdos_quick

    ! DOS properties

    function       get_Ef(E,tet,tetw,nelecs) result(Ef)
        ! adapted from qe
        implicit none
        !
        real(dp), intent(in) :: E(:,:)
        integer , intent(in) :: tet(:,:)
        real(dp), intent(in) :: tetw(:)
        real(dp), intent(in) :: nelecs
        real(dp) :: Ef
        integer  :: i
        integer  :: ntets
        integer  :: nbands
        integer  :: nkpts
        real(dp) :: Emin    ! absolute lowest band energy
        real(dp) :: Emax    ! absolute highest band energy
        real(dp) :: elw     ! lower bound for EF
        real(dp) :: eup     ! upper bound for EF
        real(dp) :: idos_up ! number of states for Ef=eup
        real(dp) :: idos_lw ! number of states for Ef=elw
        real(dp) :: idos_ep ! number of states for Ef=(elw+eup)/2
        real(dp) :: try_Ef  ! new Ef to try
        real(dp) :: crit    ! stopping criterion
        ! get bands, kpts, tetrahedra
        nbands = size(E,1)
        nkpts  = size(E,2)
        ntets  = sizE(tet,2)
        ! get lower and upper bracketing bounds
        Emin = minval(E)
        Emax = maxval(E)
        elw = Emin
        eup = Emax
        ! initailize bracketing 
        idos_up = get_tet_idos_engine(Ep=eup,E=E,tet=tet,tetw=tetw)
        idos_lw = get_tet_idos_engine(Ep=elw,E=E,tet=tet,tetw=tetw)
        crit = 1.0d+10
        ! search for Ef by gradualling reducing bracketted region
        search : do i = 1, 1000 ! hard-coded 1000 maximum iterations
            try_Ef = (eup+elw) / 2.0_dp
            idos_ep = get_tet_idos_engine(Ep=try_Ef,E=E,tet=tet,tetw=tetw)
            if (abs(idos_ep-nelecs).lt.crit) then
                crit = abs(idos_ep-nelecs)
                Ef = try_Ef
            endif
            ! converged
            if (abs(idos_ep-nelecs).lt.tiny) then
                exit search
            elseif((idos_ep-nelecs).lt.-tiny) then
                elw = try_Ef
            else
                eup = try_Ef
            endif
        enddo search
        ! check convergence
        if (abs(idos_ep-nelecs).gt.tiny) stop 'ERROR [XX] failed to converge on Ef.'
        ! check that Ef < highest band
        if (Ef.gt.Emax) stop 'ERROR [XX]: Ef is above the highest band energy'
        !
    end function   get_Ef

    ! tetrahedra methods

    function       get_dos_tet(Ep,E,tet,tetw) result(D)
        !
        implicit none
        !
        real(dp), intent(in) :: Ep(:) ! probing energies
        real(dp), intent(in) :: E(:,:) ! E(nbands,nkpts) band energies
        integer , intent(in) :: tet(:,:) ! tet(:,ntets) tetrahedra conenctivity
        real(dp), intent(in) :: tetw(:) ! tetw(ntets) weight of each tetrahedron
        real(dp), allocatable :: D(:)
        integer :: nEs
        integer :: i
        ! get number of probing energies
        nEs = size(Ep)
        ! allocate space for dos
        allocate(D(nEs))
        !$OMP PARALLEL PRIVATE(i) SHARED(D,nEs,E,Ep,tet,tetw)
        !$OMP DO
        do i = 1, nEs
            ! get dos
            D(i) = get_tet_dos_engine(Ep=Ep(i),E=E,tet=tet,tetw=tetw)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    end function   get_dos_tet

    function       get_idos_tet(Ep,E,tet,tetw) result(iD)
        !
        implicit none
        !
        real(dp), intent(in) :: Ep(:) ! probing energies
        real(dp), intent(in) :: E(:,:) ! E(nbands,nkpts) band energies
        integer , intent(in) :: tet(:,:) ! tet(:,ntets) tetrahedra conenctivity
        real(dp), intent(in) :: tetw(:) ! tetw(ntets) weight of each tetrahedron
        real(dp), allocatable :: iD(:)
        integer :: nEs
        integer :: i
        ! get number of probing energies
        nEs = size(Ep)
        ! allocate space for dos
        allocate(iD(nEs))
        !$OMP PARALLEL PRIVATE(i) SHARED(iD,nEs,E,Ep,tet,tetw)
        !$OMP DO
        do i = 1, nEs
            ! get dos
            iD(i) = get_tet_idos_engine(Ep=Ep(i),E=E,tet=tet,tetw=tetw)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    end function   get_idos_tet

    function       get_pdos_tet(Ep,E,tet,tetw,weight) result(pD)
        !
        implicit none
        !
        real(dp), intent(in) :: Ep(:) ! probing energies
        real(dp), intent(in) :: E(:,:) ! E(nbands,nkpts) band energies
        integer , intent(in) :: tet(:,:) ! tet(:,ntets) tetrahedra conenctivity
        real(dp), intent(in) :: tetw(:) ! tetw(ntets) weight of each tetrahedron
        real(dp), intent(in) :: weight(:,:,:) ! weights(nprojections,nbands,nkpts), weights (per band per kpoint) for projected dos: can be absolute square of TB or BvK eigenvectors
        real(dp), allocatable :: pD(:,:)
        integer  :: nEs
        integer  :: nbands
        integer  :: ntets
        integer  :: nprojections
        real(dp) :: wc(4) ! corner tetrahedron weights
        integer  :: i,j,k,m
        ! get number of probing energies
        nEs = size(Ep)
        ! get n of band
        nbands = size(E,1)
        ! get number of tetrahedra
        ntets = size(tet,2)
        ! get number of projections
        nprojections = size(weight,1)
        ! allocate space 
        allocate(pD(nprojections,nEs))
        ! loop over energies, bands, tetrahedra
        !$OMP PARALLEL PRIVATE(i,j,k,m,wc) SHARED(pD,E,Ep,nEs,tetw,weight)
        !$OMP DO
        do i = 1, nEs
            ! initialize D(i)
            pD(:,i) = 0.0_dp
            do j = 1, nbands
            do k = 1, ntets
                ! get tetrahedron corner weights
                wc = get_delta_wc(Ep=Ep(i),Ec=E(j,tet(:,k)))
                ! loop over projections, increment pDOS with contributions from band j in tetrahedron k
                do m = 1, nprojections
                    pD(m,i) = pD(m,i) + tetw(k) * sum(wc * weight(m,j,tet(:,k)) )
                enddo
            enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    end function   get_pdos_tet

    pure function  get_theta_wc(Ep,Ec,ntets) result(wc)
        ! get linear tetrahedron corner weights for delta function with Blochl corrections
        implicit none
        !
        real(dp), intent(in) :: Ep    ! probing energy
        real(dp), intent(in) :: Ec(4) ! corner energies
        integer , intent(in) :: ntets
        real(dp) :: wc(4)
        real(dp) :: dosEp
        real(dp) :: etot
        real(dp) :: c1,c2,c3,c4
        real(dp) :: e1,e2,e3,e4
        integer  :: k1,k2,k3,k4
        integer  :: inds(4)
        ! rank corner weights in increasing order
        inds = fourrank(Ec)
        ! sort weights
        e1 = Ec(inds(1))
        e2 = Ec(inds(2))
        e3 = Ec(inds(3))
        e4 = Ec(inds(4))
        ! k1-k4 are the irreducible k-points corresponding to e1-e4
        k1 = inds(1)
        k2 = inds(2)
        k3 = inds(3)
        k4 = inds(4)
        ! calculate weights wc
        if     (Ep.le.e1) then
            ! Eq B1 from Blochl's PhysRevB.49.16223
            wc(k1)=0.0_dp
            wc(k2)=0.0_dp
            wc(k3)=0.0_dp
            wc(k4)=0.0_dp
        elseif (Ep.le.e2) then
            ! Eq B6 from Blochl's PhysRevB.49.16223
            c4=0.25_dp/ntets*(Ep-e1)**3/(e2-e1)/(e3-e1)/(e4-e1)
            ! Eq C2 from Blochl's PhysRevB.49.16223
            dosEp=3.0_dp/ntets*(Ep-e1)**2/(e2-e1)/(e3-e1)/(e4-e1)
            ! shortcut for Blochl corrections (Eq.22 PhysRevB.49.16223)
            etot=e1+e2+e3+e4
            ! Eq B2 from Blochl's PhysRevB.49.16223
            wc(k1)=c4*(4.0_dp-(Ep-e1)*(1.0_dp/(e2-e1)+1.0_dp/(e3-e1)+1.0_dp/(e4-e1)))+dosEp*(etot-4.0_dp*e1)*0.025_dp
            ! Eq B3 from Blochl's PhysRevB.49.16223
            wc(k2)=c4*(Ep-e1)/(e2-e1)+dosEp*(etot-4.0_dp*e2)*0.025_dp
            ! Eq B4 from Blochl's PhysRevB.49.16223
            wc(k3)=c4*(Ep-e1)/(e3-e1)+dosEp*(etot-4.0_dp*e3)*0.025_dp
            ! Eq B5 from Blochl's PhysRevB.49.16223
            wc(k4)=c4*(Ep-e1)/(e4-e1)+dosEp*(etot-4.0_dp*e4)*0.025_dp
        elseif (Ep.le.e3) then
            ! Eq B11 from Blochl's PhysRevB.49.16223
            c1=0.25_dp/ntets*(Ep-e1)**2/(e4-e1)/(e3-e1)
            ! Eq B12 from Blochl's PhysRevB.49.16223
            c2=0.25_dp/ntets*(Ep-e1)*(Ep-e2)*(e3-Ep)/(e4-e1)/(e3-e2)/(e3-e1)
            ! Eq B13 from Blochl's PhysRevB.49.16223
            c3=0.25_dp/ntets*(Ep-e2)**2*(e4-Ep)/(e4-e2)/(e3-e2)/(e4-e1)
            ! Eq C3 from Blochl's PhysRevB.49.16223
            dosEp=1.0_dp/ntets/(e3-e1)/(e4-e1)*(3.0_dp*(e2-e1)+6.0_dp*(Ep-e2)-3.0_dp*(e3-e1+e4-e2)*(Ep-e2)**2/(e3-e2)/(e4-e2))
            ! shortcut for Blochl corrections (Eq.22 PhysRevB.49.16223)
            etot=e1+e2+e3+e4
            ! Eq B7 from Blochl's PhysRevB.49.16223
            wc(k1)=c1+(c1+c2)*(e3-Ep)/(e3-e1)+(c1+c2+c3)*(e4-Ep)/(e4-e1)+dosEp*(etot-4.0_dp*e1)*0.025_dp
            ! Eq B8 from Blochl's PhysRevB.49.16223
            wc(k2)=c1+c2+c3+(c2+c3)*(e3-Ep)/(e3-e2)+c3*(e4-Ep)/(e4-e2)+dosEp*(etot-4.0_dp*e2)*0.025_dp
            ! Eq B9 from Blochl's PhysRevB.49.16223
            wc(k3)=(c1+c2)*(Ep-e1)/(e3-e1)+(c2+c3)*(Ep-e2)/(e3-e2)+dosEp*(etot-4.0_dp*e3)*0.025_dp
            ! Eq B10 from Blochl's PhysRevB.49.16223
            wc(k4)=(c1+c2+c3)*(Ep-e1)/(e4-e1)+c3*(Ep-e2)/(e4-e2)+dosEp*(etot-4.0_dp*e4)*0.025_dp
        elseif (Ep.le.e4) then
            ! Eq B18 from Blochl's PhysRevB.49.16223
            c4=0.25_dp/ntets*(e4-Ep)**3/(e4-e1)/(e4-e2)/(e4-e3)
            ! Eq C4 from Blochl's PhysRevB.49.16223
            dosEp=3.0_dp/ntets*(e4-Ep)**2/(e4-e1)/(e4-e2)/(e4-e3)
            ! shortcut for Blochl corrections (Eq.22 PhysRevB.49.16223)
            etot=e1+e2+e3+e4
            ! Eq B14 from Blochl's PhysRevB.49.16223
            wc(k1)=0.25_dp/ntets-c4*(e4-Ep)/(e4-e1)+dosEp*(etot-4.0_dp*e1)*0.025_dp
            ! Eq B15 from Blochl's PhysRevB.49.16223
            wc(k2)=0.25_dp/ntets-c4*(e4-Ep)/(e4-e2)+dosEp*(etot-4.0_dp*e2)*0.025_dp
            ! Eq B16 from Blochl's PhysRevB.49.16223
            wc(k3)=0.25_dp/ntets-c4*(e4-Ep)/(e4-e3)+dosEp*(etot-4.0_dp*e3)*0.025_dp
            ! Eq B17 from Blochl's PhysRevB.49.16223
            wc(k4)=0.25_dp/ntets-c4*(4.0_dp-(e4-Ep)*(1.0_dp/(e4-e1)+1.0_dp/(e4-e2)+1.0_dp/(e4-e3)))+dosEp*(etot-4.0_dp*e4)*0.025_dp
        elseif (Ep.ge.e4) then
            ! Eq B19 from Blochl's PhysRevB.49.16223
            wc(k1)=0.25_dp/ntets
            wc(k2)=0.25_dp/ntets
            wc(k3)=0.25_dp/ntets
            wc(k4)=0.25_dp/ntets
        endif
        ! 
        wc = wc * ntets
        !
    end function   get_theta_wc

    pure function  get_delta_wc(Ep,Ec) result(wc)
        ! get linear tetrahedron corner weights for delta function with Blochl corrections
        implicit none
        !
        real(dp), intent(in) :: Ep    ! probing energy
        real(dp), intent(in) :: Ec(4) ! corner energies
        real(dp) :: wc(4)
        real(dp) :: dosEp
        real(dp) :: e1,e2,e3,e4
        integer  :: k1,k2,k3,k4
        real(dp) :: o13,f12,f13,f14,f21,f23,f31,f32,f34,f24,f41,f42
        integer  :: inds(4)
        ! initialize shortcut
        o13  = 1.0_dp/3.0_dp
        ! rank corner weights in increasing order
        inds = fourrank(Ec)
        ! sort weights
        e1 = Ec(inds(1))
        e2 = Ec(inds(2))
        e3 = Ec(inds(3))
        e4 = Ec(inds(4))
        ! k1-k4 are the irreducible k-points corresponding to e1-e4
        k1 = inds(1)
        k2 = inds(2)
        k3 = inds(3)
        k4 = inds(4)
        ! calculate weights wc
        if     (Ep.le.e1) then
            ! Eq B1 from Blochl's PhysRevB.49.16223
            wc(k1)=0.0_dp
            wc(k2)=0.0_dp
            wc(k3)=0.0_dp
            wc(k4)=0.0_dp
        elseif (Ep.lt.e2) then
            f12 = (Ep-e2)/(e1-e2)
            f21 = 1.0_dp - f12
            f13 = (Ep-e3)/(e1-e3)
            f31 = 1.0_dp - f13
            f14 = (Ep-e4)/(e1-e4)
            f41 = 1.0_dp - f14
            dosEp  = 3.0_dp * f21 * f31 * f41 / (Ep-e1)
            wc(k1) = o13 * (f12 + f13 + f14)
            wc(k2) = o13 * f21
            wc(k3) = o13 * f31
            wc(k4) = o13 * f41
            wc = wc * dosEp
        elseif (Ep.lt.e3) then
            f13 = (Ep-e3)/(e1-e3)
            f31 = 1.0_dp - f13
            f14 = (Ep-e4)/(e1-e4)
            f41 = 1.0_dp-f14
            f23 = (Ep-e3)/(e2-e3)
            f32 = 1.0_dp - f23
            f24 = (Ep-e4)/(e2-e4)
            f42 = 1.0_dp - f24
            dosEp  = 3.0_dp * (f23*f31 + f32*f24)
            wc(k1) = f14 * o13 + f13*f31*f23 / dosEp
            wc(k2) = f23 * o13 + f24*f24*f32 / dosEp
            wc(k3) = f32 * o13 + f31*f31*f23 / dosEp
            wc(k4) = f41 * o13 + f42*f24*f32 / dosEp
            dosEp  = dosEp / (e4-e1)
            wc = wc * dosEp
        elseif (Ep.lt.e4) then
            f14 = (Ep-e4)/(e1-e4)
            f24 = (Ep-e4)/(e2-e4)
            f34 = (Ep-e4)/(e3-e4)
            dosEp  = 3.0_dp * f14 * f24 * f34 / (e4-Ep)
            wc(k1) = f14 * o13
            wc(k2) = f24 * o13
            wc(k3) = f34 * o13
            wc(k4) = (3.0_dp - f14 - f24 - f34 ) * o13
            wc = wc * dosEp
        elseif (Ep.gt.e4) then
            wc(k1)=0.0_dp
            wc(k2)=0.0_dp
            wc(k3)=0.0_dp
            wc(k4)=0.0_dp
       endif
       !
    end function   get_delta_wc

    pure function  get_tet_idos_engine(Ep,E,tet,tetw) result(idos)
        ! return integrated DOS at Ep
        implicit none
        !
        real(dp), intent(in) :: Ep
        real(dp), intent(in) :: E(:,:)
        integer , intent(in) :: tet(:,:)
        real(dp), intent(in) :: tetw(:) ! weight of irreducible tetrahedron
        real(dp) :: idos
        real(dp) :: Ec(4)
        real(dp) :: e1,e2,e3,e4
        integer  :: i,j
        integer  :: nbands
        integer  :: ntets
        ! get bands, kpts, tetrahedra
        nbands = size(E,1)
        ntets  = sizE(tet,2)
        !
        idos = 0.0_dp
        do j = 1, ntets
            do i = 1, nbands
                ! get energies at the vertices of the j-th tethedron in ascending order
                Ec(1:4) = foursort(E(i,tet(:,j)))
                ! map to 
                e1 = Ec(1)
                e2 = Ec(2)
                e3 = Ec(3)
                e4 = Ec(4)
                ! calculate sum over k of the integrated charge
                if     (Ep.le.e1) then
                    ! Eq A1 from Blochl's PhysRevB.49.16223
                    ! do nothing
                elseif (Ep.le.e2) then
                    ! Eq A2 from Blochl's PhysRevB.49.16223
                    idos = idos + tetw(j)*(Ep-e1)**3/(e2-e1)/(e3-e1)/(e4-e1)
                elseif (Ep.le.e3) then
                    ! Eq A3 from Blochl's PhysRevB.49.16223
                    idos = idos + tetw(j)/(e3-e1)/(e4-e1)*((e2-e1)**2+3.0_dp*(e2-e1)*(Ep-e2)+3.0_dp*(Ep-e2)**2.0_dp-(e3-e1+e4-e2)/(e3-e2)/(e4-e2)*(Ep-e2)**3.0_dp)
                elseif (Ep.le.e4) then
                    ! Eq A4 from Blochl's PhysRevB.49.16223
                    idos = idos + tetw(j)*(1.0_dp-(e4-Ep)**3.0_dp/(e4-e1)/(e4-e2)/(e4-e3))
                elseif (Ep.ge.e4) then
                    ! Eq A5 from Blochl's PhysRevB.49.16223
                    idos = idos + tetw(j) 
                endif
            enddo
        enddo
    end function   get_tet_idos_engine

    pure function  get_tet_dos_engine(Ep,E,tet,tetw) result(dosEp)
        ! returns DOS at Ep
        implicit none
        !
        real(dp), intent(in) :: Ep
        real(dp), intent(in) :: E(:,:)
        integer , intent(in) :: tet(:,:)
        real(dp), intent(in) :: tetw(:) ! weight of irreducible tetrahedron
        real(dp) :: dosEp
        real(dp) :: Ec(4)
        real(dp) :: e1,e2,e3,e4
        integer :: i,j
        integer :: nbands
        integer :: ntets
        ! get bands, kpts, tetrahedra
        nbands = size(E,1)
        ntets  = sizE(tet,2)
        !
        dosEp = 0.0_dp
        do j = 1, ntets
            do i = 1, nbands
                ! get energies at the verticies of the j-th tethedron in ascending order
                Ec(1:4) = foursort(E(i,tet(:,j)))
                ! map to 
                e1 = Ec(1)
                e2 = Ec(2)
                e3 = Ec(3)
                e4 = Ec(4)
                ! calculate sum over k of the integrated charge
                if     (Ep.le.e1) then
                    ! Eq C1 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp*1.0_dp
                    ! do nothing
                elseif (Ep.le.e2) then
                    ! Eq C2 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp + tetw(j)*3.0_dp*(Ep-e1)**2.0_dp/(e2-e1)/(e3-e1)/(e4-e1)
                elseif (Ep.le.e3) then
                    ! Eq C3 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp + tetw(j)/(e3-e1)/(e4-e1)*(3.0_dp*(e2-e1)+6.0_dp*(Ep-e2)-3.0_dp*(e3-e1+e4-e2)/(e3-e2)/(e4-e2)*(Ep-e2)**2.0_dp)
                elseif (Ep.le.e4) then
                    ! Eq C4 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp + tetw(j)*(3.0_dp*(e4-Ep)**2.0_dp/(e4-e1)/(e4-e2)/(e4-e3))
                elseif (Ep.ge.e4) then
                    ! Eq C1 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp*1.0_dp
                    ! do nothing
                endif
            enddo
        enddo
        !
    end function   get_tet_dos_engine

    ! $$\      $$\ $$$$$$\  $$$$$$\   $$$$$$\  $$$$$$$$\ $$\       $$\        $$$$$$\  $$\   $$\ $$$$$$$$\  $$$$$$\  $$\   $$\  $$$$$$\
    ! $$$\    $$$ |\_$$  _|$$  __$$\ $$  __$$\ $$  _____|$$ |      $$ |      $$  __$$\ $$$\  $$ |$$  _____|$$  __$$\ $$ |  $$ |$$  __$$\
    ! $$$$\  $$$$ |  $$ |  $$ /  \__|$$ /  \__|$$ |      $$ |      $$ |      $$ /  $$ |$$$$\ $$ |$$ |      $$ /  $$ |$$ |  $$ |$$ /  \__|
    ! $$\$$\$$ $$ |  $$ |  \$$$$$$\  $$ |      $$$$$\    $$ |      $$ |      $$$$$$$$ |$$ $$\$$ |$$$$$\    $$ |  $$ |$$ |  $$ |\$$$$$$\
    ! $$ \$$$  $$ |  $$ |   \____$$\ $$ |      $$  __|   $$ |      $$ |      $$  __$$ |$$ \$$$$ |$$  __|   $$ |  $$ |$$ |  $$ | \____$$\
    ! $$ |\$  /$$ |  $$ |  $$\   $$ |$$ |  $$\ $$ |      $$ |      $$ |      $$ |  $$ |$$ |\$$$ |$$ |      $$ |  $$ |$$ |  $$ |$$\   $$ |
    ! $$ | \_/ $$ |$$$$$$\ \$$$$$$  |\$$$$$$  |$$$$$$$$\ $$$$$$$$\ $$$$$$$$\ $$ |  $$ |$$ | \$$ |$$$$$$$$\  $$$$$$  |\$$$$$$  |\$$$$$$  |
    ! \__|     \__|\______| \______/  \______/ \________|\________|\________|\__|  \__|\__|  \__|\________| \______/  \______/  \______/

    pure subroutine uc2ws(K,M,tiny)

        implicit none

        integer :: j
        real(dp), allocatable, intent(inout) :: K(:,:)
        real(dp), intent(in) :: M(3,3)
        real(dp), intent(in) :: tiny
        real(dp) :: G(3,26) , G2(26)

        ! generate mesh
        G(1,:) = real([-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1])
        G(2,:) = real([-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1])
        G(3,:) = real([-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        G = matmul(M,G)
        G2 = sum(G**2,1)

        ! perform calc
        ! do j = 1, size(K,2)
        do concurrent (j = 1:size(K,2))
            call uc2ws_engine(K(:,j),G,G2,tiny)
        enddo
    end subroutine  uc2ws

    pure subroutine uc2ws_engine(K,G,G2,tiny)

        implicit none

        real(dp), intent(inout) :: K(3)
        real(dp), intent(in) :: G(3,26), G2(26), tiny
        real(dp)  :: P
        logical :: go
        integer :: i

        go = .true.
        do while (go)
            go = .false.
            do i = 1, 26
                P = 2*(K(1)*G(1,i) + K(2)*G(2,i) + K(3)*G(3,i)) - G2(i)
                if (P .gt. 1D-8) then
                    K = K - G(:,i)
                    go = .true.
                endif
            enddo
        enddo
    end subroutine  uc2ws_engine

end module
%}
            
            % mex interfaces
            i=i+1; f{i} = 'uc2ws_mex'; fid=fopen([f{i},'.f90'],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); 
%{
#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    use am_mex_lib , only : uc2ws

    implicit none

    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs
    mwPointer mxGetPr
    mwPointer mxCreateDoubleMatrix
    mwPointer mxGetM, mxGetN

    real*8, allocatable :: K(:,:)
    real*8 :: tiny, G(3,26) , G2(26) , M(3,3)

    ! [K] = uc2ws(K,M,tiny)
    if(nrhs .ne. 3) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','Three inputs required.')
    if(nlhs .gt. 1) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','One output produced.')
    if(mxGetM(prhs(1)) .ne. 3) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','Input 1 requires three-dimensional column vectors.')
    if(mxGetM(prhs(2)) .ne. 3) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','Input 2 requires a three-dimensional transformation matrix.')
    if(mxGetN(prhs(2)) .ne. 3) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','Input 2 requires a three-dimensional transformation matrix.')
    if(mxGetM(prhs(3)) .ne. 1 .or. &
     & mxGetN(prhs(3)) .ne. 1) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','Input 3 requires a scalar.')

    ! inputs 
    ! K
    allocate( K(3,mxGetN(prhs(1))) )
    call mxCopyPtrToReal8(mxGetPr(prhs(1)),K,mxGetM(prhs(1))*mxGetN(prhs(1)))
    ! M(3,3)
    call mxCopyPtrToReal8(mxGetPr(prhs(2)),M,mxGetM(prhs(2))*mxGetN(prhs(2)))
    ! tiny
    call mxCopyPtrToReal8(mxGetPr(prhs(3)),tiny,mxGetM(prhs(3)))

    ! execute code
    call uc2ws(K,M,tiny)

    ! output results
    plhs(1) = mxCreateDoubleMatrix(mxGetM(prhs(1)),mxGetN(prhs(1)),0)
    call mxCopyReal8ToPtr(K,mxGetPr(plhs(1)),mxGetM(prhs(1))*mxGetN(prhs(1)))
end subroutine mexFunction
%}
            i=i+1; f{i} = 'get_dos_tet_mex'; fid=fopen([f{i},'.f90'],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); 
%{
#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    use am_mex_lib , only : get_dos_tet

    implicit none

    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs
    mwPointer mxGetPr
    mwPointer mxCreateDoubleMatrix
    mwPointer mxGetM, mxGetN

    real*8, allocatable :: Ep(:)
    real*8, allocatable :: E(:,:)
    real*8, allocatable :: tetw(:)
    real*8, allocatable :: tet(:,:)
    real*8, allocatable :: D(:)

    integer :: i

    ! [D] = get_dos_tet(Ep,E,tet,tetw)
    if(nrhs .ne. 4) call mexErrMsgIdAndTxt ('MATLAB:get_dos_tet_mex','Four inputs required.')
    if(nlhs .gt. 1) call mexErrMsgIdAndTxt ('MATLAB:get_dos_tet_mex','One output produced.')
    if(mxGetM(prhs(3)) .ne. 4) call mexErrMsgIdAndTxt ('MATLAB:get_dos_tet_mex','Connectivity list requires four numbers in each column.')
    if(mxGetN(prhs(3)) .ne. mxGetN(prhs(4))) call mexErrMsgIdAndTxt ('MATLAB:get_dos_tet_mex','Mismatched number of tetrahedra in weights and connectivity list.')

    ! inputs 
    i=0
    ! 1) Ep
    i=i+1
    allocate(   Ep( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),Ep, mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 2) E(nbands,nkpts)
    i=i+1
    allocate(   E( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),E,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 3) tet(1:4,ntets)
    i=i+1
    allocate( tet( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tet,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 4) tetw(ntets)
    i=i+1
    allocate( tetw( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tetw,mxGetM(prhs(i))*mxGetN(prhs(i)))

    ! execute code
    D = get_dos_tet(Ep,E,nint(tet),tetw)

    ! outputs
    i=0
    ! output results
    i=i+1
    plhs(i) = mxCreateDoubleMatrix(mxGetM(prhs(1)),mxGetN(prhs(1)),0)
    call mxCopyReal8ToPtr(D,mxGetPr(plhs(1)),mxGetM(prhs(1))*mxGetN(prhs(1)))
end subroutine mexFunction
%}
            i=i+1; f{i} = 'get_idos_tet_mex'; fid=fopen([f{i},'.f90'],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); 
%{
#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    use am_mex_lib , only : get_idos_tet

    implicit none

    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs
    mwPointer mxGetPr
    mwPointer mxCreateDoubleMatrix
    mwPointer mxGetM, mxGetN

    real*8, allocatable :: Ep(:)
    real*8, allocatable :: E(:,:)
    real*8, allocatable :: tetw(:)
    real*8, allocatable :: tet(:,:)
    real*8, allocatable :: D(:)

    integer :: i

    ! [D] = get_dos_tet(Ep,E,tet,tetw)
    if(nrhs .ne. 4) call mexErrMsgIdAndTxt ('MATLAB:get_idos_tet_mex','Four inputs required.')
    if(nlhs .gt. 1) call mexErrMsgIdAndTxt ('MATLAB:get_idos_tet_mex','One output produced.')
    if(mxGetM(prhs(3)) .ne. 4) call mexErrMsgIdAndTxt ('MATLAB:get_idos_tet_mex','Connectivity list requires four numbers in each column.')
    if(mxGetN(prhs(3)) .ne. mxGetN(prhs(4))) call mexErrMsgIdAndTxt ('MATLAB:get_idos_tet_mex','Mismatched number of tetrahedra in weights and connectivity list.')

    ! inputs 
    i=0
    ! 1) Ep
    i=i+1
    allocate(   Ep( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),Ep, mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 2) E(nbands,nkpts)
    i=i+1
    allocate(   E( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),E,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 3) tet(1:4,ntets)
    i=i+1
    allocate( tet( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tet,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 4) tetw(ntets)
    i=i+1
    allocate( tetw( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tetw,mxGetM(prhs(i))*mxGetN(prhs(i)))

    ! execute code
    D = get_idos_tet(Ep,E,nint(tet),tetw)

    ! outputs
    i=0
    ! output results
    i=i+1
    plhs(i) = mxCreateDoubleMatrix(mxGetM(prhs(1)),mxGetN(prhs(1)),0)
    call mxCopyReal8ToPtr(D,mxGetPr(plhs(1)),mxGetM(prhs(1))*mxGetN(prhs(1)))
end subroutine mexFunction
%}
            i=i+1; f{i} = 'get_pdos_tet_mex'; fid=fopen([f{i},'.f90'],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); 
%{
#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    use am_mex_lib , only : get_pdos_tet

    implicit none

    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs
    mwPointer mxGetPr
    mwPointer mxCreateDoubleMatrix
    mwPointer mxGetM, mxGetN

    real*8, allocatable :: Ep(:)
    real*8, allocatable :: E(:,:)
    real*8, allocatable :: tetw(:)
    real*8, allocatable :: tet(:,:)
    real*8, allocatable :: weight(:,:,:)
    real*8, allocatable :: D(:,:)

    integer :: i

    ! [D] = get_dos_tet(Ep,E,tet,tetw)
    if(nrhs .ne. 5) call mexErrMsgIdAndTxt ('MATLAB:get_pdos_tet_mex','Five inputs required.')
    if(nlhs .gt. 1) call mexErrMsgIdAndTxt ('MATLAB:get_pdos_tet_mex','One output produced.')
    if(mxGetM(prhs(3)) .ne. 4) call mexErrMsgIdAndTxt ('MATLAB:get_pdos_tet_mex','Connectivity list requires four numbers in each column.')
    if(mxGetN(prhs(3)) .ne. mxGetN(prhs(4))) call mexErrMsgIdAndTxt ('MATLAB:get_pdos_tet_mex','Mismatched number of tetrahedra in weights and connectivity list.')
    if(mxGetN(prhs(2)) .ne. mxGetN(prhs(5))) call mexErrMsgIdAndTxt ('MATLAB:get_pdos_tet_mex','Mismatched number of kpoints in projections weights and connectivity list.')

    ! inputs 
    i=0
    ! 1) Ep
    i=i+1
    allocate(   Ep( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),Ep, mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 2) E(nbands,nkpts)
    i=i+1
    allocate(   E( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),E,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 3) tet(1:4,ntets)
    i=i+1
    allocate( tet( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tet,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 4) tetw(ntets)
    i=i+1
    allocate( tetw( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tetw,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 5) weight(nprojections,nbands,nkpts)
    i=i+1 
    allocate( weight( mxGetM(prhs(i)) , mxGetM(prhs(2)) , mxGetN(prhs(2)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),weight,mxGetM(prhs(i))*mxGetM(prhs(2))*mxGetN(prhs(2)))

    ! execute code
    D = get_pdos_tet(Ep,E,nint(tet),tetw,weight)

    ! outputs
    i=0
    ! output results
    i=i+1
    plhs(i) = mxCreateDoubleMatrix(mxGetM(prhs(5)),mxGetM(prhs(1))*mxGetN(prhs(1)),0)
    call mxCopyReal8ToPtr(D,mxGetPr(plhs(1)),mxGetM(prhs(5))*mxGetM(prhs(1))*mxGetN(prhs(1)))
end subroutine mexFunction
%}
            



            % compile everything
            
            fprintf('Building mex library:\n')
                if system([am_lib.FC,' ',am_lib.FFLAGS,' ',am_lib.DEBUG,' ',am_lib.LIBS,' -c ',flib,'.f90'])
                    error(' ... %s (failed)\n',flib)
                else
                    fprintf(' ... %s (succeeded)\n',flib);
                end
            
            fprintf('Building mex interfaces:\n')
            for i = 1:numel(f)
                if system([am_lib.FC,' ',am_lib.FFLAGS,' ',am_lib.DEBUG,' ',am_lib.LIBS,' ',flib,'.o ',f{i},'.f90 -o ',f{i},am_lib.EXT])
                    warning(' ... %s (failed)\n',f{i})
                else
                    fprintf(' ... %s (succeeded)\n',f{i});
                end
            end
        
        end
        
        function [K] = uc2ws(K,M)
            % uc2ws uses M real (reciprocal) lattice vectors to reduces K(1:3,:) vectors 
            % in cartesian (reciprocal) coordinates to the definiging Wigner-Seitz cell.

            if  and( am_lib.usemex , which('uc2ws_mex') )
                % use mex function if available
                K = uc2ws_mex(K,M,am_lib.eps);
            else
                % generate mesh
                G = M * [-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1;
                         -1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1;
                         -1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1];
                G2=sum(G.^2,1);
                % call engine
                for j = 1:size(K,2); K(:,j) = uc2ws_engine(K(:,j),G,G2,am_lib.eps); end
            end

            function [K] = uc2ws_engine(K,G,G2,tiny)
                go = true;
                while go; go=false;
                    for i = 1:26
                        P = 2*K.'*G(:,i)-G2(i);
                        if P > + tiny
                            K = K - G(:,i); go=true;
                        end
                    end
                end
            end
        end

        function [D] = get_dos_tet(Ep,E,tet,tetw)
            %
            if  and( am_lib.usemex , which('get_dos_tet_mex') )
                % use mex function if available
                D = get_dos_tet_mex(Ep,E,tet,tetw);
            else
                nEps = numel(Ep); D = zeros(nEps,1);
                for m = 1:nEps
                    D(m) = get_dos_tet_engine(Ep(m),E,tet,tetw);
                end
                D = reshape(D,size(Ep));
            end

            function [dosEp] = get_dos_tet_engine(Ep,E,tet,tetw)
                nbands = size(E,1);
                ntets = size(tet,2);
                dosEp = 0;
                for j = 1:ntets
                    for i = 1:nbands
                        % get energies at the verticies of the j-th tethedron in ascending order
                        Ec(1:4) = sort(E(i,tet(:,j)));
                        % calculate sum over k of the integrated charge
                        if     Ep<=Ec(1)
                            % Eq C1 from Blochl's PhysRevB.49.16223
                            % do nothing
                        elseif Ep<=Ec(2)
                            % Eq C2 from Blochl's PhysRevB.49.16223
                            dosEp = dosEp + tetw(j)*3.0*(Ep-Ec(1)).^2/(Ec(2)-Ec(1))/(Ec(3)-Ec(1))/(Ec(4)-Ec(1));
                        elseif Ep<=Ec(3)
                            % Eq C3 from Blochl's PhysRevB.49.16223
                            dosEp = dosEp + tetw(j)/(Ec(3)-Ec(1))/(Ec(4)-Ec(1))*(3.0*(Ec(2)-Ec(1))+6.0*(Ep-Ec(2))-3.0*(Ec(3)-Ec(1)+Ec(4)-Ec(2))/(Ec(3)-Ec(2))/(Ec(4)-Ec(2))*(Ep-Ec(2)).^2);
                        elseif Ep<=Ec(4)
                            % Eq C4 from Blochl's PhysRevB.49.16223
                            dosEp = dosEp + tetw(j)*(3.0*(Ec(4)-Ep).^2/(Ec(4)-Ec(1))/(Ec(4)-Ec(2))/(Ec(4)-Ec(3)));
                        elseif Ep>=Ec(4)
                            % Eq C1 from Blochl's PhysRevB.49.16223
                            % do nothing
                        end
                    end
                end
            end
        end

        function [D] = get_idos_tet(Ep,E,tet,tetw)
            %
            if  and( am_lib.usemex , which('get_idos_tet_mex') )
                % use mex function if available
                D = get_idos_tet_mex(Ep,E,tet,tetw);
            else
                nEps = numel(Ep); D = zeros(nEps,1);
                for m = 1:nEps
                    D(m) = get_dos_tet_engine(Ep(m),E,tet,tetw);
                end
                D = reshape(D,size(Ep));
            end

            function [idos] = get_dos_tet_engine(Ep,E,tet,tetw)
                nbands = size(E,1);
                ntets = size(tet,2);
                idos = 0;
                for j = 1:ntets
                    for i = 1:nbands
                        % get energies at the verticies of the j-th tethedron in ascending order
                        Ec(1:4) = sort(E(i,tet(:,j)));
                        % calculate sum over k of the integrated charge
                        if     Ep<=Ec(1)
                            % Eq C1 from Blochl's PhysRevB.49.16223
                            % do nothing
                        elseif Ep<=Ec(2)
                            % Eq A2 from Blochl's PhysRevB.49.16223
                            idos = idos + tetw(j)*(Ep-Ec(1)).^3/(Ec(2)-Ec(1))/(Ec(3)-Ec(1))/(Ec(4)-Ec(1));
                        elseif Ep<=Ec(3)
                            % Eq A3 from Blochl's PhysRevB.49.16223
                            idos = idos + tetw(j)/(Ec(3)-Ec(1))/(Ec(4)-Ec(1))*((Ec(2)-Ec(1)).^2+3.0*(Ec(2)-Ec(1))*(Ep-Ec(2))+3.0*(Ep-Ec(2)).^2.0-(Ec(3)-Ec(1)+Ec(4)-Ec(2))/(Ec(3)-Ec(2))/(Ec(4)-Ec(2))*(Ep-Ec(2)).^3.0);
                        elseif Ep<=Ec(4)
                            % Eq A4 from Blochl's PhysRevB.49.16223
                            idos = idos + tetw(j)*(1.0-(Ec(4)-Ep).^3.0/(Ec(4)-Ec(1))/(Ec(4)-Ec(2))/(Ec(4)-Ec(3)));
                        elseif Ep>=Ec(4)
                            % Eq A5 from Blochl's PhysRevB.49.16223
                            idos = idos + tetw(j);
                        end
                    end
                end
            end
        end

        function [D] = get_pdos_tet(Ep,E,tet,tetw)
            %
            if  and( am_lib.usemex , which('get_pdos_tet_mex') )
                % use mex function if available
                D = get_pdos_tet_mex(Ep,E,tet,tetw);
            else
                error('matlab version not yet implemented. mex is required for pdos');
            end
        end
        
    end
   
end








