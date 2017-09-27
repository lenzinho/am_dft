classdef am_dft

    properties (Constant)
        tiny      = 1E-4; % precision of atomic coordinates
        eps       = 1E-8; % numerical precision
        units_eV  = 0.06465555; % sqrt( [eV/Ang^2] * [1/amu] ) --> 0.06465555 [eV]
        units_THz = 98.22906;   % sqrt( [eV/Ang^2] * [1/amu] ) --> 98.22906 [THz=1/ps]
        units_GHz = 98229.06;   % sqrt( [eV/Ang^2] * [1/amu] ) --> 98229.06 [GHz=1/fs]

        usemex    = false;
        potdir    = '/Users/lenzinho/Linux/vasp.5.4/potcars/PBE.54/';
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

            import am_lib.*  am_dft.*

            % check inputs
            for f = {opts.fposcar,opts.fforce_position}
            if ~exist(f{:},'file'); error('File not found: %s',f{:}); end
            end

            % file name to save to/load from
            sname = sprintf('%s_%s_%0.2f_%0.2f.mat',...
                opts.fposcar, opts.fforce_position, opts.cutoff2, opts.cutoff3);
            if and(exist(sname,'file'),opts.continue); load(sname); else

                % get cells
                [uc,pc] = get_cell('poscar',opts.fposcar);
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
            %
            % defaults that should work for almost anything:
            % opts.continue=true;
            % opts.fposcar='POSCAR';
            % opts.feigenval='EIGENVAL';
            % opts.spdf=get_vasp('potcar:vrhfin:orbital');
            % opts.nskips=0;
            % opts.Ef=get_vasp('outcar:fermi');
            % opts.cutoff2=3;
            %

            import am_lib.* am_dft.*

            % check inputs
            for f = {opts.fposcar,opts.feigenval}
            if ~exist(f{:},'file'); error('File not found: %s',f{:}); end
            end

            % file name to save to/load from
            sname = sprintf('%s_%s_%02i_%0.2f_%0.2f_%s.mat',...
                opts.fposcar, opts.feigenval, opts.nskips, ...
                opts.Ef, opts.cutoff2, strjoin(opts.spdf,'+') );
            if and(exist(sname,'file'),opts.continue); load(sname); else

                % get cells
                [uc,pc] = get_cell('poscar',opts.fposcar);
                % load dft
                [dft]   = load_eigenval(opts.feigenval,opts.Ef);
                % construct sc
                [sc,sc.u2p,sc.p2u] = get_supercell(pc,diag([5,5,5])); sc.i2u = sc.p2u(pc.i2p); sc.u2i = pc.p2i(sc.u2p); 
                % get tb
                [tb,pp] = get_tb(pc,sc,dft,opts.cutoff2,opts.spdf,opts.nskips);
                % save results
                save(sname,'uc','pc','dft','tb','pp');

            end

            % plot results
                % plot correlation for dft vs bvk forces on atoms
                figure('color','white'); plot_tb_vs_dft(tb,dft); drawnow;
        end

    end
    
    
    % vasp
    
    methods (Static)

        % io
        
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

        function [dft]   = load_procar(Ef)

            import am_dft.* am_lib.*

            if nargin < 1; Ef = get_vasp('outcar:fermi'); end

            fprintf(' ... loading PROCAR'); tic;
                dft.nks    = get_vasp('procar:nkpts');
                dft.nbands = get_vasp('procar:nbands');
                dft.natoms = get_vasp('procar:natoms');
                dft.norbitals = get_vasp('procar:norbitals');
                dft.nspins = get_vasp('procar:ispin');
                dft.k = get_vasp('procar:kpoint').';
                dft.w = get_vasp('procar:weight').';
                dft.E = get_vasp('procar:energy'); dft.E = dft.E-Ef;
                dft.f = get_vasp('procar:occupation');
                dft.orbital = get_vasp('procar:orbital');
                dft.lmproj = get_vasp('procar:lmproj');
            fprintf('(%.f secs)\n',toc);

            fprintf('     %-15s = %-8.3f [eV] \n','fermi energy',Ef);
            fprintf('     %-15s = %i\n','nkpoints',dft.nks);
            fprintf('     %-15s = %i\n','nbands',dft.nbands);
            fprintf('     %-15s = %i\n','natoms',dft.natoms);
            fprintf('     %-15s = %i\n','norbitals',dft.norbitals);
            fprintf('     %-15s = %i\n','nspins',dft.nspins);

            % normalize projections by summing over atoms and characters?
            % dft.lmproj = dft.lmproj ./ sum(sum(dft.lmproj,2),3);

        end

        function [dos]   = load_doscar(Ef)
            
            import am_lib.*
            
            [str,~] = load_file_('DOSCAR');
                % get energies, density of states, and integration DOS
                t=sscanf(str{1},'%i'); natoms = t(1);
                t=sscanf(str{6},'%f'); nEs = round(t(3)); if nargin<1; Ef = t(4); end; 
                nspins = round((numel(strsplit(strtrim(str{7}),' '))-1)/2);
                t=sscanf(sprintf('%s\n',str{6+[1:nEs]}),'%f'); t=reshape(t,2*nspins+1,nEs).';
                dos.E = t(:,1)-Ef; dos.D = sum(t(:,2*[1:nspins]),2); dos.iD = sum(t(:,1+2*[1:nspins]),2);
            % proj(nEs,nspins,norbitals,natoms) ==> proj(nspins,norbitals,natoms,nbands,nkpts,nEs)
            norbitals = 9;
            proj = zeros(nEs,nspins,norbitals,natoms);
            for i = [1:natoms]
                t = sscanf(sprintf('%s\n',str{6+i*(1+nEs)+[1:nEs]}),'%f');
                t = reshape(t,1+nspins*9,nEs).'; % nEs x spin x orbitals
                t = reshape(t(:,2:end),nEs,nspins,norbitals);
                proj(:,:,:,i) = t;
            end
            % proj(nspins,norbitals,natoms,nbands,nkpts,nEs)
            dos.proj = permute(proj,[2,3,4,5,6,1]);
            % normalize projections
            n = sum_(dos.proj,[1,2,3,4,5]); n(eq_(n,0))=1; dos.proj=dos.proj./n; 
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

            import am_lib.* am_dft.*

            % count number of lines in file and check that all runs completed properly
            nlines = count_lines_(fbands); if mod(nlines,nbands)~=0; error('lines appear to be missing.'); end

            % open file and parse
            nsteps=nlines/nbands; fid=fopen(fbands); en=reshape(fscanf(fid,'%f'),nbands,nsteps); fclose(fid);
        end

        function [incar] = load_incar(fincar)
            import am_lib.* am_dft.*
            % load incar file into memory as string
            str = load_file_(fincar);
            % generate base incar and get tokens
            incar = generate_incar('base'); token_list = fields(incar);
            for i = 1:numel(token_list)
                token = token_list{i};
                % parse
                    t = extract_token_(strtrim(str),token);
                    t = strsplit(t,'%');
                    t = strrep(t{1},'=','');
                    t = strtrim(t);
                incar.(token) = t;
            end
        end

        function incar   = generate_incar(flag,nondefault_incar_opts_)
            
            import am_dft.*
            
            if isstruct(flag)
                incar = flag;
            else
                switch flag
                    case 'base'
                        % lepsilon
                        opts_ = { ...
                        'istart'   , '' , 'icharg' , '' , ... % restarting
                        'ispin'    , '' , 'gga'    , '' , 'ivdw'   , '' , ... % exchange correlations stuff
                        'encut'    , '' , 'prec'   , '' , 'algo'   , '' , 'ialgo'  , '' , 'amix'   , '' , 'ediff'    , '' , 'nelmdl' , '' , ...  % precision/convergence stuff
                        'nbands'   , '' , 'ismear' , '' , 'sigma'  , '' , 'lorbit' , '' , ...
                        'ldau'     , '' ,'ldautype', '' , 'ldaul' , '' , 'ldauu'   , '' , 'ldauj'     , '' , 'ldauprint' , '' , 'lmaxmix'   , '' , ... ... % DFT+U stuff
                        'ibrion'   , '' , 'isif'   , '' , 'ediffg', '' , 'potim'   , '' , 'tebeg'     , '' , 'nblock'    , '' , 'nsw'     , '' , 'smass'  , '' , 'addgrid' , '' , ... % aimd stuff
                        'ncore'    , '' , 'npar'   , '' , 'nsim'   , '' , ... % cluster stuff 
                        'lreal'    , '' , 'lwave'     , '' , 'lcharg'  , '' , 'lvtot'  , '' ... % save stuff
                        };
                        incar = generate_incar('',opts_);
                    case 'defaults'
                        opts_ = { ...
                        'istart' , '0'       , 'icharg'    , '2'       , 'gga'     , 'AM'      , 'prec'   , 'Accurate' , ...
                        'ediff'  , '1e-9'    , 'algo'      , 'VeryFast', 'encut'  , '650'      , ...
                        'ismear' , '0'       , 'sigma'     , '0.2'     , ...
                        'isif'   , '3'       , 'lorbit'    , '11'      , 'nblock'  , '1'       , 'addgrid', '.FALSE.'   , ...
                        'lreal'  , '.FALSE.' , 'lwave'     , '.FALSE.' , 'lcharg'  , '.FALSE.' , 'lvtot'  , '.FALSE.' };
                        incar = generate_incar('base',opts_);
                    case 'rlx'
                        % Setting the convergence criterion based on forces
                        % being below < 0.005 eV/Ang (based on G. Kresse's
                        % PHYSICAL REVIEW B 78, 104116 2008) works only if
                        % there is at least one ionic degree of freedom.
                        % For NaCl all ions are constrained at Wyckoff
                        % positions; therefore, ediffg = 10*ediff is more
                        % appropriate.
                        opts_ = { ...
                            'ediffg' , '1e-8'    , 'ibrion' , '2'       , 'potim'  , '0.5'    , ...
                            'nsw'    , '100'     , 'smass'  , '0'       , 'addgrid', '.TRUE.' , ...
                        };
                        incar = generate_incar('defaults',opts_);
                    case 'scf'
                        opts_ = { ...
                            'ibrion' , '-1'       , 'nsw'    , '0'       , 'lorbit' , '11'      , ...
                            'addgrid', '.FALSE.'  , ... % takes too long with addgrid
                            'lwave'  , '.TRUE.'   , 'lcharg' , '.TRUE.'  , 'lvtot'  , '.TRUE.'  , ...
                            };
                        incar = generate_incar('defaults',opts_);
                    case {'nscf','evk'}
                        opts_ = { ...
                            'icharg' , '11'       ,  ...
                            'lwave'  , '.FALSE.'  , 'lcharg' , '.FALSE.' , 'lvtot'  , '.FALSE.'  , ...
                            };
                        incar = generate_incar('scf',opts_);
                    case 'eps'
                        opts_ = {'lepsilon', '.TRUE.'};
                        incar = generate_incar('nscf',opts_);
                    case 'opt'
                        opts_ = {'loptics','.TRUE.','nedos','10000'};
                        incar = generate_incar('nscf',opts_);
                    case 'aimd'
                        % similar to rlx but with ibrion=0, a lower energy convergence, no convergence on ionic motion, nose thermostat, and more steps 
                        opts_ = {...
                            'ibrion' , '0' , 'ediff','1e-6', 'ediffg', '', 'smass', '0', 'nsw', '1000', ...
                            };
                end
            end
            
            if nargin > 1
                n = numel(nondefault_incar_opts_);
                if n > 1; for i = 1:2:n
                    incar.(nondefault_incar_opts_{i})=nondefault_incar_opts_{i+1};
                end; end
            end
            
        end

        
        % for automating vasp simulations

        function           write_poscar(uc,fposcar)

            import am_lib.* am_dft.*

            % write_poscar(uc,fposcar)
            if nargin < 2; fposcar='POSCAR'; end

            if det(uc.bas)<0; error('vasp requires basis vectors with positive determinants'); end

            n = size(uc.tau,3);
            for i = 1:n
                % file name
                fname = sprintf('%s',fposcar);
                if n ~= 1; fname = sprintf('%s%s_%06i',fname,fposcar,i); end
                % open file
                fid=fopen(fname,'w');
                    % header
                    fprintf(fid,'%s',get_formula(uc));
                    if n ~= 1; fprintf(fid,' %i of %i',i,n); end
                    fprintf(fid,'\n');
                    % print body
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

        function           write_potcar(uc,potdir)
            import am_dft.*

            % default
            if nargin<2;potdir=am_dft.potdir; end
            % database
            potcar_database = {...
                '  ' ,'  ' ,'Li_sv' ,'Be_sv' ,'  ' ,'  ' ,'  ' ,'O_s_GW' , ... %  h     he    li    be    b     c     n     o
                '  ' ,'  ' ,'Na_pv' ,'Mg_pv' ,'Al_GW' ,'Si' ,'  ' ,'  ' , ... %  f     ne    na    mg    al    si    p     s
                '  ' ,'  ' ,'K_sv' ,'Ca_sv' ,'Sc_sv_GW' ,'Ti_pv' ,'  ' ,'  ' , ... %  cl    ar    k     ca    sc    ti    v     cr
                'Mn_GW' ,'Fe_GW' ,'  ' ,'Ni_sv_GW' ,'  ' ,'  ' ,'Ga_sv_GW' ,'  ' , ... %  mn    fe    co    ni    cu    zn    ga    ge
                '  ' ,'  ' ,'  ' ,'  ' ,'Rb_sv' ,'Sr_sv' ,'Y_sv' ,'  ' , ... %  as    se    br    kr    rb    sr    y     zr
                '  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' , ... %  nb    mo    tc    ru    rh    pd    ag    cd
                'In_sv_GW' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'Cs_sv' ,'Ba_sv' , ... %  in    sn    sb    te    i     xe    cs    ba
                'La_GW' ,'Ce_3' ,'Pr_3' ,'Nd_3' ,'Pm_3' ,'Sm_3' ,'Eu_3' ,'Gd_3' , ... %  la    ce    pr    nd    pm    sm    eu    gd
                'Tb_3' ,'Dy_3' ,'Ho_3' ,'Er_3' ,'Tm_3' ,'Yb_3' ,'Lu' ,'  ' , ... %  tb    dy    ho    er    tm    yb    lu    hf
                '  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' , ... %  ta    w     re    os    ir    pt    au    hg
                '  ' ,'  ' ,'Bi_GW' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' , ... %  tl    pb    bi    po    at    rn    fr    ra
                '  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' , ... %  ac    th    pa    u     np    pu    am    cm
                '  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' , ... %  bk    cf    es    fm    md    no    lr    rf
                '  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' };    %  db    sg    bh    hs    mt    ds    rg    uub
            Z = get_atomic_number(uc.symb); nZs = numel(Z);
            for i = 1:nZs
                if isempty(potcar_database{Z(i)}); error('potcar not defined'); end
                if i == 1
                    cmd = sprintf('cat %s/%s/POTCAR > POTCAR' ,potdir,potcar_database{Z(i)});
                else
                    cmd = sprintf('cat %s/%s/POTCAR >> POTCAR',potdir,potcar_database{Z(i)});
                end
                system(cmd);
            end
        end

        function           write_incar(incar,fincar)
            if nargin < 2; fincar='INCAR'; end
            token_list = fields(incar); fclose('all');
            fid = fopen(fincar,'w');
                for i = 1:numel(token_list)
                    token = token_list{i};
                    if ~isempty(incar.(token))
                        fprintf(fid,'%-10s ',token);
                        fprintf(fid,'= %-10s',incar.(token));
                        fprintf(fid,' %% %-10s',get_comments_(token));
                        fprintf(fid,'\n');
                    end
                end
            fclose(fid);

            function c = get_comments_(token)
                switch token
                % starting
                case 'istart'; 		c='(0) new  (1) cont (2) samecut';
                case 'icharg'; 		c='(1) file (2) atom (3) const <scf | nscf> (11)-chgcar';
                % electronic
                case 'gga';    		c='exchannge correlation function';
                case 'nelmdl'; 		c='non-selfconsistent iterations: (<0) only once, (>0) at every ionic step';
                case 'outcar:ispin';  		c='(1) non-spin-polarized (2) spin-polarized';
                case 'prec';   		c='Accurate, Normal';
                case 'ivdw';   		c='vdW';
                case 'algo';   		c='(VeryFast) RMM-DIIS (Fast) Davidson/RMM-DIIS (Normal) Davidson';
                case 'ialgo';  		c='(38) Davidson (48) RMM-DIIS algorithm';
                case 'outcar:nbands'; 		c='number of bands';
                case 'encut';  		c='cutoff';
                case 'ismear'; 		c='(0) gauss (1) mp (-5) tetrah';
                case 'sigma';  		c='maximize with entropy below 0.1 meV / atom';
                case 'ediff';  		c='toten convergence';
                % dft+U
                case 'ldau'; 		c='Activates DFT+U calculation';
                case 'ldautype'; 	c='(1) Liechtenstein (2) Dudarev Ueff = U - J (3) adds Ueff on the atomic sites';
                case 'ldaul'; 	    c='one entry per species: (>0) quantum number to which U is applied (-1) U not applied to this species';
                case 'ldauu'; 	    c='(for LDAUTYPE = 1,2) Hubbard U per species  (for LDAUTYPE = 3) Ueff acting on spin up';
                case 'ldauj'; 	    c='(for LDAUTYPE = 1,2) Exchange J per species (for LDAUTYPE = 3) Ueff acting on spin down';
                case 'ldauprint';   c='(0) silent (1) write occupancy matrix OUTCAR (2) occupancy matrix OUTCAR and potential matrix STDOUT';
                case 'lmaxmix'; 	c='(4) d-electrons, (6) f-electrons';
                % ionic
                case 'ediffg';      c='(<0) relaxation will stop if all forces are less than |EDIFFG| (>0) relaxation stops if energy change is less than EDIFFG';
                case 'tebeg';  		c='temperature of the MD';
                case 'ibrion'; 		c='(-1) no update (0) MD (1) RMM-DIIS quasi-Newton (2) conjugate-gradient';
                case 'potim';  		c='Timestep in femtoseconds';
                case 'isif';   		c='Relax (2) ions (3) ions,volume,shape (4) ions,shape';
                case 'nblock'; 		c='write after every iteration';
                case 'nsw';    		c='number of ionic steps';
                case 'smass';  		c='(-3) NVE (=>0) Nose thermostat';
                % parallelization
                case 'ncore';  		c='approx SQRT( number of cores)';
                case 'npar';   		c='Parallelization';
                case 'nsim';   		c='Number of bands which are optimized by RMMS-DIIS algorithm simultaneously';
                % oribtal info
                case 'lorbit'; 		c='(11) lm-proj DOS and PROCAR (0) DOS';
                % projection scheme
                case 'lreal';  		c='(.FALSE.) for <20 atoms (Auto) for >20 atoms';
                % save
                case 'lwave';  		c='write wavefunctions';
                case 'lcharg'; 		c='write charge density';
                case 'lvtot';  		c='write local charge potential';
                otherwise;     		c='';
                end
            end
        end

        function           write_kpoints(bz,fkpoints)
            if nargin < 2; fkpoints='KPOINTS'; end
            if     isnumeric(bz)
                if all(bz==1)
                    % simple list [nx,ny,nz] kpoint point density
                    fid = fopen(fkpoints,'w');
                        fprintf(fid,'%s \n','Automatically generated mesh');
                        fprintf(fid,'%s \n','0');
                        fprintf(fid,'%s \n','Gamma');
                        fprintf(fid,' %i %i %i \n',[1 1 1]);
                        fprintf(fid,' %i %i %i \n',[0 0 0]);
                   fclose(fid);
                else
                    % simple list [nx,ny,nz] kpoint point density
                    fid = fopen(fkpoints,'w');
                        fprintf(fid,'%s \n','Automatically generated mesh');
                        fprintf(fid,'%s \n','0');
                        fprintf(fid,'%s \n','Monkhorst-Pack');
                        fprintf(fid,' %i %i %i \n',ceil(bz));
                        fprintf(fid,' %i %i %i \n',[0 0 0]);
                   fclose(fid);
                end
            elseif isstruct(bz)
                % still need to implement.
            end
        end

        function           write_qsub(hardware,account,nnodes,hrs,job)
            % defaults
            if nargin < 2; account = []; end
            if nargin < 3; nnodes = 1; end
            if nargin < 4; hrs = 12; end
            if nargin < 5; job = 'job'; end

            fid = fopen('QSUB','w');
                switch hardware
                    case 'taub'
                        fprintf(fid,'#PBS -l nodes=%i:ppn=12\n',nnodes);
                        fprintf(fid,'#PBS -N %s \n',job);
                        fprintf(fid,'#PBS -l walltime=%i:00:00\n',hrs);
                        fprintf(fid,'#PBS -q cse \n');
                        fprintf(fid,'#PBS -V \n');
                        fprintf(fid,'module load openmpi/1.4-gcc\n');
                        fprintf(fid,'module load intel/11.1\n');
                        fprintf(fid,'cd \\${PBS_O_WORKDIR}\n');
                        fprintf(fid,'mpiexec /home/amei2/vasp.5.3.3 > output\n');
                    case 'trio'
                        if isempty(account); error('Account is required for running on triolith'); end
                        fprintf(fid,'#!/bin/bash \n');
                        fprintf(fid,'#SBATCH -N %i \n',nnodes);
                        fprintf(fid,'#SBATCH -t %i:00:00 \n',hrs);
                        fprintf(fid,'#SBATCH -J %s \n',job);
                        fprintf(fid,'#SBATCH --exclusive \n');
                        fprintf(fid,'#SBATCH -A %s \n',account);
                        fprintf(fid,'module load vasp/5.4.4-18Apr17 \n');
                        fprintf(fid,'mpprun vasp_std > output \n');
                        fprintf(fid,'#mpprun vasp-gamma > output \n');
                end
            fclose(fid);
        end

        function [x]  =    get_vasp(flag)
            import am_dft.*
            
            % check file is present
            flist = {'PROCAR','OUTCAR','vasprun.xml'};
            for i = 1:numel(flist)
                if     contains(flag,lower(flist{i}(1:6)))
                    if exist(flist{i},'file')~=2
                        fprintf('\n');
                        fprintf('\n');
                        fprintf('ERROR: Trying to read from %s when it does not exist',flist{i});
                        fprintf('\n');
                        fprintf('\n');
                        x = [];
                        return;
                    end
                end
            end
            
            switch flag
                % OUTCAR
            	case 'outcar:site_occupancy' % for DFT+U
                    [~,m] = system('awk ''/onsite density matrix/'' OUTCAR'); m = numel(strfind(m,'matrix'));
                    [~,x] = system('awk ''/ o = /{ print $3 }'' OUTCAR');     x = sscanf(x,'%f'); x = reshape(x,[],m);
                case 'outcar:toten'
                    [~,x] = system('awk ''/TOTEN/ { print $5 }'' OUTCAR');  x = sscanf(x,'%f');
                case 'outcar:fermi'
                    [~,x] = system('awk ''/E-fermi/ { print $3 }'' OUTCAR');  x = sscanf(x,'%f');
                case 'outcar:pressure'
                    [~,x] = system('awk ''/pressure/ { print $4 }'' OUTCAR'); x = sscanf(x,'%f');
                case 'outcar:temperature'
                    [~,x] = system('awk ''/EKIN_LAT/ { print $6 }'' OUTCAR'); x = sscanf(x,'%f');
                case 'outcar:ispin'
                    [~,x] = system('awk ''/ISPIN/ { print $3 }'' OUTCAR'); x = sscanf(x,'%f');
                case 'outcar:nbands'
                    [~,x] = system('awk ''/NBANDS/ {print $15}'' OUTCAR'); x = sscanf(x,'%f');
                case 'outcar:kpoint'
                    [~,x] = system('awk ''/ plane waves: / { print $4 "  " $5 "  " $6}'' OUTCAR'); x = sscanf(x,'%f'); x = reshape(x,3,[]).';
                case 'outcar:gamma'
                    [~,x] = system('awk ''/GAMMA/ {print $4}'' OUTCAR'); x = sscanf(x,'%f');
                case 'outcar:niterations'
                    [~,x] = system('awk ''/Iteration/'' OUTCAR'); x = numel(strfind(x,'Iteration'));
                case 'outcar:magnetization'
                    [~,x] = system('awk ''/magnetization/'' OUTCAR | awk ''/electron/ { print $6 }'''); x = sscanf(x,'%f');
                case 'outcar:total_charge' % obtained when iorbit=11 is specified
                    % x( natoms , {s,p,d,f} )
                    natoms = get_vasp('vasprun:natoms');
                    [~,x] = system(sprintf('grep -A %i ''total charge'' OUTCAR | tail -n %i | awk ''{print $2 " " $3 " " $4 " " $5}''',natoms+3,natoms)); x = sscanf(x,'%f'); x = reshape(x,[],natoms).';
                    
                % POTCAR
                case 'potcar:zval'  % valence on each atom in potcar
                    [~,x] = system('awk ''/ZVAL/{print $6}'' POTCAR'); x = sscanf(x,'%f');
                case 'potcar:vrhfin:configuration' % ionic configuration of each atom (orbitals + occupation)
                    [~,x] = system('grep ''VRHFIN'' POTCAR | sed -E ''s/([0-9])/\1,/g'' | sed -E ''s/ //g'' | sed -E ''s/,/ /g'' | sed -n -E ''s/^.*://p'' '); 
                    x = strtrim(strsplit(strtrim(x),'\n'));
                case 'potcar:vrhfin:orbital' % orbitals of each atom 
                    [~,x] = system('grep ''VRHFIN'' POTCAR | sed -E ''s/([0-9])//g'' | sed -E ''s/ //g'' | sed -n -E ''s/^.*://p'' '); 
                    x = strtrim(strsplit(strtrim(x),'\n'));
                    
                % PROCAR
                case 'procar:ispin'
                    [~,x] = system('awk ''/# of k-points:/'' PROCAR'); x = numel(strfind(x,'ions'));
                case 'procar:nkpts'
                    [~,x] = system('grep -m 1 "# of k-points:" PROCAR | awk ''{print $4}'''); x = sscanf(x,'%f');
                case 'procar:nbands'
                    [~,x] = system('grep -m 1 "# of k-points:" PROCAR | awk ''{print $8}'''); x = sscanf(x,'%f');
                case 'procar:natoms'
                    [~,x] = system('grep -m 1 "# of k-points:" PROCAR | awk ''{print $12}'''); x = sscanf(x,'%f');
                case 'procar:norbitals'
                    [~,x] = system('grep -m 1 py PROCAR'); x = numel(strsplit(strtrim(x))) - 2;
                case 'procar:orbital'
                    [~,x] = system('grep -m 1 py PROCAR'); x = strsplit(strtrim(x)); x=x(2:(end-1));
                case 'procar:kpoint'
                    [~,x] = system('grep weight PROCAR | sed ''s/-/ -/g'' | awk ''{print $5 " " $6 " " $7}'''); x = sscanf(x,'%f'); x = reshape(x,3,[]).';
                case 'procar:weight'
                    [~,x] = system('grep weight PROCAR | grep ''k-point'' | sed ''s/-/ -/g'' | awk ''{print $10}''' ); x = sscanf(x,'%f');
                case 'procar:energy'
                    nkpts     = get_vasp('procar:nkpts');
                    nbands    = get_vasp('procar:nbands');
                    nspins    = get_vasp('procar:ispin');
                    [~,x] = system('awk ''/energy/{print $5}'' PROCAR'); x = sscanf(x,'%f');
                    % E(nbands,nkpts,nspins)
                    x = reshape(x,nbands,nkpts,nspins);
                case 'procar:occupation'
                    nkpts     = get_vasp('procar:nkpts');
                    nbands    = get_vasp('procar:nbands');
                    nspins    = get_vasp('procar:ispin');
                    [~,x] = system('awk ''/energy/{print $8}'' PROCAR'); x = sscanf(x,'%f');
                    % E(nbands,nkpts,nspins)
                    x = reshape(x,nbands,nkpts,nspins);
                case 'procar:lmproj'
                    natoms    = get_vasp('procar:natoms');
                    norbitals = get_vasp('procar:norbitals');
                    nkpts     = get_vasp('procar:nkpts');
                    nbands    = get_vasp('procar:nbands');
                    nspins    = get_vasp('procar:ispin');
                    switch norbitals
                        case 1+3+5
                            [~,x]=system(sprintf('grep -A%i py PROCAR | grep -v py | awk ''{print $2 " " $3 " " $4 " " $5 " " $6 " " $7 " " $8 " " $9 " " $10}''',natoms));
                        case 1+3+5+7
                            [~,x]=system(sprintf('grep -A%i py PROCAR | grep -v py | awk ''{print $2 " " $3 " " $4 " " $5 " " $6 " " $7 " " $8 " " $9 " " $10 " " $11 " " $12 " " $13 " " $14 " " $15 " " $16 " " $17}''',natoms));
                    end
                    x = strrep(x,'--',''); x = sscanf(x,'%f'); x = reshape(x,norbitals,natoms,nbands,nkpts,nspins); % last is ispin
                    % lmproj(nspins,norbitals,nions,nbands,nkpts)
                    x = permute(x,[5,1,2,3,4]);
                % VASPRUN
                case 'vasprun:scf_energy'
                    [~,x] = system('grep ''e_fr_energy'' vasprun.xml | sed ''s/\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*/ 1e8 /g'' | awk ''{print $3}'''); x = sscanf(x,'%f');
                case 'vasprun:natoms'
                    [~,x] = system('grep ''<atoms>'' vasprun.xml | awk ''{print $2}'''); x = sscanf(x,'%f');
            end
        end
            
        function           vasp_fix_input(flags)
            
            import am_dft.* am_lib.*
            
            % no flag
            if nargin < 1; flags=''; end
            
            % defaults
            if contains(flags,'dryrun'); isdryrun = true; else; isdryrun=false; end

            % PRIMITIVE BASIS UNREASONABLY LARGE
            uc = load_poscar('POSCAR'); d = uc.natoms/abs(det(uc.bas*0.1));
            if d<0.001 
                fprintf('DIR: %s \n',pwd); 
                fprintf('       ERROR: at least of the lattice vectors is too long. \n');
                fprintf('       bas: %10.5f %10.5f %10.5f \n',uc.bas(:,1));
                fprintf('            %10.5f %10.5f %10.5f \n',uc.bas(:,2));
                fprintf('            %10.5f %10.5f %10.5f \n',uc.bas(:,3));
                target_density = 50; % [atoms/nm^3 ~ 5E22 atoms/cm3]
                uc.bas = (d/target_density)^(1/3)*uc.bas;
                fprintf('       FIX: scaling lattice vectors to get a reasonable density of ~ %i atoms/nm^3 \n',target_density)
                fprintf('       bas: %10.5f %10.5f %10.5f \n',uc.bas(:,1));
                fprintf('            %10.5f %10.5f %10.5f \n',uc.bas(:,2));
                fprintf('            %10.5f %10.5f %10.5f \n',uc.bas(:,3));
                if ~isdryrun
                    write_poscar(uc);
                end
            end

            % SCF DIVERGED
            b = get_vasp('vasprun:scf_energy');
            if ~isempty(b) && any(abs(b)>1E7)
                fprintf('DIR: %s \n',pwd);
                fprintf('       ERROR: electronic cycle divreged (d eps = %g). \n',b(end));
                % fix is to try a lower AMIX and set SCF algo to Normal
                incar = load_incar('INCAR');
                % reduce amix by 1/2
                if isfield(incar,'amix') && ~isempty(incar.amix); amix = str2num(incar.amix); else; amix = 0.4; end
                amix = amix/2;
                incar = generate_incar(incar,{'algo','Normal','amix',num2str(amix)});
                % print the fix
                fprintf('       FIX: set algo = normal and reduce amix to %0.2f\n',amix)
                if ~isdryrun
                    write_incar(incar);
                end
                return;
            end
            
            % NOT ENOUGH BANDS
            if grepcontains_('Your highest band is occupied at some k-points!','output')
                fprintf('DIR: %s \n',pwd); 
                nbands = get_vasp('outcar:nbands');
                fprintf('       ERROR: not enough bands (nbands = %g). \n',nbands);
                % read incar
                incar = load_incar('INCAR');
                if isfield(incar,'ncore'); ncore = str2double(incar.ncore); else; ncore = 1; end
                % update nbands
                nbands = ceil(nbands*1.2); nbands = ncore-mod(nbands,ncore)+nbands;
                fprintf('       FIX: increase nbands by ~1.2x to %i\n',nbands)
                % fix is to try a lower AMIX and set SCF algo to Normal
                incar = generate_incar(incar,{'nbands',nbands});
                if ~isdryrun
                    write_incar(incar);
                end
            end

            % NONHERMITIAN MATRIX
            if grepcontains_('WARNING: Sub-Space-Matrix is not hermitian!','output')
                fprintf('DIR: %s \n',pwd);
                fprintf('       ERROR: Sub-Space-Matrix is not hermitian. \n');
                fprintf('       This is usually caused by some other problem. Investigate the output manually. \n')
            end

            
            function b = grepcontains_(str,file)
                [~,b]=system(sprintf('grep ''%s'' %s',str,file));
                b = ~isempty(b);
            end
        end
        
        
    end
    
    
    % core

    methods (Static)
        
        % test suite
        
        function test
            
            import am_dft.* am_lib.*
            
            % generate and test point groups
            for i = 1:32; pg_id(i) = identify_pointgroup(generate_pg(i,false)); end
            criteria = pg_id==[1:32];
            test_(all(criteria),'point group generation',sprintf('failed to generate pointgroups:%s',sprintf(' %i',find(criteria))))
            
            % test winger j=1 rotation matrices
            R = generate_pg(32,false);
            W = get_wigner(1,R,'tesseral');
            criteria = squeeze(all(all(eq_(W,R),1),2));
            test_(all(criteria),'wigner SO(3)',sprintf('failed to generate consistent SO(3) rotations:%s',sprintf(' %i',find(criteria))))
            
%             % test winger rotation axes
%             R = generate_pg(32,false);
%             W = get_wigner(1,R,'spherical');
%             (get_wigner_axis(DG)~=0)-(round_(R_axis_(R))~=0)
%             criteria = squeeze(all(all(eq_(W,R),1),2));
%             test_(all(criteria),'wigner SO(3)',sprintf('failed to generate consistent SO(3) rotations:%s',sprintf(' %i',find(criteria))))
            
            function test_(logical,test_name,fail_msg)
                if logical
                    fprintf('      %s: pass\n',test_name);
                else
                    fprintf('      %s: %s\n',test_name,fail_msg);
                end 
            end
        end
        


        % symmetry

%         function [sg,pg]      = get_groups(pc, tol)
%
%             % getting point groups
%             [T,H,S,R]    = get_symmetries(pc, tol)
%
%             fprintf(' ... getting group properties')
%
%             % get multiplication table
%             [MT,E,I] = get_multiplication_table(S);
%             if all(all(MT-MT.'==0))
%                 fprintf(', abelian' );
%             else
%                 fprintf(', non-abelian' );
%             end
%
%             % get generators
%             G = get_generators(MT);
%             fprintf(', %i generators (', numel(G)); fprintf(' %i', G); fprintf(')');
%
%             % character table
%
%         end

        function [T,H,S,R]    = get_symmetries(pc, tol)
            % [T,H,S,R] = get_symmetries(pc, tol=am_dft.tiny)
            % T = all possible translations which restore the crystal to iteself
            % H = holohogries (all possible rotations which restore the bravais lattice onto iteself)
            % S = space group symmetries
            % R = point group symmetries

            import am_lib.* am_dft.*

            % set default numerical tolerance
            if nargin < 2; tol = am_lib.tiny; end

            % define function to check first two dimensions
            check_ = @(A) all(all(abs(A)<tol,1),2);

            % define function to sort atoms and species into a unique order (reference)
            X_ = @(species,tau) sortc_([species;mod_(tau)]); X = X_(pc.species,pc.tau(:,:,1));

            % get vectors that preserve periodic boundary conditions
            N = 1; V=mod_(pc.tau(:,pc.species==pc.species(N))-pc.tau(:,N), tol);  nVs=size(V,2); ex_=false(1,nVs);
            for j = 1:nVs; ex_(j) = check_( X_(pc.species,pc.tau(1:3,:,1)-V(:,j))-X ); end
            T=[V(:,ex_),eye(3)]; T=T(:,rankc_(normc_(T)));

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
                if check_( X_(pc.species,H(:,:,i)*pc.tau+V(:,j)) - X ); nSs=nSs+1; S(1:3,1:4,nSs)=[ H(:,:,i), V(:,j) ]; end
            end; end; S = S(:,:,1:nSs);

            % set identity first
            id = member_(flatten_(eye(4)),reshape(S,4^2,[])); S(:,:,[1,id])=S(:,:,[id,1]);

            if nargout == 3; return; end

            % get point symmetries
            R  = reshape(uniquec_( reshape(S(1:3,1:3,:),[9,nSs]) ),3,3,[]);

            % set identity first
            id = member_(flatten_(eye(3)),reshape(R,3^2,[])); R(:,:,[1,id])=R(:,:,[id,1]);
        end

        function [MT,E,I]     = get_multiplication_table(S,tol)
            % [MT,E,I] = get_multiplication_table(S)
            % get multiplication table: S(:,:,i)*S(:,:,j) = S(:,:,MT(i,j))
            
            import am_lib.* am_dft.*

            if nargin<2; tol=am_dft.tiny; end
            
            s = size(S);
            
            if   s(1)==1
                % seitz operator combined with permutation (represented as a two-part cell)
                
                s = size(S{1});
                
                if     s(1)==4 && s(2)==4
                    % seitz operator (applies mod to translational components)
                    md_ = @(X) [X(1:12,:);mod_(X(13:15,:),tol);X(16:end,:)];
                    rs_ = @(X) md_(reshape(X,s(1)*s(2),[]));
                    nsyms = s(3);

                    ref = [rs_(S{1});reshape(S{2},size(S{2},1),[])];
                    opr = [rs_(matmulp_(S{1},permute(S{1},[1,2,4,3]))); reshape(operm_(S{2},S{2}),size(S{2},1),[])];

                    MT = reshape( member_( opr , ref, tol), s(3), s(3)); 
                elseif s(1)==3 && s(2)==3
                    % point operator
                    rs_ = @(X) reshape(X,s(1)*s(2),[]);
                    nsyms = s(3);

                    ref = [rs_(S{1});reshape(S{2},size(S{2},1),[])];
                    opr = [rs_(matmulp_(S{1},permute(S{1},[1,2,4,3]))); reshape(operm_(S{2},S{2}),size(S{2},1),[])];

                    MT = reshape( member_( opr , ref, tol), s(3), s(3));
                end
                
                
            elseif s(1)==4 && s(2)==4
                % seitz operator (applies mod to translational components)
                md_ = @(X) [X(1:12,:);mod_(X(13:15,:),tol);X(16:end,:)];
                rs_ = @(X) md_(reshape(X,s(1)*s(2),[]));
                nsyms = s(3);
                
                MT = reshape( member_( rs_(matmulp_(S,permute(S,[1,2,4,3]))), rs_(S), tol) , nsyms, nsyms);

            elseif s(1)==3 && s(3)==1
                % translation operator
                S = mod_(S);
                nsyms = s(2);
                
                MT = reshape( member_( mod_(reshape(osum_(S, S, 2),3,[]), tol), S, tol), nsyms, nsyms);
            else
                % point operator
                % s(1)==3 && s(2)==3
                % but also others...
                rs_ = @(X) reshape(X,s(1)*s(2),[]);
                nsyms = s(3);

                MT = reshape( member_( rs_(matmulp_(S,permute(S,[1,2,4,3]))), rs_(S), tol ) , nsyms, nsyms);
            end

            if any(sum(MT,1)~=sum([1:nsyms])) || any(sum(MT,2)~=sum([1:nsyms]))
                error('MT is incorrect. Check for mistakes in the symmetry and ensure that the symmetry is in the primitive basis');
            end

            % find identity
            if nargout>1; E = find(all(MT==[1:nsyms].',1)); end

            % get inverse indicies
            if nargout>2; I = [MT==E]*[1:nsyms].'; end
        end
        
        function [MT]         = relabel_multiplication_table(MT,fwd)
            % relabel multiplication table
            rev(fwd) = [1:size(MT,1)]; MT = rev(MT(fwd,fwd));
        end
        
        function [G,u]        = get_generators(MT)
            % [G,u] = get_generators(MT)
            % u reindexes symmetries based on the generators
            import am_lib.*
            % [gen_id] = get_generators(MT)
            nsyms = size(MT,1); flatten_ = @(x) x(:);
            for ngens = 1:nsyms
                % generate a list of all possible generators (ignores generator order for speed)
                genlist = nchoosek_(nsyms,ngens);
                % loop over the list one set of generators at a time, checking whether the entire multiplication table can be generated
                for i = 1:size(genlist,2)
                    % initialize, include identity
                    x = genlist(:,i);
                    % expand multiplication table
                    while true
                        u = unique([flatten_(MT(x,x));x],'stable');
                        switch numel(u)
                            case nsyms
                                % found!
                                G = genlist(:,i);
                                return;
                                %
                            case numel(x)
                                % exausted all possibiltiies
                                break;
                            otherwise
                                % update
                                x = u;
                        end
                    end
                end
            end
        end

        function [Glist,ulist]= get_generators_list(MT,ngens)
            % returns all possible sets of generators with ngens elements
            import am_lib.*
            % [gen_id] = get_generators(MT)
            nsyms = size(MT,1); flatten_ = @(x) x(:);
            % generate a list of all possible generators (incldues order of generator for completness)
            Glist = perm_norep_(nsyms,ngens); ngenlists = size(Glist,2);
            ex_ = false(1,ngenlists); ulist = zeros(nsyms,ngenlists);
            % loop over the list one set of generators at a time, checking whether the entire multiplication table can be generated
            for i = 1:ngenlists
                % initialize
                x = Glist(:,i);
                % expand multiplication table
                while true
                    u = unique([flatten_(MT(x,x));x],'stable');
                    switch numel(u)
                        case nsyms
                            % found!
                            ex_(i) = true;
                            ulist(:,i) = u;
                            break;
                            %
                        case numel(x)
                            ex_(i) = false;
                            % exausted all possibiltiies
                            break;
                        otherwise
                            % update
                            x = u;
                    end
                end
            end
            Glist = Glist(:,ex_); ulist = ulist(:,ex_);
        end

        function                plot_cayley_graph(MT)
            %
            [G] = get_generators(MT);
            i=[]; j=[]; x = MT(1,:);
            for k = 1:numel(G)
                i = [i;flatten_(x)];
                j = [j;flatten_(MT(G(k),x))];
            end
            plot(digraph(i,j));
        end

        function [IR]         = get_irreducible_representations(S)
            % s2c = identifies the class to which symmetries belong
            % CT = character table
            % irreps
            %
            % MATHEMATICS OF COMPUTATION, VOLUME 24, NUMBER 111, JULY, 1970
            % Computing Irreducible Representations of Groups
            % By John D. Dixon

            import am_lib.* am_dft.*

            % get regular rep G by putting identity along diagonal of multiplciation table
            [MT,~,I] = get_multiplication_table(S); nSs = size(MT,2);
            RR = double(accessc_(MT,I)==permute([1:nSs],[1,3,2]));
            
            % initialize decomposition loop
            [d,~,nSs] = size(RR); U = eye(nSs); inds = ones(nSs,1); ninds = 1;

            % loop until irreps are fully decomposed
            while true
                % loop over cycle structures
                for j = 1:max(inds)
                    ex_ = inds==j;
                    H = dixon_decomposition_( RR(ex_,ex_,:) );
                    [Vp,E] = eig(H,'vector'); [Vp] = orth_(Vp,E);
                    for ig = 1:nSs
                        RR(ex_,ex_,ig) = Vp\RR(ex_,ex_,ig)*Vp;
                    end
                    U(:,ex_) = (U(:,ex_)*Vp);
                end
                inds = [1:d]*merge_(double(sum(abs(RR),3)>am_lib.eps));
                if ninds == max(inds); break; else; ninds = max(inds); end
            end
            
            % get irreducible representations
            nirreps = max(inds); IR = cell(1,nirreps);
            for i = 1:nirreps
                IR{i} = RR(inds==i,inds==i,1:nSs);
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

        function [CT,cc_id,ir_id] = get_character_table(IR)
            % R=generate_pg(32,true);
            % IR = get_irreducible_representations(R);
            % [CT,cc_id,ir_id] = get_character_table(IR);
            % print_character_table(R,CT,cc_id)
            
            import am_lib.*
            
            %
            nirreps = numel(IR);
            nsyms   = size(IR{1},3);
            
            % get character table
            CT = zeros(nirreps,nsyms);
            for i = 1:nsyms; for j = 1:nirreps
                CT(j,i) = trace(IR{j}(:,:,i));
            end; end
        
            % correct numerical error
            CT = wdv_(CT);

            % get irreducible irreps
            [CT,~,cc_id]=uniquec_(CT);

            % get irreducible classes
            [CT,~,ir_id]=uniquec_(CT.'); CT=CT.';
        
        end
        
        function                print_character_table(R,CT,cc_id)
            % clear;clc
            % R=generate_pg(32,true);
            % RR = get_regular_representation(R);
            % IR = get_irreducible_representations(RR);
            % [CT,cc_id,ir_id] = get_character_table(IR);
            % print_character_table(R,CT,cc_id)

            import am_dft.* am_lib.*
            
            % identify the prototypical symmetries
            [~,unique_c_id]=unique(cc_id);

            % get number of classes
            nclasses=numel(unique_c_id);

            fprintf('      '); fprintf('%9s',repmat('---------',1,nclasses)); fprintf('\n');
            % class label 
            fprintf('      '); for i = 1:nclasses; fprintf('%9s',['#',num2str(i)]); end; fprintf('\n');
            fprintf('      '); fprintf('%9s',repmat(' --------',1,nclasses)); fprintf('\n');
            % print class name 
            ps_name = decode_ps(identify_point_symmetries(R(:,:,unique_c_id)));
            fprintf('      '); for i = 1:nclasses; fprintf('%9s',ps_name{i}); end; fprintf('\n');
            % print class elements
            cl_size = sum(cc_id.'==cc_id(unique_c_id),2);
            fprintf('      '); for i = 1:nclasses; fprintf('%9i',cl_size(i)); end; fprintf('\n');
            % bar
            fprintf('      '); fprintf('%9s',repmat(' --------',1,nclasses)); fprintf('\n');
            % print character table
            for l = 1:size(CT,1)
                chi = CT(l,:); chi = wdv_(chi);
                fprintf('      '); fprintf('%9.4g',chi); fprintf('\n');
            end
            % bar
            fprintf('      '); fprintf('%9s',repmat(' --------',1,nclasses)); fprintf('\n');
            % print SO(3) and SU(2) characters single and double groups
            character_ = @(l,alpha) sin((l+1/2)*alpha)./sin(alpha/2);
            alpha = 2*pi./get_order(R);
            for l = [0:.5:4]
                chi = character_(l,alpha(:,unique_c_id)); chi = wdv_(chi);
                fprintf('      '); fprintf('%9.4g',chi); fprintf('\n');
            end
            fprintf('      '); fprintf('%9s',repmat('---------',1,nclasses)); fprintf('\n');
            
            % symmetries in each class
            fprintf('Symmetries in each classes:\n')
            ps_name_long = get_long_ps_name(R);
            for i = 1:numel(unique_c_id)
                fprintf('%9s:%s\n',['#',num2str(i)], sprintf(' %-12s',ps_name_long{cc_id==i}));
            end
        end

        function cc_id        = identify_classes(MT)
            %
            % for AX = XB, if elements A and B are conjugate pairs for some other element X in the group,  they are in the same class
            %
            import am_lib.* am_dft.*

            nsyms = size(MT,2);
            % allocate space for conjugacy class
            cc_id = zeros(nsyms,1);
            % allocate space for conjugate elements
            conjugates = zeros(nsyms,1);
            % get inverse indicies
            I = (MT==1) * [1:size(MT,1)].';
            % determine conjugacy classes
            k = 0;
            for i = 1:nsyms
            if cc_id(i)==0
                k=k+1;
                % conjugate each element with all other group elements
                % A = X(j) * B * X(j)^-1
                for j = 1:nsyms
                    conjugates(j) = MT(j,MT(i,I(j)));
                end
                % for each subgroup element created by conjugation find the corresponding index of the element in the group
                % in order to save the class class_id number
                for j = 1:nsyms
                    cc_id( conjugates(j) ) = k;
                end
            end
            end
            % relabel classes based on how many elements each class has
            cc_id = reindex_using_occurances(cc_id);
        end

        function [ps_id,tr,dt]= identify_point_symmetries(R)
            nsyms=size(R,3); ps_id = zeros(1,nsyms); tr = zeros(1,nsyms); dt = zeros(1,nsyms);
            for i = 1:nsyms
                % get trace and determinant (fractional)
                tr(i) = trace(R(1:3,1:3,i)); dt(i) = det(R(1:3,1:3,i));
                if     (tr(i)==+3 && dt(i)==+1); ps_id(i) = 1;  % 'e'
                elseif (tr(i)==-1 && dt(i)==+1); ps_id(i) = 2;  % 'c_2'
                elseif (tr(i)==+0 && dt(i)==+1); ps_id(i) = 3;  % 'c_3'
                elseif (tr(i)==+1 && dt(i)==+1); ps_id(i) = 4;  % 'c_4'
                elseif (tr(i)==+2 && dt(i)==+1); ps_id(i) = 5;  % 'c_6'
                elseif (tr(i)==-3 && dt(i)==-1); ps_id(i) = 6;  % 'i'
                elseif (tr(i)==+1 && dt(i)==-1); ps_id(i) = 7;  % 's_2'
                elseif (tr(i)==+0 && dt(i)==-1); ps_id(i) = 8;  % 's_6'
                elseif (tr(i)==-1 && dt(i)==-1); ps_id(i) = 9;  % 's_4'
                elseif (tr(i)==-2 && dt(i)==-1); ps_id(i) = 10; % 's_3'
                else;                            ps_id(i) = 0;  % unknown
                end
            end
        end

        function bv_code      = identify_bravais_lattice(bas, tol,algo)
            % bv_code = identify_bravais_lattice(bas, tol,algo)
            %
            % Identifies the type of crystal system given a basis. 
            %
            % bas     : unit cell column basis vectors 
            % tol     : numerical tolerence
            % algo (1): de Graef p 86
            %      (2): Tinkham  p 61
            % 
            % Antonio Mei Sep 20 2017
            
            import am_dft.bas2abc am_lib.eq_ am_dft.identify_bravais_lattice
            
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
                    t  = numel(uniquetol(ii, tol));
                    o  = numel(uniquetol(ij, tol));
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
                    abc = bas2abc(bas);
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
                    elseif ~eq_(abc(1),abc(2), tol) && ~eq_(abc(1),abc(3), tol) && ~eq_(abc(2),abc(3), tol)
                            bv_code = 4; % tetragonal   (a=b!=c, alpha=beta=gamma!=90)
                    else
                            bv_code = 2; % hexagonal    (a=b!=c, alpha=beta=90,gamma=120)
                    end
            end
            
            % if an algo is not solicited, run both and check for match.
            if nargin < 3; if identify_bravais_lattice(bas, tol,1)~=identify_bravais_lattice(bas, tol,2)
                warning('Algorithm mismatch. \n')
            end; end
        end

        function pg_code      = identify_pointgroup(R)
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
            import am_lib.* am_dft.*
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
                switch ps_code(i)
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
            end
        end
        
        function ps_name_long = get_long_ps_name(R)
            import am_dft.* am_lib.*
            nsyms=size(R,3); 
            ps_id=identify_point_symmetries(R);
            ps_name=decode_ps(ps_id);
            ps_axis=round_(R_axis_(R));
            syms x y z
            for i = 1:nsyms
                ps_name_long{i} = sprintf('%s',strtrim(ps_name{i}));
                if ps_id(i)~=1 && ps_id(i)~=6
                    ax = ps_axis(:,i).'; if sum(lt_(ax,0))>sum(gt_(ax,0)); ax=-ax; end
                    ax_name = sprintf('%s',ax*[x;y;z]); ax_name=strrep(ax_name,' ',''); ax_name=strrep(ax_name,'+',''); 
                    ps_name_long{i} = sprintf('%s(%s)',ps_name_long{i}, ax_name );
                end
            end
        end
        
        function ss_name_long = get_long_ss_name(S)
            import am_dft.* am_lib.*
            ss_name_long = get_long_ps_name(S(1:3,1:3,:));
            nsyms=size(S,3); E = find(all(all(eq_(S,eye(4)),1),2));
            for i = 1:nsyms
                if i==E; continue; end
                t = S(1:3,4,i);
                t_name = sprintf('%s,',sym(t)); t_name=strrep(t_name,' ',''); t_name=t_name(1:(end-1)); 
                ss_name_long{i} = sprintf('%s [%s]',ss_name_long{i}, t_name );
            end
        end
            
        function order        = get_order(S, tol)

            import am_lib.* am_dft.*

            if nargin<2; tol = am_dft.tiny; end
            
            % define comparison function
            check_ = @(x) any(~eq_(x(:),0,tol));

            % get sizes
            s = size(S); nSs = s(3);

            if     s(1) == 3 && s(2) == 3       % point symmetry
                order = ones(1,nSs); 
                for i = 1:nSs
                    X=S(1:3,1:3,i);
                    while check_(X - eye(3))
                        order(i) = order(i) + 1;
                        X = X*S(:,:,i);
                    end
                end
            elseif s(1) == 4 && s(2) == 4       % space symmetry
                order = ones(1,nSs); 
                for i = 1:nSs
                    X=S(1:4,1:4,i);
                    while check_(X - eye(4))
                        order(i) = order(i) + 1;
                        X = X*S(:,:,i); X(1:3,4)=mod_(X(1:3,4));
                    end
                end
            end
        end

        function S            = generate_sg(sg_code,from_memory)

            import am_lib.* am_dft.*

            if nargin<2; from_memory=true; end

            if from_memory
                % ~ 500x faster for large space groups
                switch sg_code
                case 1; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1];
                case 2; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1];
                case 3; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1];
                case 4; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,0,1];
                case 5; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1];
                case 6; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1];
                case 7; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1];
                case 8; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1];
                case 9; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 10; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1];
                case 11; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,0,1];
                case 12; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1];
                case 13; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1];
                case 14; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1];
                case 15; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 16; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1];
                case 17; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1];
                case 18; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1];
                case 19; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1];
                case 20; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1];
                case 21; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1];
                case 22; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1];
                case 23; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1];
                case 24; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1];
                case 25; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1];
                case 26; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1];
                case 27; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1];
                case 28; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,0,1];
                case 29; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1];
                case 30; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1];
                case 31; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1];
                case 32; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1];
                case 33; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 34; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 35; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1];
                case 36; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1];
                case 37; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 38; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1];
                case 39; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1];
                case 40; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 41; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1];
                case 42; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1];
                case 43; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/4,1/4,1/4,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/4,3/4,3/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,3/4,1/4,3/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,3/4,3/4,1/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/4,1/4,1/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,3/4,3/4,1/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,3/4,1/4,3/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/4,3/4,3/4,1];
                case 44; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 45; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1];
                case 46; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1];
                case 47; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1];
                case 48; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 49; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1];
                case 50; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1];
                case 51; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,0,1];
                case 52; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1];
                case 53; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1];
                case 54; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1];
                case 55; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1];
                case 56; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1];
                case 57; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,0,1];
                case 58; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 59; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1];
                case 60; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1];
                case 61; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1];
                case 62; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 63; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1];
                case 64; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1];
                case 65; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1];
                case 66; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 67; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1];
                case 68; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1];
                case 69; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1];
                case 70; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/4,1/4,1/4,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/4,3/4,3/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,3/4,1/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,3/4,3/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,3/4,1/4,3/4,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/4,1/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/4,1/4,1/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,3/4,3/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,3/4,1/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,3/4,3/4,1/4,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/4,3/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/4,3/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/4,3/4,3/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,3/4,3/4,1/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,3/4,1/4,3/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/4,1/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1];
                case 71; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 72; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1];
                case 73; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1];
                case 74; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 75; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1];
                case 76; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,3/4,1];
                case 77; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1];
                case 78; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,3/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/4,1];
                case 79; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 80; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,1/2,1/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,0,3/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,0,3/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,1/2,1/4,1];
                case 81; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1];
                case 82; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1];
                case 83; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1];
                case 84; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1];
                case 85; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1];
                case 86; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1];
                case 87; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1];
                case 88; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,1/2,1/4,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,0,3/4,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,3/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,0,3/4,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,3/4,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,1/2,1/4,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/4,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1];
                case 89; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1];
                case 90; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1];
                case 91; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/4,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/4,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,3/4,1];
                case 92; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/4,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,3/4,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1];
                case 93; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1];
                case 94; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1];
                case 95; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,3/4,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,3/4,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,1/4,1];
                case 96; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,3/4,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,3/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/4,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1];
                case 97; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1];
                case 98; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,1/2,1/4,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,3/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,0,3/4,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,0,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/4,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,1/2,1/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,3/4,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1];
                case 99; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1];
                case 100; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,0,1];
                case 101; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1];
                case 102; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1];
                case 103; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1];
                case 104; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 105; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1];
                case 106; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 107; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 108; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,0,1];
                case 109; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,1/2,1/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,0,3/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,0,3/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,1/2,1/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,1/2,1/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,1/2,1/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,0,3/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,0,3/4,1];
                case 110; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,1/2,1/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,0,3/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,0,3/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,1/2,3/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,1/2,3/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,1/2,1/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,0,1/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,0,1/4,1];
                case 111; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1];
                case 112; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1];
                case 113; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,0,1];
                case 114; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 115; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1];
                case 116; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1];
                case 117; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,0,1];
                case 118; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1];
                case 119; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1];
                case 120; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,0,1];
                case 121; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 122; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,3/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/4,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,3/4,1,0,1,0,0,1,0,0,0,0,0,1,0,0,1/2,1/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,0,3/4,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/4,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,0,3/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,1/2,1/4,1];
                case 123; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1];
                case 124; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1];
                case 125; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,0,1];
                case 126; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 127; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,0,1];
                case 128; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 129; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,0,1];
                case 130; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 131; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1];
                case 132; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1];
                case 133; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1];
                case 134; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1];
                case 135; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 136; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1];
                case 137; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 138; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,0,1];
                case 139; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 140; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,0,1];
                case 141; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,1/2,1/4,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,3/4,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,0,3/4,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/4,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,3/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,0,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/4,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,3/4,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,1/2,1/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,3/4,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/4,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,0,3/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,1/2,1/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,0,3/4,1,0,1,0,0,1,0,0,0,0,0,1,0,0,1/2,1/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1];
                case 142; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,1/2,1/4,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/4,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,0,3/4,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,3/4,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,3/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,0,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,3/4,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,3/4,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,1/2,1/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/4,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/4,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,0,1/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,1/2,3/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,0,1/4,1,0,1,0,0,1,0,0,0,0,0,1,0,0,1/2,3/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1];
                case 143; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1];
                case 144; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,1/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,2/3,1];
                case 145; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,2/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/3,1];
                case 146; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,2/3,1/3,1/3,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/3,2/3,2/3,1,0,1,0,0,-1,-1,0,0,0,0,1,0,2/3,1/3,1/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,1/3,2/3,2/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,2/3,1/3,1/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,1/3,2/3,2/3,1];
                case 147; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1];
                case 148; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,2/3,1/3,1/3,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/3,2/3,2/3,1,0,-1,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,2/3,1/3,1/3,1,1,0,0,0,0,1,0,0,0,0,1,0,1/3,2/3,2/3,1,0,1,0,0,-1,-1,0,0,0,0,1,0,2/3,1/3,1/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,-1,0,1/3,2/3,2/3,1,1,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,-1,0,2/3,1/3,1/3,1,0,1,0,0,-1,-1,0,0,0,0,1,0,1/3,2/3,2/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,2/3,1/3,1/3,1,1,1,0,0,-1,0,0,0,0,0,-1,0,1/3,2/3,2/3,1,1,1,0,0,-1,0,0,0,0,0,-1,0,2/3,1/3,1/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,1/3,2/3,2/3,1];
                case 149; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1];
                case 150; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1];
                case 151; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,1/3,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,2/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,2/3,1,1,1,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,-1,0,0,0,1/3,1];
                case 152; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,1/3,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,2/3,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,1/3,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,2/3,1];
                case 153; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,2/3,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/3,1,1,1,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,-1,0,0,0,2/3,1];
                case 154; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,2/3,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/3,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,2/3,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,1/3,1];
                case 155; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,2/3,1/3,1/3,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/3,2/3,2/3,1,0,1,0,0,-1,-1,0,0,0,0,1,0,2/3,1/3,1/3,1,0,1,0,0,1,0,0,0,0,0,-1,0,2/3,1/3,1/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/3,2/3,2/3,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,1/3,2/3,2/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,2/3,1/3,1/3,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,2/3,1/3,1/3,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,2/3,1/3,1/3,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,1/3,2/3,2/3,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,1/3,2/3,2/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,1/3,2/3,2/3,1];
                case 156; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1];
                case 157; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1];
                case 158; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,1/2,1];
                case 159; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,-1,-1,0,0,0,0,1,0,0,0,1/2,1];
                case 160; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,2/3,1/3,1/3,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/3,2/3,2/3,1,0,1,0,0,-1,-1,0,0,0,0,1,0,2/3,1/3,1/3,1,0,-1,0,0,-1,0,0,0,0,0,1,0,2/3,1/3,1/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,1/3,2/3,2/3,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/3,2/3,2/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,2/3,1/3,1/3,1,1,1,0,0,0,-1,0,0,0,0,1,0,2/3,1/3,1/3,1,-1,0,0,0,1,1,0,0,0,0,1,0,2/3,1/3,1/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,1/3,2/3,2/3,1,1,1,0,0,0,-1,0,0,0,0,1,0,1/3,2/3,2/3,1,-1,0,0,0,1,1,0,0,0,0,1,0,1/3,2/3,2/3,1];
                case 161; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,2/3,1/3,1/3,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/3,2/3,2/3,1,0,1,0,0,-1,-1,0,0,0,0,1,0,2/3,1/3,1/3,1,0,-1,0,0,-1,0,0,0,0,0,1,0,2/3,1/3,5/6,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,-1,-1,0,0,0,0,1,0,1/3,2/3,2/3,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/3,2/3,1/6,1,-1,-1,0,0,1,0,0,0,0,0,1,0,2/3,1/3,1/3,1,1,1,0,0,0,-1,0,0,0,0,1,0,2/3,1/3,5/6,1,-1,0,0,0,1,1,0,0,0,0,1,0,2/3,1/3,5/6,1,-1,-1,0,0,1,0,0,0,0,0,1,0,1/3,2/3,2/3,1,1,1,0,0,0,-1,0,0,0,0,1,0,1/3,2/3,1/6,1,-1,0,0,0,1,1,0,0,0,0,1,0,1/3,2/3,1/6,1];
                case 162; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1];
                case 163; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,-1,0,0,0,1,1,0,0,0,0,-1,0,0,0,1/2,1,1,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,-1,-1,0,0,0,0,1,0,0,0,1/2,1];
                case 164; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1];
                case 165; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,1/2,1,1,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,1/2,1];
                case 166; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,2/3,1/3,1/3,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/3,2/3,2/3,1,0,-1,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,2/3,1/3,1/3,1,1,0,0,0,0,1,0,0,0,0,1,0,1/3,2/3,2/3,1,0,1,0,0,-1,-1,0,0,0,0,1,0,2/3,1/3,1/3,1,0,1,0,0,1,0,0,0,0,0,-1,0,2/3,1/3,1/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/3,2/3,2/3,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,-1,0,1/3,2/3,2/3,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/3,2/3,2/3,1,1,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,2/3,1/3,1/3,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,-1,0,2/3,1/3,1/3,1,0,1,0,0,-1,-1,0,0,0,0,1,0,1/3,2/3,2/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,2/3,1/3,1/3,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,2/3,1/3,1/3,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,2/3,1/3,1/3,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,1/3,2/3,2/3,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,1/3,2/3,2/3,1,1,1,0,0,-1,0,0,0,0,0,-1,0,1/3,2/3,2/3,1,1,1,0,0,0,-1,0,0,0,0,1,0,1/3,2/3,2/3,1,-1,0,0,0,1,1,0,0,0,0,1,0,1/3,2/3,2/3,1,1,1,0,0,0,-1,0,0,0,0,1,0,2/3,1/3,1/3,1,-1,0,0,0,1,1,0,0,0,0,1,0,2/3,1/3,1/3,1,1,1,0,0,-1,0,0,0,0,0,-1,0,2/3,1/3,1/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,1/3,2/3,2/3,1];
                case 167; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,2/3,1/3,1/3,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/3,2/3,2/3,1,0,-1,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,2/3,1/3,1/3,1,1,0,0,0,0,1,0,0,0,0,1,0,1/3,2/3,2/3,1,0,1,0,0,-1,-1,0,0,0,0,1,0,2/3,1/3,1/3,1,0,1,0,0,1,0,0,0,0,0,-1,0,2/3,1/3,5/6,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/3,2/3,1/6,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,1,1,0,0,0,0,-1,0,1/3,2/3,2/3,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/3,2/3,1/6,1,1,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,2/3,1/3,5/6,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,1,1,0,0,0,0,-1,0,2/3,1/3,1/3,1,0,1,0,0,-1,-1,0,0,0,0,1,0,1/3,2/3,2/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,2/3,1/3,1/3,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,2/3,1/3,5/6,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,2/3,1/3,5/6,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,1/3,2/3,1/6,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,1/3,2/3,1/6,1,1,1,0,0,-1,0,0,0,0,0,-1,0,1/3,2/3,2/3,1,1,1,0,0,0,-1,0,0,0,0,1,0,1/3,2/3,1/6,1,-1,0,0,0,1,1,0,0,0,0,1,0,1/3,2/3,1/6,1,1,1,0,0,0,-1,0,0,0,0,1,0,2/3,1/3,5/6,1,-1,0,0,0,1,1,0,0,0,0,1,0,2/3,1/3,5/6,1,1,1,0,0,-1,0,0,0,0,0,-1,0,2/3,1/3,1/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,1/3,2/3,2/3,1];
                case 168; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1];
                case 169; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,1/3,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,2/3,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,5/6,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/6,1];
                case 170; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,2/3,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/3,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,1/6,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,5/6,1];
                case 171; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,2/3,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/3,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,2/3,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/3,1];
                case 172; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,1/3,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,2/3,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,1/3,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,2/3,1];
                case 173; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,1/2,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1];
                case 174; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1];
                case 175; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1];
                case 176; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,1/2,1,1,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,-1,0,0,0,1/2,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1];
                case 177; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1];
                case 178; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,1/3,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,1/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,2/3,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,5/6,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,2/3,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,5/6,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/6,1,1,1,0,0,0,-1,0,0,0,0,-1,0,0,0,1/6,1,-1,0,0,0,1,1,0,0,0,0,-1,0,0,0,1/2,1];
                case 179; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,2/3,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,2/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/3,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,1/6,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,1/3,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/6,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,5/6,1,1,1,0,0,0,-1,0,0,0,0,-1,0,0,0,5/6,1,-1,0,0,0,1,1,0,0,0,0,-1,0,0,0,1/2,1];
                case 180; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,2/3,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,2/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/3,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,2/3,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,1/3,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,2/3,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/3,1,1,1,0,0,0,-1,0,0,0,0,-1,0,0,0,1/3,1,-1,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1];
                case 181; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,1/3,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,1/3,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,2/3,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,1/3,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,2/3,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/3,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,2/3,1,1,1,0,0,0,-1,0,0,0,0,-1,0,0,0,2/3,1,-1,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1];
                case 182; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,1,1,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,-1,0,0,0,1,1,0,0,0,0,-1,0,0,0,1/2,1];
                case 183; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1];
                case 184; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,1/2,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,-1,-1,0,0,0,0,1,0,0,0,1/2,1];
                case 185; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,1/2,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,1/2,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1];
                case 186; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,1/2,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,-1,-1,0,0,0,0,1,0,0,0,1/2,1];
                case 187; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1];
                case 188; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,-1,0,0,0,1/2,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,1,1,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1];
                case 189; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1];
                case 190; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,-1,0,0,0,1/2,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,-1,-1,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,-1,-1,0,0,0,0,1,0,0,0,1/2,1];
                case 191; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1];
                case 192; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,1/2,1,1,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,1/2,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,-1,0,0,0,1,1,0,0,0,0,-1,0,0,0,1/2,1,-1,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,-1,-1,0,0,0,0,1,0,0,0,1/2,1];
                case 193; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,1/2,1,1,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,-1,0,0,0,1/2,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,1/2,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,1,1,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,-1,-1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1];
                case 194; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,1,1,0,0,0,0,1,0,0,0,1/2,1,-1,-1,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,-1,-1,0,0,0,0,-1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,-1,-1,0,0,0,0,-1,0,0,0,1/2,1,1,1,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,1,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,1,1,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,-1,0,0,0,1,1,0,0,0,0,-1,0,0,0,1/2,1,-1,-1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,-1,-1,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,-1,-1,0,0,0,0,1,0,0,0,1/2,1];
                case 195; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1];
                case 196; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,0,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1];
                case 197; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1];
                case 198; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1];
                case 199; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1/2,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,0,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,1/2,1];
                case 200; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,1,0,0,0,0,0,1];
                case 201; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/2,1/2,1/2,1];
                case 202; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,0,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,0,1,0,0,0,0,-1,0,1,0,0,0,0,1/2,1/2,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,-1,0,0,0,0,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/2,0,1/2,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/2,1/2,0,1];
                case 203; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/4,1/4,1/4,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/4,3/4,3/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,3/4,1/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,3/4,3/4,1/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,3/4,1/4,3/4,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/4,1/4,1/4,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/4,1/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/4,1/4,1/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,3/4,3/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,3/4,1/4,3/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,3/4,3/4,1/4,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/4,3/4,3/4,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/4,3/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/4,3/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,0,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,3/4,1/4,3/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/4,3/4,3/4,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,1,0,0,0,0,-1,0,1,0,0,0,3/4,3/4,1/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,3/4,3/4,1/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,3/4,1/4,3/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,3/4,1/4,3/4,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,3/4,3/4,1/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/4,3/4,3/4,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/4,1/4,1/4,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/4,1/4,1/4,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/4,1/4,1/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/4,1/4,1/4,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/4,1/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,1,0,0,0,0,-1,0,1,0,0,0,3/4,1/4,3/4,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,-1,0,0,0,0,1,0,1,0,0,0,3/4,3/4,1/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/4,3/4,3/4,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/4,3/4,3/4,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/4,3/4,3/4,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,1,0,0,0,0,1,0,-1,0,0,0,3/4,3/4,1/4,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,3/4,1/4,3/4,1,0,1,0,0,0,0,1,0,-1,0,0,0,3/4,1/4,3/4,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,1,0,0,3/4,3/4,1/4,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,3/4,1/4,3/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,3/4,3/4,1/4,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/4,3/4,3/4,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/4,1/4,1/4,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/4,1/4,1/4,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/4,1/4,1/4,1,0,0,-1,0,1,0,0,0,0,1,0,0,3/4,1/4,3/4,1,0,0,1,0,1,0,0,0,0,-1,0,0,3/4,3/4,1/4,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/4,3/4,3/4,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/4,3/4,3/4,1,0,0,1,0,-1,0,0,0,0,1,0,0,3/4,3/4,1/4,1,0,0,1,0,-1,0,0,0,0,1,0,0,3/4,1/4,3/4,1];
                case 204; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/2,1/2,1/2,1];
                case 205; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/2,1/2,0,1];
                case 206; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1/2,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,0,0,-1,0,1,0,0,0,0,1/2,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/2,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,0,0,1/2,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,0,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,0,1/2,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/2,0,0,1,0,0,1,0,-1,0,0,0,0,1,0,0,0,0,1/2,1];
                case 207; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1];
                case 208; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,1/2,1];
                case 209; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,1,0,0,1,0,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,1,0,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,0,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,0,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,0,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/2,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,0,1/2,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,0,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,0,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,0,1/2,1/2,1,1,0,0,0,0,0,1,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/2,0,1/2,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/2,0,1/2,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/2,0,1/2,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/2,1/2,0,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/2,1/2,0,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,0,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,0,1];
                case 210; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,3/4,1/4,3/4,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,3/4,3/4,1/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/4,1/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/4,1/4,1/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,3/4,3/4,1/4,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,0,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,3/4,3/4,1/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/4,3/4,3/4,1,1,0,0,0,0,0,-1,0,0,1,0,0,3/4,1/4,3/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/4,3/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/4,3/4,3/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,3/4,1/4,3/4,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,1,0,0,1,0,0,-1,0,0,0,3/4,1/4,3/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/4,1/4,1/4,1,1,0,0,0,0,0,-1,0,0,1,0,0,3/4,3/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,3/4,1/4,3/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/4,3/4,3/4,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/4,3/4,3/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,3/4,3/4,1/4,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/4,1/4,1/4,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/4,3/4,3/4,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/4,1/4,1/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,3/4,1/4,3/4,1,-1,0,0,0,0,0,1,0,0,1,0,0,3/4,3/4,1/4,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/4,1/4,1/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/4,1/4,1/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/4,3/4,3/4,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/4,3/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,3/4,3/4,1/4,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/4,1/4,1/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,3/4,1/4,3/4,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/4,1/4,1/4,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/4,3/4,3/4,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,3/4,3/4,1/4,1,-1,0,0,0,0,0,1,0,0,1,0,0,3/4,1/4,3/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/4,3/4,3/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/4,1/4,1/4,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,3/4,3/4,1/4,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,3/4,1/4,3/4,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/4,1/4,1/4,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/4,3/4,3/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,3/4,1/4,3/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,3/4,3/4,1/4,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,3/4,1/4,3/4,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,3/4,3/4,1/4,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/4,3/4,3/4,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/4,1/4,1/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,3/4,3/4,1/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,3/4,1/4,3/4,1];
                case 211; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,1/2,1];
                case 212; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/4,3/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/4,1/4,1/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,3/4,1/4,3/4,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,3/4,1/4,3/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,3/4,3/4,1/4,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/4,3/4,3/4,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,3/4,3/4,1/4,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/4,1/4,1/4,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/4,3/4,3/4,1,-1,0,0,0,0,0,1,0,0,1,0,0,3/4,1/4,3/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/4,1/4,1/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,3/4,3/4,1/4,1];
                case 213; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,3/4,1/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,3/4,3/4,3/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/4,3/4,1/4,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/4,3/4,1/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/4,1/4,3/4,1,1,0,0,0,0,0,-1,0,0,1,0,0,3/4,1/4,1/4,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/4,1/4,3/4,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,3/4,3/4,3/4,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,3/4,1/4,1/4,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/4,3/4,1/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,3/4,3/4,3/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/4,1/4,3/4,1];
                case 214; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,3/4,1/4,1/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/4,3/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,3/4,3/4,3/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/4,3/4,1/4,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/4,3/4,1/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/4,1/4,3/4,1,1,0,0,0,0,0,-1,0,0,1,0,0,3/4,1/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/4,1/4,1/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,3/4,1/4,3/4,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,0,1,0,0,-1,0,0,0,3/4,1/4,3/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,3/4,3/4,1/4,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/4,3/4,3/4,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/4,1/4,3/4,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,3/4,3/4,3/4,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,3/4,1/4,1/4,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/4,3/4,1/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,3/4,3/4,3/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/4,1/4,3/4,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,1/2,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,3/4,3/4,1/4,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/4,1/4,1/4,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,0,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/4,3/4,3/4,1,-1,0,0,0,0,0,1,0,0,1,0,0,3/4,1/4,3/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,1/2,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/4,1/4,1/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,3/4,3/4,1/4,1];
                case 215; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,0,0,0,1,-1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,-1,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,1];
                case 216; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,1,0,0,1,0,0,1,0,0,0,0,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,1/2,1/2,1,1,0,0,0,0,0,1,0,0,1,0,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,0,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,0,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,0,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,0,0,0,1,-1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,-1,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,0,-1,0,0,-1,0,0,0,0,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,-1,0,0,0,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,0,1,0,0,-1,0,0,0,0,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,-1,0,0,-1,0,0,1,0,0,0,0,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/2,0,1/2,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/2,0,1/2,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/2,0,1/2,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,0,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/2,1/2,0,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/2,1/2,0,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,0,1];
                case 217; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,0,0,0,1,-1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,-1,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,1/2,1];
                case 218; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,1/2,1];
                case 219; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,1/2,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,1/2,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,0,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,1,0,0,1,0,0,1,0,0,0,0,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,1/2,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,1/2,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/2,0,0,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/2,0,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/2,0,0,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/2,0,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/2,0,0,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,0,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,0,-1,0,0,-1,0,0,0,0,1/2,0,1,-1,0,0,0,0,0,1,0,0,-1,0,0,0,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,0,1,0,0,-1,0,0,0,0,1/2,0,1,-1,0,0,0,0,0,-1,0,0,1,0,0,0,1/2,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,0,-1,0,0,1,0,0,0,0,1/2,0,1,1,0,0,0,0,0,-1,0,0,-1,0,0,0,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,0,0,1/2,1,-1,0,0,0,0,0,1,0,0,-1,0,0,0,0,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,0,1,0,0,-1,0,0,0,0,0,1/2,1,-1,0,0,0,0,0,-1,0,0,1,0,0,0,0,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,0,-1,0,0,1,0,0,0,0,0,1/2,1,1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,1/2,1];
                case 220; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/4,1/4,1/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,3/4,3/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/4,3/4,3/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,3/4,3/4,1/4,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1/4,1/4,1/4,1,0,-1,0,0,1,0,0,0,0,0,-1,0,3/4,1/4,3/4,1,1,0,0,0,0,0,1,0,0,1,0,0,1/4,1/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,3/4,1/4,1/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/4,1/4,3/4,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,0,1,0,0,1,0,0,0,3/4,3/4,3/4,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/4,3/4,1/4,1,1,0,0,0,0,0,1,0,0,1,0,0,3/4,3/4,3/4,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/4,3/4,3/4,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/4,3/4,3/4,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,0,1,0,0,-1,0,0,0,3/4,3/4,1/4,1,-1,0,0,0,0,0,-1,0,0,1,0,0,3/4,3/4,1/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,0,-1,0,0,1,0,0,0,3/4,1/4,3/4,1,1,0,0,0,0,0,-1,0,0,-1,0,0,3/4,1/4,3/4,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,1/2,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,3/4,1/4,1/4,1,-1,0,0,0,0,0,1,0,0,-1,0,0,3/4,1/4,1/4,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,0,0,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/4,1/4,3/4,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/4,1/4,3/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,1/2,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/4,3/4,1/4,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/4,3/4,1/4,1];
                case 221; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,1];
                case 222; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,0,0,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,0,0,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,1/2,1];
                case 223; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,1/2,1];
                case 224; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,0,-1,0,0,1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,-1,0,0,1,0,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,1,0,0,-1,0,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,1];
                case 225; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,1,0,0,1,0,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,1,0,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,0,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,0,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,0,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/2,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,0,1/2,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,0,1,0,0,0,0,-1,0,1,0,0,0,0,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,1/2,1/2,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,-1,0,0,0,0,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,-1,0,0,-1,0,0,1,0,0,0,0,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,-1,0,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/2,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/2,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,0,1/2,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/2,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,0,1/2,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/2,0,1/2,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,0,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,0,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,0,1/2,1/2,1,1,0,0,0,0,0,1,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/2,0,1/2,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/2,0,1/2,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/2,0,1/2,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/2,1/2,0,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/2,1/2,0,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,0,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,0,1,0,0,-1,0,0,0,0,1/2,1/2,1,1,0,0,0,0,0,1,0,0,1,0,0,0,1/2,1/2,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,0,-1,0,0,-1,0,0,0,0,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,1,0,0,1,0,0,1,0,0,0,0,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/2,0,1/2,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,0,1/2,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/2,0,1/2,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,0,1/2,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/2,1/2,0,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,0,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,1/2,0,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,0,1];
                case 226; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/2,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,0,0,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/2,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,1/2,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,1/2,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,0,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,1,0,0,1,0,0,-1,0,0,0,0,1/2,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,1/2,0,1,1,0,0,0,0,0,-1,0,0,1,0,0,0,1/2,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,0,0,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,1,0,0,0,0,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,0,1,0,0,0,0,-1,0,1,0,0,0,0,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,0,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,0,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/2,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,0,0,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/2,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/2,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,0,1/2,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/2,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,1/2,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,0,-1,0,0,1,0,0,0,0,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,1/2,0,1,-1,0,0,0,0,0,1,0,0,-1,0,0,0,1/2,0,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,1/2,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,0,-1,0,0,1,0,0,0,0,0,1/2,1,0,0,1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,-1,0,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/2,0,0,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,0,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/2,0,0,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/2,0,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/2,0,0,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/2,0,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,0,1/2,0,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,0,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,0,1/2,0,1,-1,0,0,0,0,0,1,0,0,1,0,0,0,1/2,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,0,1/2,0,1,1,0,0,0,0,0,1,0,0,-1,0,0,0,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,0,0,1/2,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,0,0,1/2,1,-1,0,0,0,0,0,1,0,0,1,0,0,0,0,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,0,0,1/2,1,1,0,0,0,0,0,1,0,0,-1,0,0,0,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/2,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,0,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/2,0,0,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,0,0,1,0,0,1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,0,0,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/2,0,0,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,0,1,0,0,-1,0,0,0,0,1/2,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,1/2,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,0,-1,0,0,-1,0,0,0,0,1/2,0,1,1,0,0,0,0,0,-1,0,0,-1,0,0,0,1/2,0,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,1,0,0,1,0,0,1,0,0,0,0,1/2,0,1,-1,0,0,0,0,0,-1,0,0,1,0,0,0,1/2,0,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,0,1,0,0,-1,0,0,0,0,0,1/2,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,1/2,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,0,0,1/2,1,1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,1/2,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,1/2,1,-1,0,0,0,0,0,-1,0,0,1,0,0,0,0,1/2,1];
                case 227; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,3/4,1/4,3/4,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/4,1/4,1/4,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,3/4,3/4,1/4,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/4,3/4,3/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/4,1/4,1/4,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,3/4,1/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/4,1/4,1/4,1,1,0,0,0,0,1,0,0,0,0,-1,0,3/4,1/4,3/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,3/4,3/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/4,3/4,3/4,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,0,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,3/4,3/4,1/4,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/4,1/4,1/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/4,3/4,3/4,1,1,0,0,0,0,0,-1,0,0,1,0,0,3/4,1/4,3/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/4,3/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,3/4,3/4,1/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/4,3/4,3/4,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,3/4,3/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/4,3/4,3/4,1,1,0,0,0,0,1,0,0,0,0,-1,0,3/4,3/4,1/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,3/4,1/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/4,1/4,1/4,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,1,0,0,1,0,0,-1,0,0,0,3/4,1/4,3/4,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/4,3/4,3/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/4,1/4,1/4,1,1,0,0,0,0,0,-1,0,0,1,0,0,3/4,3/4,1/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/4,1/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,3/4,1/4,3/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,3/4,1/4,3/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/4,3/4,3/4,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/4,3/4,3/4,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,3/4,1/4,3/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,3/4,3/4,1/4,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/4,1/4,1/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,3/4,3/4,1/4,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/4,3/4,3/4,1,0,1,0,0,0,0,-1,0,1,0,0,0,3/4,1/4,3/4,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/4,1/4,1/4,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/4,3/4,3/4,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,3/4,1/4,3/4,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/4,3/4,3/4,1,-1,0,0,0,0,0,1,0,0,1,0,0,3/4,3/4,1/4,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/4,1/4,1/4,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/4,1/4,1/4,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,3/4,3/4,1/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/4,1/4,1/4,1,0,1,0,0,0,0,1,0,-1,0,0,0,3/4,3/4,1/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/4,3/4,3/4,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/4,1/4,1/4,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/2,0,1/2,1,0,1,0,0,0,0,1,0,-1,0,0,0,3/4,1/4,3/4,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/4,3/4,3/4,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/2,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,0,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/4,3/4,3/4,1,-1,0,0,0,0,0,1,0,0,-1,0,0,0,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,3/4,1/4,3/4,1,0,-1,0,0,0,0,1,0,1,0,0,0,3/4,3/4,1/4,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,3/4,3/4,1/4,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/4,1/4,1/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,3/4,1/4,3/4,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/4,1/4,1/4,1,0,1,0,0,0,0,-1,0,1,0,0,0,3/4,3/4,1/4,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/4,3/4,3/4,1,0,1,0,0,1,0,0,0,0,0,1,0,0,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,0,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,3/4,3/4,1/4,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/4,1/4,1/4,1,-1,0,0,0,0,0,1,0,0,1,0,0,3/4,1/4,3/4,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/4,3/4,3/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/4,1/4,1/4,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/4,3/4,3/4,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/4,1/4,1/4,1,-1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,3/4,1/4,3/4,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,3/4,3/4,1/4,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,3/4,1/4,3/4,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/4,1/4,1/4,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/4,3/4,3/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,3/4,1/4,3/4,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/4,3/4,3/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,3/4,3/4,1/4,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,3/4,1/4,3/4,1,0,0,-1,0,0,-1,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/4,1/4,1/4,1,0,0,-1,0,0,-1,0,0,1,0,0,0,0,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,3/4,1/4,3/4,1,0,0,-1,0,1,0,0,0,0,1,0,0,3/4,1/4,3/4,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/2,0,1/2,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,3/4,3/4,1/4,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/4,3/4,3/4,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/4,3/4,3/4,1,0,0,1,0,0,-1,0,0,-1,0,0,0,0,1/2,1/2,1,0,0,1,0,0,-1,0,0,-1,0,0,0,0,0,0,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/4,1/4,1/4,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,0,1/2,1,1,0,0,0,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,3/4,3/4,1/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,3/4,1/4,3/4,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,3/4,3/4,1/4,1,0,0,1,0,-1,0,0,0,0,1,0,0,3/4,3/4,1/4,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,3/4,3/4,1/4,1,0,0,1,0,-1,0,0,0,0,1,0,0,3/4,1/4,3/4,1,-1,0,0,0,0,0,-1,0,0,1,0,0,0,1/2,1/2,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/4,3/4,3/4,1,0,0,-1,0,0,1,0,0,-1,0,0,0,0,1/2,1/2,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/2,0,1/2,1,1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/4,3/4,3/4,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/2,0,1/2,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,1,0,0,0,0,1,0,0,3/4,3/4,1/4,1,1,0,0,0,0,0,1,0,0,1,0,0,0,1/2,1/2,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,0,1/2,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/4,1/4,1/4,1,0,0,1,0,1,0,0,0,0,-1,0,0,3/4,1/4,3/4,1,-1,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/4,1/4,1/4,1,0,0,-1,0,0,1,0,0,-1,0,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/4,1/4,1/4,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,0,1/2,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,1/2,0,1];
                case 228; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,3/4,1/4,3/4,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,3/4,3/4,3/4,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,3/4,3/4,1/4,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,3/4,1/4,1/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,0,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/4,1/4,1/4,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/4,3/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/4,1/4,1/4,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/4,3/4,1/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,3/4,3/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,3/4,1/4,1/4,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,0,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,3/4,3/4,1/4,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,3/4,3/4,3/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/4,3/4,3/4,1,1,0,0,0,0,0,-1,0,0,1,0,0,3/4,1/4,3/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,3/4,1/4,1/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/4,1/4,3/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,1/2,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/4,3/4,3/4,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/4,1/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/4,3/4,3/4,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/4,1/4,3/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,3/4,1/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,3/4,3/4,3/4,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,1,0,0,1,0,0,-1,0,0,0,3/4,1/4,3/4,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,3/4,1/4,1/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/4,1/4,1/4,1,1,0,0,0,0,0,-1,0,0,1,0,0,3/4,3/4,1/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,3/4,3/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/4,3/4,1/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,3/4,1/4,3/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/4,3/4,3/4,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/4,3/4,3/4,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/4,3/4,1/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,3/4,3/4,1/4,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/4,1/4,1/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/4,1/4,3/4,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/4,3/4,3/4,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/4,3/4,1/4,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/4,1/4,1/4,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,3/4,1/4,1/4,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,3/4,1/4,3/4,1,0,-1,0,0,0,0,1,0,1,0,0,0,3/4,1/4,1/4,1,-1,0,0,0,0,0,1,0,0,1,0,0,3/4,3/4,1/4,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,3/4,3/4,3/4,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/4,1/4,1/4,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/4,1/4,3/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/4,1/4,1/4,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/4,1/4,3/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/4,3/4,3/4,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,3/4,3/4,3/4,1,0,0,-1,0,0,-1,0,0,1,0,0,0,0,1/2,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/4,3/4,1/4,1,0,1,0,0,0,0,-1,0,1,0,0,0,3/4,1/4,1/4,1,0,0,-1,0,0,-1,0,0,1,0,0,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,1/2,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/4,3/4,3/4,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/2,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/4,3/4,1/4,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/4,1/4,3/4,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,0,0,1,-1,0,0,0,0,0,1,0,0,-1,0,0,0,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,3/4,3/4,1/4,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/4,1/4,1/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,3/4,1/4,3/4,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/4,1/4,1/4,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/4,1/4,3/4,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/4,3/4,3/4,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,1/2,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,3/4,3/4,1/4,1,0,-1,0,0,0,0,1,0,1,0,0,0,3/4,3/4,3/4,1,-1,0,0,0,0,0,1,0,0,1,0,0,3/4,1/4,3/4,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/4,3/4,3/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/4,1/4,1/4,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,3/4,1/4,1/4,1,0,1,0,0,0,0,-1,0,1,0,0,0,3/4,3/4,3/4,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/4,3/4,1/4,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,-1,0,0,0,0,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,3/4,3/4,1/4,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,3/4,1/4,3/4,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/4,1/4,1/4,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/4,3/4,3/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,3/4,1/4,3/4,1,0,1,0,0,0,0,1,0,-1,0,0,0,3/4,1/4,1/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,3/4,3/4,1/4,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/4,3/4,1/4,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,0,0,1,0,-1,0,0,0,3/4,3/4,3/4,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/2,0,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,3/4,1/4,3/4,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/4,3/4,1/4,1,0,0,-1,0,0,1,0,0,-1,0,0,0,0,0,1/2,1,0,0,-1,0,0,1,0,0,-1,0,0,0,0,1/2,0,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,3/4,3/4,1/4,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/4,3/4,3/4,1,0,0,1,0,1,0,0,0,0,-1,0,0,3/4,1/4,1/4,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/2,0,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/4,1/4,1/4,1,1,0,0,0,0,0,-1,0,0,-1,0,0,0,1/2,0,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,0,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,3/4,3/4,1/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,3/4,1/4,3/4,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/4,1/4,3/4,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/4,1/4,3/4,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,1/2,1/2,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,0,0,1,-1,0,0,0,0,0,-1,0,0,1,0,0,0,0,1/2,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/4,1/4,3/4,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/4,3/4,1/4,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/2,0,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,0,0,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,3/4,1/4,1/4,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/2,0,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,0,1/2,0,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,1,0,0,3/4,1/4,1/4,1,-1,0,0,0,0,0,-1,0,0,1,0,0,0,1/2,0,1,1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/4,1/4,3/4,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,1/2,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,3/4,3/4,3/4,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/4,3/4,1/4,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,3/4,3/4,3/4,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,1,0,0,3/4,3/4,3/4,1,0,0,1,0,0,1,0,0,1,0,0,0,0,1/2,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,1/2,1];
                case 229; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,1,0,0,0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,-1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,1,0,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,-1,0,0,0,0,0,0,1,1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/2,1/2,1/2,1];
                case 230; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,1/2,1,-1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,1/2,1,-1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,1/2,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,3/4,1/4,1/4,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,1/2,1,1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,1/2,1,1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,1/2,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,1,0,1/4,3/4,3/4,1,-1,0,0,0,0,-1,0,0,0,0,1,0,0,1/2,0,1,-1,0,0,0,0,1,0,0,0,0,-1,0,1/2,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1/2,1/2,1/2,1,0,1,0,0,1,0,0,0,0,0,-1,0,1/4,3/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,1/2,1/2,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1/2,0,1/2,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,3/4,3/4,3/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,0,1/2,1/2,1,0,1,0,0,-1,0,0,0,0,0,1,0,1/4,3/4,1/4,1,0,-1,0,0,0,0,-1,0,1,0,0,0,1/2,1/2,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,0,1/4,3/4,1/4,1,0,-1,0,0,-1,0,0,0,0,0,1,0,3/4,1/4,1/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,1/4,1/4,3/4,1,1,0,0,0,0,0,-1,0,0,1,0,0,3/4,1/4,1/4,1,1,0,0,0,0,1,0,0,0,0,-1,0,0,1/2,0,1,1,0,0,0,0,-1,0,0,0,0,1,0,1/2,0,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,1,0,0,0,0,1,0,1/2,1/2,0,1,0,1,0,0,0,0,-1,0,1,0,0,0,1/2,0,1/2,1,0,1,0,0,1,0,0,0,0,0,1,0,1/4,1/4,1/4,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,1/2,1/2,1,0,-1,0,0,1,0,0,0,0,0,-1,0,3/4,1/4,3/4,1,0,1,0,0,0,0,1,0,-1,0,0,0,1/2,1/2,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,1,0,0,0,3/4,1/4,3/4,1,0,1,0,0,-1,0,0,0,0,0,-1,0,3/4,3/4,1/4,1,-1,0,0,0,0,0,1,0,0,-1,0,0,1/4,3/4,3/4,1,1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,1/2,1,0,-1,0,0,0,0,1,0,-1,0,0,0,0,1/2,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/4,1/4,1/4,1,0,1,0,0,0,0,-1,0,-1,0,0,0,1/2,0,0,1,0,1,0,0,-1,0,0,0,0,0,1,0,3/4,1/4,3/4,1,0,-1,0,0,0,0,-1,0,1,0,0,0,0,0,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,1/2,1/2,1/2,1,0,0,1,0,0,1,0,0,-1,0,0,0,3/4,1/4,3/4,1,0,-1,0,0,1,0,0,0,0,0,1,0,3/4,3/4,1/4,1,1,0,0,0,0,0,-1,0,0,1,0,0,1/4,3/4,3/4,1,0,0,1,0,-1,0,0,0,0,-1,0,0,1/2,0,1/2,1,0,0,1,0,0,-1,0,0,1,0,0,0,1/4,1/4,3/4,1,0,1,0,0,1,0,0,0,0,0,1,0,3/4,3/4,3/4,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,3/4,3/4,3/4,1,0,0,-1,0,-1,0,0,0,0,1,0,0,0,1/2,1/2,1,0,0,-1,0,0,1,0,0,1,0,0,0,3/4,1/4,1/4,1,0,-1,0,0,1,0,0,0,0,0,-1,0,1/4,3/4,1/4,1,-1,0,0,0,0,0,1,0,0,1,0,0,1/4,3/4,1/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,1/2,1/2,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,3/4,3/4,3/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,1/4,1/4,3/4,1,0,0,-1,0,0,-1,0,0,1,0,0,0,1/4,3/4,1/4,1,0,1,0,0,-1,0,0,0,0,0,-1,0,1/4,1/4,3/4,1,-1,0,0,0,0,0,1,0,0,-1,0,0,3/4,1/4,1/4,1,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,1/2,1,0,1,0,0,0,0,-1,0,1,0,0,0,0,1/2,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,1/2,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,0,0,1/2,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,1/2,0,1/2,1,0,0,-1,0,0,1,0,0,-1,0,0,0,3/4,3/4,1/4,1,1,0,0,0,0,0,1,0,0,1,0,0,1/4,1/4,1/4,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,1/2,1/2,1,0,0,1,0,0,-1,0,0,-1,0,0,0,1/4,3/4,3/4,1,1,0,0,0,0,0,-1,0,0,-1,0,0,3/4,1/4,3/4,1,0,0,1,0,-1,0,0,0,0,1,0,0,1/2,1/2,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1/4,1/4,1/4,1,-1,0,0,0,0,0,-1,0,0,1,0,0,3/4,3/4,1/4,1,0,0,1,0,-1,0,0,0,0,-1,0,0,0,1/2,0,1,0,0,1,0,0,-1,0,0,1,0,0,0,3/4,3/4,1/4,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/4,1/4,1/4,1,0,0,-1,0,-1,0,0,0,0,1,0,0,1/2,0,0,1,0,0,-1,0,0,1,0,0,1,0,0,0,1/4,3/4,3/4,1,-1,0,0,0,0,0,1,0,0,1,0,0,3/4,1/4,3/4,1,0,0,-1,0,1,0,0,0,0,-1,0,0,0,0,1/2,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/4,1/4,1/4,1,1,0,0,0,0,0,1,0,0,-1,0,0,3/4,3/4,1/4,1,0,0,-1,0,0,1,0,0,-1,0,0,0,1/4,1/4,3/4,1,1,0,0,0,0,0,1,0,0,1,0,0,3/4,3/4,3/4,1,0,0,1,0,0,-1,0,0,-1,0,0,0,3/4,1/4,1/4,1,1,0,0,0,0,0,-1,0,0,-1,0,0,1/4,3/4,1/4,1,0,0,1,0,0,1,0,0,1,0,0,0,3/4,3/4,3/4,1,-1,0,0,0,0,0,-1,0,0,1,0,0,1/4,1/4,3/4,1,0,0,1,0,-1,0,0,0,0,1,0,0,0,0,1/2,1,0,0,-1,0,1,0,0,0,0,1,0,0,0,1/2,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,1/2,0,0,1];
                case 231; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1];
                case 232; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,0,0,1];
                case 233; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,1];
                case 234; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1];
                case 235; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,1/2,1];
                case 236; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1];
                case 237; S=[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1/2,1/2,1/2,1,0,-1,0,0,0,0,-1,0,-1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1/2,1/2,1/2,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,0,1/2,1/2,1/2,1,-1,0,0,0,0,0,-1,0,0,-1,0,0,1/2,1/2,1/2,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1/2,1/2,1/2,1,1,0,0,0,0,0,1,0,0,1,0,0,1/2,1/2,1/2,1];
                end
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

        function R            = generate_pg(pg_code,from_memory)
            
            import am_lib.* am_dft.*
            
            if nargin<2; from_memory=true; end

            if from_memory
                % generated with this code:
                % import am_dft.*
                % for i = 1:32
                %     R{i} = round(generate_pg(i,false));
                %     x(i) = identify_pointgroup(R{i});
                %     d{i} = decode_pg(x(i));
                %     fprintf('case %i; R = [',i); fprintf('%g',R{i}(1)); fprintf(',%g',R{i}(2:end)); fprintf(']; %% %s\n',d{i});
                % end
                
                switch pg_code
                case 1;  R = [1,0,0,0,1,0,0,0,1]; % c_1
                case 2;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1]; % s_2
                case 3;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1]; % c_2
                case 4;  R = [1,0,0,0,1,0,0,0,1,1,0,0,0,-1,0,0,0,1]; % c_1h
                case 5;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,-1,0,0,0,-1,0,0,0,-1,1,0,0,0,-1,0,0,0,1]; % c_2h
                case 6;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,1,0,0,0,-1,0,0,0,-1]; % d_2
                case 7;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,1,0,0,0,-1,0,0,0,1,-1,0,0,0,1,0,0,0,1]; % c_2v
                case 8;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,-1,0,0,0,-1,0,0,0,-1,1,0,0,0,-1,0,0,0,-1,1,0,0,0,1,0,0,0,-1,1,0,0,0,-1,0,0,0,1,-1,0,0,0,1,0,0,0,1]; % d_2h
                case 9;  R = [1,0,0,0,1,0,0,0,1,0,1,0,-1,-1,0,0,0,1,-1,-1,0,1,0,0,0,0,1]; % c_3
                case 10; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,0,1,0,-1,-1,0,0,0,1,0,-1,0,1,1,0,0,0,-1,-1,-1,0,1,0,0,0,0,1,1,1,0,-1,0,0,0,0,-1]; % s_6
                case 11; R = [1,0,0,0,1,0,0,0,1,0,1,0,1,0,0,0,0,-1,0,1,0,-1,-1,0,0,0,1,1,0,0,-1,-1,0,0,0,-1,-1,-1,0,0,1,0,0,0,-1,-1,-1,0,1,0,0,0,0,1]; % d_3
                case 12; R = [1,0,0,0,1,0,0,0,1,0,-1,0,-1,0,0,0,0,1,0,1,0,-1,-1,0,0,0,1,-1,0,0,1,1,0,0,0,1,1,1,0,0,-1,0,0,0,1,-1,-1,0,1,0,0,0,0,1]; % c_3v
                case 13; R = [1,0,0,0,1,0,0,0,1,0,-1,0,-1,0,0,0,0,-1,-1,0,0,0,-1,0,0,0,-1,0,1,0,-1,-1,0,0,0,1,0,1,0,1,0,0,0,0,1,-1,0,0,1,1,0,0,0,-1,0,-1,0,1,1,0,0,0,-1,1,1,0,0,-1,0,0,0,-1,-1,-1,0,1,0,0,0,0,1,1,0,0,-1,-1,0,0,0,1,-1,-1,0,0,1,0,0,0,1,1,1,0,-1,0,0,0,0,-1]; % d_3d
                case 14; R = [1,0,0,0,1,0,0,0,1,0,1,0,-1,0,0,0,0,1,-1,0,0,0,-1,0,0,0,1,0,-1,0,1,0,0,0,0,1]; % c_4
                case 15; R = [1,0,0,0,1,0,0,0,1,0,-1,0,1,0,0,0,0,-1,-1,0,0,0,-1,0,0,0,1,0,1,0,-1,0,0,0,0,-1]; % s_4
                case 16; R = [1,0,0,0,1,0,0,0,1,0,1,0,-1,0,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,-1,0,0,0,-1,0,0,0,1,0,-1,0,1,0,0,0,0,-1,0,-1,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,-1,0,1,0,-1,0,0,0,0,-1]; % c_4h
                case 17; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,1,0,-1,0,0,0,0,1,0,1,0,1,0,0,0,0,-1,0,-1,0,-1,0,0,0,0,-1,-1,0,0,0,-1,0,0,0,1,0,-1,0,1,0,0,0,0,1,1,0,0,0,-1,0,0,0,-1]; % d_4
                case 18; R = [1,0,0,0,1,0,0,0,1,0,1,0,-1,0,0,0,0,1,1,0,0,0,-1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,0,1,0,1,0,0,0,0,1,0,-1,0,-1,0,0,0,0,1,0,-1,0,1,0,0,0,0,1,-1,0,0,0,1,0,0,0,1]; % c_4v
                case 19; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,-1,0,1,0,0,0,0,-1,0,-1,0,-1,0,0,0,0,1,0,1,0,1,0,0,0,0,1,-1,0,0,0,-1,0,0,0,1,0,1,0,-1,0,0,0,0,-1,1,0,0,0,-1,0,0,0,-1]; % d_2d
                case 20; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,1,0,-1,0,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,0,1,0,1,0,0,0,0,-1,1,0,0,0,-1,0,0,0,1,0,-1,0,-1,0,0,0,0,-1,-1,0,0,0,-1,0,0,0,1,0,-1,0,1,0,0,0,0,-1,0,-1,0,1,0,0,0,0,1,1,0,0,0,-1,0,0,0,-1,0,-1,0,-1,0,0,0,0,1,0,1,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,-1,0,1,0,-1,0,0,0,0,-1,-1,0,0,0,1,0,0,0,1]; % d_4h
                case 21; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,0,1,0,-1,-1,0,0,0,1,0,-1,0,1,1,0,0,0,1,-1,-1,0,1,0,0,0,0,1,1,1,0,-1,0,0,0,0,1]; % c_6
                case 22; R = [1,0,0,0,1,0,0,0,1,1,0,0,0,1,0,0,0,-1,0,1,0,-1,-1,0,0,0,1,0,1,0,-1,-1,0,0,0,-1,-1,-1,0,1,0,0,0,0,1,-1,-1,0,1,0,0,0,0,-1]; % c_3h
                case 23; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,0,1,0,-1,-1,0,0,0,1,1,0,0,0,1,0,0,0,-1,0,-1,0,1,1,0,0,0,1,0,-1,0,1,1,0,0,0,-1,-1,-1,0,1,0,0,0,0,1,0,1,0,-1,-1,0,0,0,-1,1,1,0,-1,0,0,0,0,1,1,1,0,-1,0,0,0,0,-1,-1,-1,0,1,0,0,0,0,-1]; % c_6h
                case 24; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,0,1,0,1,0,0,0,0,-1,0,1,0,-1,-1,0,0,0,1,0,-1,0,-1,0,0,0,0,-1,0,-1,0,1,1,0,0,0,1,1,0,0,-1,-1,0,0,0,-1,-1,-1,0,0,1,0,0,0,-1,-1,-1,0,1,0,0,0,0,1,-1,0,0,1,1,0,0,0,-1,1,1,0,0,-1,0,0,0,-1,1,1,0,-1,0,0,0,0,1]; % d_6
                case 25; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,0,-1,0,-1,0,0,0,0,1,0,1,0,-1,-1,0,0,0,1,0,1,0,1,0,0,0,0,1,0,-1,0,1,1,0,0,0,1,-1,0,0,1,1,0,0,0,1,1,1,0,0,-1,0,0,0,1,-1,-1,0,1,0,0,0,0,1,1,0,0,-1,-1,0,0,0,1,-1,-1,0,0,1,0,0,0,1,1,1,0,-1,0,0,0,0,1]; % c_6v
                case 26; R = [1,0,0,0,1,0,0,0,1,1,0,0,0,1,0,0,0,-1,0,-1,0,-1,0,0,0,0,1,0,1,0,-1,-1,0,0,0,1,0,-1,0,-1,0,0,0,0,-1,0,1,0,-1,-1,0,0,0,-1,-1,0,0,1,1,0,0,0,1,1,1,0,0,-1,0,0,0,1,-1,-1,0,1,0,0,0,0,1,-1,0,0,1,1,0,0,0,-1,1,1,0,0,-1,0,0,0,-1,-1,-1,0,1,0,0,0,0,-1]; % d_3h
                case 27; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,0,1,0,1,0,0,0,0,-1,0,1,0,-1,-1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,0,-1,0,-1,0,0,0,0,-1,0,-1,0,1,1,0,0,0,1,1,0,0,0,1,0,0,0,-1,1,0,0,-1,-1,0,0,0,-1,0,-1,0,-1,0,0,0,0,1,-1,-1,0,0,1,0,0,0,-1,-1,-1,0,1,0,0,0,0,1,0,-1,0,1,1,0,0,0,-1,-1,0,0,1,1,0,0,0,-1,0,1,0,1,0,0,0,0,1,1,1,0,0,-1,0,0,0,-1,1,1,0,-1,0,0,0,0,1,0,1,0,-1,-1,0,0,0,-1,-1,0,0,1,1,0,0,0,1,1,1,0,0,-1,0,0,0,1,1,1,0,-1,0,0,0,0,-1,1,0,0,-1,-1,0,0,0,1,-1,-1,0,0,1,0,0,0,1,-1,-1,0,1,0,0,0,0,-1]; % d_6h
                case 28; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,1,0,0,0,1,1,0,0,0,1,0,0,0,-1,-1,0,0,0,-1,0,0,0,1,-1,0,0,0,0,1,1,0,0,0,1,0,0,-1,0,0,0,-1,1,0,0,0,0,-1,-1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,1,0,0,0,-1,0,-1,0,0,0,-1,0,0,0,1,1,0,0,0,-1,0,0,0,-1]; % t
                case 29; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,1,0,0,0,1,1,0,0,-1,0,0,0,-1,0,0,0,-1,0,1,0,0,0,-1,-1,0,0,1,0,0,0,-1,0,0,0,1,0,-1,0,0,0,1,-1,0,0,0,0,1,1,0,0,0,1,0,0,-1,0,0,0,-1,-1,0,0,0,-1,0,0,0,-1,1,0,0,0,0,-1,-1,0,0,0,1,0,0,-1,0,0,0,1,1,0,0,0,0,1,-1,0,0,0,-1,0,0,1,0,0,0,-1,1,0,0,0,0,-1,1,0,0,0,-1,0,0,0,-1,-1,0,0,0,-1,0,0,1,0,0,0,1,-1,0,0,0,0,1,1,0,0,0,-1,0,-1,0,0,0,-1,0,0,0,1,0,0,-1,1,0,0,0,1,0,1,0,0,0,-1,0,0,0,-1,0,0,1,-1,0,0,0,1,0,1,0,0,0,1,0,0,0,-1,-1,0,0,0,1,0,0,0,1]; % t_h
                case 30; R = [1,0,0,0,1,0,0,0,1,0,1,0,0,0,1,1,0,0,0,1,0,-1,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,0,1,0,-1,0,1,0,0,-1,0,0,0,0,1,0,1,0,-1,0,0,0,-1,0,0,0,1,1,0,0,0,0,-1,0,1,0,0,-1,0,1,0,0,0,0,1,0,-1,0,0,0,-1,1,0,0,0,0,1,0,1,0,-1,0,0,0,-1,0,0,0,1,-1,0,0,0,0,-1,0,1,0,1,0,0,0,0,-1,-1,0,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,-1,1,0,0,0,-1,0,0,1,0,0,0,-1,-1,0,0,0,0,1,-1,0,0,0,-1,0,0,1,0,1,0,0,0,0,-1,-1,0,0,0,1,0,0,0,-1,1,0,0,0,-1,0,0,0,-1,0,0,-1,0,-1,0,-1,0,0,-1,0,0,0,0,-1,0,-1,0,0,-1,0,-1,0,0,0,0,-1]; % o
                case 31; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1,1,0,0,1,0,0,0,-1,0,0,0,-1,0,-1,0,-1,0,0,0,0,1,0,-1,0,0,0,1,-1,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,-1,-1,0,0,0,-1,0,1,0,0,0,0,-1,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,-1,1,0,0,0,0,1,0,1,0,1,0,0,0,0,1,1,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1,0,-1,0,-1,0,0,0,0,1,-1,0,0,0,-1,0,-1,0,0,0,0,-1,0,1,0,0,0,-1,0,1,0,-1,0,0,0,0,-1,-1,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,-1,0,-1,0,1,0,0,0,0,-1,1,0,0,0,-1,0]; % t_d
                case 32; R = [1,0,0,0,1,0,0,0,1,0,1,0,0,0,1,1,0,0,0,1,0,-1,0,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,0,0,1,1,0,0,0,1,0,0,0,1,0,-1,0,1,0,0,0,-1,0,0,0,-1,-1,0,0,-1,0,0,0,0,1,0,1,0,-1,0,0,0,-1,0,0,0,1,0,-1,0,1,0,0,0,0,-1,1,0,0,0,0,-1,0,1,0,0,0,-1,-1,0,0,0,-1,0,0,-1,0,1,0,0,0,0,1,0,-1,0,0,0,-1,1,0,0,0,0,-1,0,1,0,-1,0,0,0,0,1,0,1,0,-1,0,0,1,0,0,0,0,-1,0,-1,0,0,-1,0,0,0,1,-1,0,0,1,0,0,0,1,0,0,0,-1,0,0,-1,0,1,0,1,0,0,0,0,-1,-1,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,1,0,0,0,0,1,0,-1,0,0,1,0,-1,0,0,0,0,-1,0,0,-1,1,0,0,0,-1,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,-1,0,0,0,0,-1,0,-1,0,1,0,0,0,0,1,-1,0,0,0,-1,0,0,1,0,0,0,-1,1,0,0,0,1,0,1,0,0,0,0,-1,0,0,1,0,-1,0,-1,0,0,-1,0,0,0,1,0,0,0,-1,0,0,1,1,0,0,0,-1,0,-1,0,0,0,0,-1,0,1,0,1,0,0,0,-1,0,0,0,-1,0,0,1,-1,0,0,0,1,0,0,0,-1,0,-1,0,-1,0,0,0,-1,0,0,0,1,1,0,0,-1,0,0,0,0,-1,0,-1,0,0,0,-1,1,0,0,0,1,0,0,-1,0,-1,0,0,0,0,1,1,0,0,0,-1,0,0,0,1,0,-1,0,-1,0,0,0,0,-1,-1,0,0,0,1,0,0,0,1,0,0,1,0,1,0,1,0,0,1,0,0,0,0,1,0,1,0,0,1,0,1,0,0,0,0,1]; % o_h
                end
                R = reshape(R,3,3,[]);
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
            if s(1)==4 && s(2)==4
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

        function abc_angles   = bas2abc(bas)
            % [a,b,c,alpha,beta,gamma] = bas2abc(bas)
            M = bas.' * bas;
            % de Graef p 86
            a = sqrt(M(1,1));
            b = sqrt(M(2,2));
            c = sqrt(M(3,3));
            alpha= acosd(M(2,3)/(b*c));
            beta = acosd(M(1,3)/(a*c));
            gamma= acosd(M(1,2)/(a*b));
            abc_angles = [a,b,c,alpha,beta,gamma];
        end

        function bas          = abc2bas(abc,brav)
            
            import am_dft.abc2bas
            
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

        function [bjc,gen]    = find_isomorphic_bijection(G,H)
            % bjc = find_isomorphic_bijection(g,h), g and h are multiplication tables 
            % isomorphism: g <==> h(:,:,inds) , returns inds = [] if groups are not isomorphic
            % % check with:
            % G = get_multiplication_table(g);
            % H = get_multiplication_table(h);
            % i=1; G - relabel_multiplication_table(H,bjc(:,i))

            import am_lib.* am_dft.*

            % number of symmetries
            nsyms = size(G,1);
            
            % sorting if desired (changes order of bjc but does not affect final result)
            Gfwd = [1:nsyms]; Hfwd = [1:nsyms];

            % get generators
            [gg,Gg] = get_generators(G); ngenerators = numel(gg);
            [hg,Hg] = get_generators_list(H,ngenerators);
            
            % compare multiplication tables, looking for bijections between symmetries g and h, i.e. isomorphism between G and H
            j=0; bjc=[]; gen=[];
            for i = 1:size(Hg,2)
                % fwd = relabeling of H to match G
                fwd(Gg) = Hg(:,i);
                if all(all( G == relabel_multiplication_table(H,fwd) ))
                    % isomorphism found,  save bijection which takes symmetries in H to G in original order
                    j=j+1; bjc(Gfwd,j) = Hfwd(fwd); gen(:,j) = hg(:,i);
                    % save
                end
            end
            nbjcs = j;

            % relabel generators for h in terms of g
            for j = 1:nbjcs
                bak(bjc(:,j)) = [1:nsyms]; gen(:,j) = bak(gen(:,j));
            end

        end

        function [fwd,M]      = find_pointgroup_transformation(g,h,algo)
            % finds the permutation fwd and rotation M which transforms point group g into h:
            % matmul_(matmul_(M(:,:,k),g(:,:,:)),inv(M(:,:,k))) == h(:,:,bjc(:,k)) for any k

            import am_dft.* am_lib.*
            
            % algo = 1, returns only 1 bijection of many possibilities
            if nargin<3; algo=0; end

            % find bijection bjc: h -> g and generators for h
            G=get_multiplication_table(g);
            H=get_multiplication_table(h);
            [bjc,gen] = find_isomorphic_bijection(G,H);
            srt_ = rankc_(bjc); bjc = bjc(:,srt_); gen = gen(:,srt_);

            % check diagonalize symmetry values to make sure the symmetries are of the correct type
            % this removes the possibility that multiplication table may be obeyed with symmetries of differnet types
            nbjcs=size(bjc,2); ex_ = true(1,nbjcs); 
            g_ps_id = identify_point_symmetries(g); 
            h_ps_id = identify_point_symmetries(h);
            for j = 1:nbjcs;   ex_ = all(h_ps_id(bjc(:,j))==g_ps_id,1); end
            bjc = bjc(:,ex_); gen = gen(:,ex_); nbjcs=size(bjc,2);

            % find integer transformation matrix (try all possible bijections and
            % integer transformations with nonzero-determinants); in reality, not all
            % matrices are tried, instead: 
            %       first try the identity
            %       then try all matrices with elements in [-1,0,1] with determinant 1 
            %       then try all matrices with elements in [-2,-1,0,1,2] with determinant 1 
            q = 0; ngens = size(gen,1); 
            for m = 0:1
                % generate the possible transformation matrices
                switch m
                    case 0
                        % first try the identity
                        X = eye(3); nXs = 1;
                    otherwise
                        % then try other integer matrices with nonzero determinants
                        N=9; Q=[-m:m]; nQs=numel(Q);[Y{N:-1:1}]=ndgrid(1:nQs); X=reshape(Q(reshape(cat(N+1,Y{:}),[],N)).',3,3,[]);
                        ex_ = false(1,size(X,3)); for i = 1:size(X,3); ex_(i) = abs(det(X(:,:,i)))<am_lib.eps; end
                        X = X(:,:,~ex_); nXs = size(X,3);
                        % sort it to make it nice (put positive values first)
                        X=X(:,:,rankc_( -[max(reshape(X,9,[]));min(reshape(X,9,[]));sum(reshape(X,9,[]),1)] ));
                end
                % get inverse elements
                invX = zeros(3,3,nXs); for i = 1:nXs; invX(:,:,i) = inv(X(:,:,i)); end

                % find similarity transform which converts all symmetries from g to h by
                % checking each matrix generated above
                for k = 1:nXs; for j = 1:nbjcs
                    gencheck_ = true;
                    % check generators
                    for i = 1:ngens
                        if any(any(abs( X(:,:,k) * g(:,:,gen(i,j)) * inv(X(:,:,k)) - h(:,:,bjc(gen(i,j),j)) )>am_lib.eps ))
                            gencheck_ = false; break;
                        end
                    end
                    % check whether a match has been found
                    if gencheck_
                        % at this point a bijection permutation and transformation matrix have been
                        % found such that: matmul_(matmul_(X(:,:,k),g(:,:,:)),invX(:,:,k)) == h(:,:,bjc(:,j))
                        % save it 
                        q = q+1;
                        M(:,:,q) = X(:,:,k); 
                        fwd(:,q) = bjc(:,j);
                        % only one is needed
                        if algo==1
                            return; 
                        end
                    end
                end;end
            end
        end

        
        % DEVELOPMENT
        
            function MT = get_mt_(r_mt,t_mt,rt)
                
                import am_lib.*
                 
                % define multiplication
                mult_ = @(a,b) [r_mt(a(1),b(1));t_mt(a(2),b(2))];
                %
                nsyms=size(rt,2);
                % construct multiplication table
                for j1 = 1:nsyms
                for j2 = 1:nsyms
                    MT(j1,j2) = member_(mult_(rt(:,j1),rt(:,j2)),rt);
                end
                end
            end

        % DEVELOPMENT
        

        % unit cells
        
        function [uc,pc,ic,cc]= get_cell(flag, arg, tol)
            % [uc,pc,ic,cc]   = get_cell(fposcar)
            % fposcar can be a poscar or cif file

            import am_lib.* am_dft.*
            
            % time
            fprintf(' ... getting cell'); tic
            
            % validate input
            validatestring(flag,{'poscar','cif','material'});
            switch flag
                case {'poscar','cif'};  if exist(arg,'file')~=2; fprintf('\n'); error('File does not exist: %s',arg); end
            end
            switch flag
                case 'poscar';   uc = load_poscar(arg);
                case 'cif';      uc = load_cif(arg);
                case 'material'; uc = load_material(arg);
                case 'create';   uc = am_dft.create_cell(arg{:});
            end
            
            % set default numerical tolerance
            if nargin < 3; tol = am_dft.tiny; end

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

            % save bas2pc and tau2pc to convert [uc-frac] to [pc-frac]
            uc.bas2pc = pc.bas/uc.bas; uc.tau2pc = pc.bas\uc.bas;
            
            % save bas2pc and tau2pc to convert [uc-frac] to [vc-frac]
            cc.bas2pc = pc.bas/cc.bas; cc.tau2pc = pc.bas\cc.bas;

            % print basic symmetry info
            [~,H,~,R] = get_symmetries(pc, tol);
            bv_code = identify_bravais_lattice(pc.bas, tol);
            
            % holohodry should give same info as bravais lattice
            hg_code = identify_pointgroup(H); 
            pg_code = identify_pointgroup(R); 
            sg_code = identify_spacegroup(pg_code); % BETA
            
            % atomic density [g/cm3]
            mass_density = sum(uc.mass(uc.species)) / 6.02214E23 / det(uc.bas*0.1 * 1E-7);
            number_density = uc.natoms / det(uc.bas*0.1);
            formula_density = number_density/(uc.natoms/gcd_(uc.nspecies));

            % print relevant information
            verbose = true;
            if verbose
                fprintf(' (%.3f s) \n',toc);
                fprintf('     %-15s = %s\n','formula',get_formula(uc));
                fprintf('     %-15s = %s\n','primitive',decode_bravais(bv_code));
                fprintf('     %-15s = %s\n','holohodry',decode_holohodry(hg_code));
                fprintf('     %-15s = %s\n','point group',decode_pg(pg_code));
                fprintf('     %-15s = %s\n','space group',strrep(cell2mat(join(decode_sg(sg_code),',')),' ',''));
                fprintf('     %-15s = %-8.3f [g/cm3] \n','mass density',mass_density);
                fprintf('     %-15s = %-8.3f [atoms/nm3]\n','number density',number_density);
                fprintf('     %-15s = %-8.3f [f.u./nm3]\n','formula density',formula_density);
            end
            
            % sub functions
            
            function [uc]    = load_cif(fcif)

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
                            while and(~any(strmatchi_(str{j+i},{'_','#','loop_'})),lt(i+j,nlines))
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
                [~,~,species]=unique(species_symb,'stable'); species=species(:).'; 
                % make sure there is only one atom at each site
                [~,i]=am_lib.uniquec_(tau); tau=tau(:,i); species=species(i); symb=species_symb(i); Z=am_dft.get_atomic_number(symb);

                % get space group symmetries
                position_aliases = {...
                    '_symmetry_Int_Tables_number', ... % icsd
                    '_symmetry_equiv_pos_as_xyz', ... % crystal maker
                    };
                for alias = position_aliases; ex_ = contains(str,alias{:}); if any(ex_)
                    switch alias{:}
                        case '_symmetry_equiv_pos_as_xyz'
                            j=find(ex_); i=1; syms x y z;
                            while ~isempty(strtrim(str{j+i})) && ~any(strmatchi_(str{j+i},{'_','#','loop_'})) && lt(i+j,nlines)
                                buffer = strsplit(str{j+i},'''');
                                buffer = strsplit(buffer{2},','); 
                                tmp=coeffs(eval(buffer{1}),'all'); T(1) = double(tmp(end));
                                tmp=coeffs(eval(buffer{2}),'all'); T(2) = double(tmp(end));
                                tmp=coeffs(eval(buffer{3}),'all'); T(3) = double(tmp(end));
                                S(1:3,1:3,i) = double(equationsToMatrix(...
                                    [eval(buffer{1}),eval(buffer{2}),eval(buffer{3})] ));
                                S(1:3,4:4,i) = mod_(T(:));
                                S(4:4,1:3,i) = 0;
                                S(4:4,4:4,i) = 1;
                                i=i+1;
                            end
                            break;
                        case '_symmetry_Int_Tables_number'
                            % true == generate from memory
                            sg_id = am_lib.extract_token_(str,'_symmetry_Int_Tables_number',true); S = am_dft.generate_sg(sg_id,true);
                            break;
                    end
                end; end
                nSs = size(S,3);

                % generate all positions based on symmetry operations
                seitz_apply_ = @(S,tau) am_lib.mod_(reshape(am_lib.matmul_(S(1:3,1:3,:),tau),3,[],size(S,3)) + S(1:3,4,:));
                tau = reshape(am_lib.mod_(seitz_apply_(S,tau)),3,[]); species=repmat(species,1,nSs);
                [tau, j] = am_lib.uniquec_(tau); species=species(j);
                [species,i]=sort(species); tau=tau(:,i);

                % define primitive cell creation function and make structure
                uc_ = @(bas,symb,Z,species) struct('units','frac','bas',bas, ...
                    'symb',{symb},'mass',am_dft.get_atomic_mass(Z),'nspecies',sum(unique(species).'==species,2).', ...
                    'natoms',numel(species),'tau',tau,'species',species);
                uc = uc_(bas,symb,Z,species);
            end
            
            function [uc]    = load_poscar(fposcar)
                fid=fopen(fposcar,'r');                % open file
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
                    for i=1:length(uc.nspecies); uc.mass(i)=am_dft.get_atomic_mass(am_dft.get_atomic_number(uc.symb{i})); end
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
            end

            function [uc]    = load_material(material)
                % conventional cell
                switch material
                    case 'VN'; uc = am_dft.create_cell(am_dft.abc2bas(4.1340,'cubic'),[[0;0;0], [1;1;1]/2], {'V','N'}, 225);
                    case 'Si'; uc = am_dft.create_cell(am_dft.abc2bas(5.4209,'cubic'),[0;0;0], {'Si'}, 227);
                    otherwise; error('load_material: unknown material');
                end
            end
            
        end

        function [dc,idc]     = get_displaced_cell(pc,bvk,n,kpt,amp,mode,nsteps)
            % Note: can set mode=[] for interactive selection
            % n=[4;4;4]; kpt=[0;0;1/4]; amp=10; mode=6; nsteps=51;
            % [~,md] = get_displaced_cell(pc,bvk,n,kpt,amp,mode,nsteps);
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

            import am_lib.* am_dft.*

            % initialize figure
            set(gcf,'color','w'); hold on;

            % plot atoms
            clist=color_(max(pc.species));
            for i = 1:max(pc.species)
                ex_ = pc.species==i; radius = ones(1,sum(ex_)) * get_crystal_radius(get_atomic_number(pc.symb{i})) * 5000;
                h(i) = scatter3_(pc.bas*pc.tau(:,ex_), radius,'MarkerEdgeColor','k','MarkerFaceColor',clist(i,:));
            end
            
            % plot pc boundaries
            plothull_(pc.bas*[0,1,0,1,0,1,0,1;0,0,1,1,0,0,1,1;0,0,0,0,1,1,1,1]);
            
            hold off; daspect([1 1 1]); box on;

            % legend
            legend(h,pc.symb{:},'boxoff'); axis off;

        end

        function [F]          = plot_md_cell(md, varargin)
            % n=[4;4;4]; kpt=[0;0;1/4]; amp=10; mode=6; nsteps=51;
            % [~,md] = get_displaced_cell(pc,bvk,n,kpt,amp,mode,nsteps);
            % clf; [F]=plot_md_cell(md,'view',[0;1;0]); movie(F,3); % write_poscar(md,'POSCAR_test')

            import am_lib.* am_dft.*

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

        function [T]          = coincidence_site_lattice(B1, B2, multiplicity, tol)
            % T is the matrix which transforms lattice B1 into B2 through: B2 == B1*T
            
            import am_lib.rankc_ am_lib.matmul_ am_dft.get_niggli_basis am_dft.bas2abc
            
            fprintf(' ... getting coincidence site lattice'); tic
            
            % determinant of upper triangular matrix
            det_up_tri_ = @(A) A(1,1)*A(2,2)*A(3,3);

            % generate all possible transformation matrices Q
            % requirepment 0 <= Qij < Qjj for i<j
            %  i1 i6 i5      11 12 13
            %  0  i2 i4      21 22 23
            %  0  0  i3      31 32 33
            m = multiplicity; k=0;
            for i1 = 1:m; for i2 = 1:m; for i3 = 1:m
            for i4 = 1:i3;for i5 = 1:i3;for i6 = 1:i2
                k=k+1; 
                Q(:,:,k) = [i1 i6 i5; 0 i2 i4; 0 0 i3];
                d(k) = det_up_tri_(Q(:,:,k));
            end; end; end
            end; end; end
            ex_  = (d<=multiplicity); Q=Q(:,:,ex_); d=d(ex_);
            inds = rankc_(d); Q=Q(:,:,inds); d=d(inds);

            X1=matmul_(B1,Q); X2=matmul_(B2,Q); k=0;
            n1s=size(X1,3);n2s=size(X2,3); v=zeros(3,n1s*n2s);
            % find metric tensors which match
            for i1 = 1:n1s
                M1 = bas2abc(get_niggli_basis(X1(:,:,i1)));
            for i2 = 1:n2s
                M2 = bas2abc(get_niggli_basis(X2(:,:,i2)));
                % save info
                k=k+1; v(:,k) = [sum(abs(M1-M2)),i1,i2];
            end
            end
            % sort coincidence site lattice based on mismatch
            inds = rankc_( v(1,:) ); v = v(:,inds);

            % go back to the inds and recalculate the transformation to get from B1 to B2
            ninds=numel(inds); T=zeros(3,3,ninds);
            for i = 1:5%ninds
                i1 = v(2,i); i2 = v(3,i);
                % get Niggli basis
                [N1,T1]=get_niggli_basis(X1(:,:,i1));
                [N2,T2]=get_niggli_basis(X2(:,:,i2));
                % find transformation T which takes N1 to N2 (i.e., N1 = N2*T), using:
                % N1=B1*Q(:,:,i1)*T1, N2=B2*Q(:,:,i2)*T2, that is: B2*Q(:,:,i2)*T2 = B1*Q(:,:,i1)*T1 * T
                % B2 == B1*T
                T(:,:,i) = Q(:,:,i1)*T1*(N2\N1)/T2/Q(:,:,i2);
            end

            % print results
            fprintf(' (%.3f s) \n',toc);
            fprintf('     %-10s %-7s %-7s %5s %5s %5s %5s %5s %5s %5s %5s %5s  \n','mismatch','mult1','mult2','T11','T12','T13','T21','T22','T23','T31','T32','T33');
            for i = 1:5%ninds
            fprintf('     %-10.3f %-7i %-7i %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f \n',v(1,i),d(v(2,i)),d(v(3,i)),T(:,:,i));
            end
        end
        
        function [bas,T]      = get_niggli_basis(bas, tol)
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
            
            if nargin<2; tol=am_dft.tiny; end
            
            T = eye(3); [x,y,z,a,b,c,l,m,n,bas] = update_(bas, eye(3), tol);
           
            % begin reduction procedure
            for counter = 1:100
                if (a > b + tol || (~ abs(a - b) > tol && abs(x) >  abs(y) + tol))
                    % Procedure A1
                    A = [0, -1, 0; -1, 0, 0; 0, 0, -1];
                    [x,y,z,a,b,c,l,m,n,bas] = update_(bas, A, tol); T = T * A;
                end
                if (b > c + tol || (~ abs(b - c) > tol && abs(y) >  abs(z) + tol))
                    % Procedure A2
                    A = [-1, 0, 0; 0, 0, -1; 0, -1, 0];
                    [x,y,z,a,b,c,l,m,n,bas] = update_(bas, A, tol); T = T * A;
                    continue
                end
                if l * m * n == 1
                    % Procedure A3
                    if l == -1; i = -1; else; i = 1; end
                    if m == -1; j = -1; else; j = 1; end
                    if n == -1; k = -1; else; k = 1; end
                    A = [i, 0, 0; 0, j, 0; 0, 0, k];
                    [x,y,z,a,b,c,l,m,n,bas] = update_(bas, A, tol); T = T * A;
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
                    [x,y,z,a,b,c,l,m,n,bas] = update_(bas, A, tol); T = T * A;
                end
                if ( abs(x) > b + tol || (~ abs(b - x) > tol && 2 * y < z - tol) || (~ abs(b + x) > tol && z < -tol))
                    % Procedure A5
                    A = [1, 0, 0; 0, 1, - sign(x); 0, 0, 1];
                    [x,y,z,a,b,c,l,m,n,bas] = update_(bas, A, tol); T = T * A;
                elseif ( abs(y) > a + tol || (~ abs(a - y) > tol && 2 * x < z - tol) || (~ abs(a + y) > tol && z < -tol))
                    % Procedure A6
                    A = [1, 0, - sign(y); 0, 1, 0; 0, 0, 1];
                    [x,y,z,a,b,c,l,m,n,bas] = update_(bas, A, tol); T = T * A;
                elseif ( abs(z) > a + tol || (~ abs(a - z) > tol && 2 * x < y - tol) || (~ abs(a + z) > tol && y < -tol))
                    % Procedure A7
                    A = [1, - sign(z), 0; 0, 1, 0; 0, 0, 1];
                    [x,y,z,a,b,c,l,m,n,bas] = update_(bas, A, tol); T = T * A;
                elseif (x + y + z + a + b < -tol || (~ abs(x + y + z + a + b) > tol && 2 * (a + y) + z > tol))
                    % Procedure A8
                    A = [1, 0, 1; 0, 1, 1; 0, 0, 1];
                    [x,y,z,a,b,c,l,m,n,bas] = update_(bas, A, tol); T = T * A;
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
        
        function [m]          = get_metric(bas, algo)
            % metric = get_metric(bas)
            % Compute metric tensors from unit cell column basis vector.
            
            if nargin < 2; algo = 2; end
                
            switch algo
                case 1
                    % explicit
                    m = zeros(3,3);
                    for i = 1:3
                    for j = 1:i
                        m(i,j)=dot(bas(:,i),bas(:,j)); m(j,i)=m(i,j);
                    end
                    end
                case 2
                    % matrix multiplication
                    m = bas.'*bas;
            end
        end
        
        
        % brillouin zones

        function [fbz,ibz]    = get_zones(pc,n)

            import am_lib.* am_dft.*

            % get full brillouin zone
            [fbz] = get_fbz(pc,n);

            % get irreducible zone
            [ibz,i2f,f2i] = get_ibz(fbz,pc);

            % save mapping to zones
            fbz.f2i = f2i; fbz.i2f = i2f;
            ibz.i2f = i2f; ibz.f2i = f2i;
        end

        function [bzp]        = get_bz_path(pc,n,brav)

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

            import am_lib.* am_dft.*

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

        function plot_interpolated(fbz,bzp,x, varargin)
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

        function plot_dispersion(model,bzp,flag, varargin)
            % model is either bvk or tb

            import am_lib.* am_dft.*

            % define figure properties
            fig_ = @(h)       set(h,'color','white');
            axs_ = @(h,qt,ql) set(h,'Box','on','XTick',qt,'Xticklabel',ql);
            fig_(gcf);

            switch lower(strtrim(flag))
                case 'electron'
                    % get electron band structure along path
                    bzp = get_dispersion(model,bzp,flag);
                    % plot results
                    plot(bzp.x,sort(bzp.E), varargin{:});
                    axs_(gca,bzp.qt,bzp.ql); axis tight; ylabel('Energy [eV]'); xlabel('Wavevector k');
                case 'phonon'
                    % get phonon band structure along path
                    bzp = get_dispersion(model,bzp,flag);
                    % plot results
                    plot(bzp.x,sort(real(bzp.hw)*1E3),'-k',bzp.x,-sort(abs(imag(bzp.hw))), varargin{:});
                    axs_(gca,bzp.qt,bzp.ql); axis tight; ylabel('Energy [meV]'); xlabel('Wavevector k');
            end


        end

        function plot_dispersion_orbital_character(dft,bzp)
            % FPOSCAR = 'POSCAR'; Ef = 5.0740;
            % [~,pc] = get_cells(FPOSCAR);
            % [dft]   = load_procar('evk/PROCAR',Ef);
            % [bzp]   = get_bz_path(pc,40,'sc');

            import am_lib.*

            % normalize columns of matrix
            normc_ = @(m) ones(size(m,1),1)*sqrt(1./sum(m.*m)).*m;

            % get eigenvalues and band character weights
            c = zeros(dft.nbands,dft.nks); character_labels = {'s','p','d','f'};
            for i = 1:dft.nks
                %    s     py     pz     px    dxy    dyz    dz2    dxz    dx2    f-3    f-2    f-1     f0     f1     f2     f3
                % l-PROJECTION (sum over spin, atoms, m-quantum numbers)
                % lmproj(nspins,norbitals,nions,nbands,nkpts)
                Vp(1,:) = squeeze(sum(sum(sum(dft.lmproj(:,   [1],:,:,i),1),2),3)); % s
                Vp(2,:) = squeeze(sum(sum(sum(dft.lmproj(:, [2:4],:,:,i),1),2),3)); % p
                Vp(3,:) = squeeze(sum(sum(sum(dft.lmproj(:, [5:9],:,:,i),1),2),3)); % d
                Vp(4,:) = squeeze(sum(sum(sum(dft.lmproj(:,[9:16],:,:,i),1),2),3)); % f
                c(:,i) = assign_color_(normc_(Vp));
            end

            % define figure properties
            fig_ = @(h)       set(h,'color','white');
            axs_ = @(h,qt,ql) set(h,'Box','on','XTick',qt,'Xticklabel',ql);
            fig_(gcf);

            % plot band structure
            figure(1); clf; fig_(gcf); hold on;
            for j = 1:dft.nbands; plotc_(bzp.x,dft.E(j,:),c(j,:)); end
            hold off;

            % label axes
            axs_(gca,bzp.qt,bzp.ql); axis tight; ylabel('Energy [eV]'); xlabel('Wavevector k');

            % plot fermi level
            line([0,bzp.x(end)],[0,0],'linewidth',2,'color',[1,1,1]*0.5,'linestyle',':');

            % apply color map and label axes
            colormap( get_colormap('spectral',100).^(2) ); h = colorbar; caxis([0,1]);
            cticks = assign_color_(eye(size(Vp,1))); [~,inds] = sort( cticks );
            set(h,'Ticks',cticks(inds),'TickLabels',character_labels(inds));
        end

        function plot_dispersion_atomic_character(dft,bzp,pc)
            % FPOSCAR = 'POSCAR'; Ef = 5.0740;
            % [~,pc] = get_cells(FPOSCAR);
            % [dft]   = load_procar('evk/PROCAR',Ef);
            % [bzp]   = get_bz_path(pc,40,'sc');

            import am_lib.*

            % normalize columns of matrix
            normc_ = @(m) ones(size(m,1),1)*sqrt(1./sum(m.*m)).*m;

            % get eigenvalues and band character weights
            c = zeros(dft.nbands,dft.nks); character_labels = pc.symb;
            for i = 1:dft.nks
                %    s     py     pz     px    dxy    dyz    dz2    dxz    dx2    f-3    f-2    f-1     f0     f1     f2     f3
                % atomic-PROJECTION: project on irreducible atoms
                % lmproj(nspins,norbitals,nions,nbands,nkpts)
                for j = 1:numel(unique(pc.p2i))
                    Vp(j,:) = squeeze(sum(sum(sum(dft.lmproj(:,:,j==pc.p2i,:,i),1),2),3));
                end
                c(:,i) = assign_color_(normc_(Vp));
            end

            % define figure properties
            fig_ = @(h)       set(h,'color','white');
            axs_ = @(h,qt,ql) set(h,'Box','on','XTick',qt,'Xticklabel',ql);
            fig_(gcf);

            % plot band structure
            figure(1); clf; fig_(gcf); hold on;
            for j = 1:dft.nbands; plotc_(bzp.x,dft.E(j,:),c(j,:)); end
            hold off;

            % label axes
            axs_(gca,bzp.qt,bzp.ql); axis tight; ylabel('Energy [eV]'); xlabel('Wavevector k');

            % plot fermi level
            line([0,bzp.x(end)],[0,0],'linewidth',2,'color',[1,1,1]*0.5,'linestyle',':');

            % apply color map and label axes
            colormap( get_colormap('spectral',100).^(2) ); h = colorbar; caxis([0,1]);
            cticks = assign_color_(eye(size(Vp,1))); [~,inds] = sort( cticks );
            set(h,'Ticks',cticks(inds),'TickLabels',character_labels(inds));
        end

        function plot_nesting(ibz,fbz,bzp,degauss,Ep, varargin)

            import am_lib.*
            import am_dft.*

            % plot results
            plot_interpolated(fbz,bzp, ibz2fbz(fbz,ibz,get_nesting(fbz,ibz,degauss,Ep)) , varargin{:})
        end

        function plot_bz(fbz)

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


        % phonons (harmonic)

        function [bvk,pp]     = get_bvk(pc,uc,md,cutoff2,flags)

            import am_lib.* am_dft.*

            if nargin<5; flags='-identify -model -fit -enforce'; end

            % get irreducible shells
            fprintf(' ... identifying pairs ');
            if contains(flags,'-identify'); tic;
                [ip,pp] = get_pairs(pc,uc,cutoff2);
                fprintf('(%.f secs)\n',toc);
            else; fprintf('(skipped)\n'); end
            
            % [cart] print shell results
            print_pairs(uc,pp)
            
            % [cart] print shell results
            print_ip_symmetries(ip,pp);

            % force constant model
            fprintf(' ... determining harmonic force constant interdependancies and dynamical matrix ');
            if contains(flags,'-model'); tic;
                bvk = get_bvk_model(ip,pp,uc);
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
            
            function print_ip_symmetries(ip,pp)
                bar_ = @(x) repmat('-',[1,x]); fprintf('     %s irreducible pair symmetries %s\n', bar_(25), bar_(25) );
                for i = 1:ip.nshells
                    ex_=ip.s_ck(:,i);
                    ex_bond_group_ = ex_; ex_bond_group_(ex_) = all(pp.Q{2}(:,ex_)==[1;2],1); BG=am_dft.get_long_ss_name(pp.Q{1}(:,:,ex_bond_group_)); % bound group
                    ex_revr_group_ = ex_; ex_revr_group_(ex_) = all(pp.Q{2}(:,ex_)==[2;1],1); RG=am_dft.get_long_ss_name(pp.Q{1}(:,:,ex_revr_group_)); % reversal group
                    fprintf('     irreducible pair %i\n', i);
                    fprintf('       %-10s: ','stabilizer');      
                    for i = 1:numel(BG)
                        fprintf(' %20s',BG{i}); 
                        if mod(i,3)==0 && i~=numel(BG)
                            fprintf('\n'); 
                            fprintf('       %-10s  ',' '); 
                        end
                    end
                    fprintf('\n');
                    fprintf('\n');
                    fprintf('       %-10s: ','reversal');  
                    for i = 1:numel(RG)
                        fprintf(' %20s',RG{i}); 
                        if mod(i,3)==0 && i~=numel(RG)
                            fprintf('\n'); 
                            fprintf('       %-10s  ',' '); 
                        end
                    end
                    fprintf('\n');
                    fprintf('\n');
                end
            end
        end

        function [bvk,D,Dasr] = get_bvk_model(ip,pp,uc)
            % [bvk,D,Dasr] = get_bvk_model(ip,pp,uc)
            %
            % Example for Si:
            %
            %     >> bvk.phi{1} (zeroth neighbor) ==> [0.00000    0.00000    0.00000]
            % 
            %     [ c01_11,      0,      0]
            %     [      0, c01_11,      0]
            %     [      0,      0, c01_11]
            % 
            %     >> bvk.phi{3} (first neighbor)  ==> [1.36750    1.36750    1.36750]
            % 
            %     [ c03_11, c03_21, c03_21]
            %     [ c03_21, c03_11, c03_21]
            %     [ c03_21, c03_21, c03_11]
            % 
            %     >> bvk.phi{2} (second neighbor) ==> [2.73500    2.73500    0.00000]
            % 
            %            (INCORRECT)	        (CORRECT, see MELVIN LAX p 264)
            %     [ c02_11, c02_21,      0]  ==>  [ c02_11, c02_21, c02_31]
            %     [ c02_21, c02_11,      0]  ==>  [ c02_21, c02_11, c02_31]
            %     [      0,      0, c02_33]  ==>  [-c02_31,-c02_31, c02_33]
            %
            %     This erroneous version of the force constant matrix was also
            %     published in 
            %
            %       H. M. J. Smith, Philosophical Transactions of the Royal Society a:
            %       Mathematical, Physical and Engineering Sciences 241, 105 (1948).  
            %
            %     as pointed out in
            %
            %   	M. Lax, Symmetry Principles in Solid State and Molecular Physics
            %   	(Dover Publications, 2012), p 264.
            %
            %     Also correct in 
            %
            %       R. Tubino and J. L. Birman, Phys. Rev. B 15, 5843 (1977).
            %
            %     Interestingly, the fact that second-neighbor force constant is
            %     anti-symmetric implies that the order of differentiation matters and
            %     the intrinsic transpositional symmetry does not apply.
            %
            

            import am_lib.* am_dft.*

            % set sym digits
            digits(10);
                
            % construct transpose super operator (flip symmetry)
            F = zeros(9,9); F(sub2ind([9,9],[1:9],[1,4,7,2,5,8,3,6,9])) = 1;

            % get form of force constants for irreducible prototypical bonds
            for i = 1:ip.nshells
                sym_list = find(ip.s_ck(:,i)); 
                W = kron_( pp.Q{1}(1:3,1:3,sym_list), pp.Q{1}(1:3,1:3,sym_list) );
                for j = 1:numel(sym_list)
                    if all(pp.Q{2}(:,sym_list(j))==[2;1]); W(:,:,j) = F*W(:,:,j); end
                end
                
                % use stabilzer group to determine crystallographic symmetry relations; A*B*C' equals kron(C,A)*B(:)
                W = sum(W-eye(9),3);

                %
                %    This doesn't seem to be a correct symmetry since the example in Lax p. 264
                %    shows a force constant matrix which is anti-symmetric.
                %
                % enforce intrinsic symmetry (immaterial order of differentiation: c == c.')
                % F = zeros(9,9); F(sub2ind([9,9],[1:9],[1,4,7,2,5,8,3,6,9])) = 1; W = W + F-eye(9);

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
                'nbranches',3*pp.pc_natoms,'nshells',size(sav.W,2),'s_ck',ip.s_ck,'W',{sav.W},'phi',{sav.phi},'d',ip.d,'v',ip.v,'xy',ip.xy);
            bvk = bvk_(pp,ip,sav);

            % [cart] define function to get bond vector
            vec_ = @(xy) uc2ws(uc.bas*(uc.tau(:,xy(2,:))-uc.tau(:,xy(1,:))),pp.bas);

            % construct symbolic dynamical matrix
            D=sym(zeros(bvk.nbranches,bvk.nbranches,bvk.nshells)); kvec=sym('k%d',[3,1],'real'); mass=sym('m%d',[1,numel(pp.i2u)],'positive');
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
                D(mp,np,i) = D(mp,np,i) + phi .* exp(sym(2i*pi * rij(:).','d') * kvec(:) );
            end
            end

            % simplify (speeds evaluation up significantly later)
            for i = 1:bvk.nshells; D(:,:,i) = simplify(rewrite(D(:,:,i),'cos'),'steps',20); end

            % stupid matlab, doesn't allow for sum(H,3)
            Dsum=sym(zeros(bvk.nbranches,bvk.nbranches));
            for i = 1:bvk.nshells; Dsum = Dsum + D(:,:,i); end

            % multiply by 1/sqrt(mass)
            mass = repelem(mass(pp.species(pp.p2u)),1,3); mass = 1./sqrt(mass.' * mass); Dsum = Dsum .* mass;

            % attach symbolic dynamical matrix to bvk
            bvk.D = matlabFunction(Dsum);

            % enforce acoustic sum rule algebraically?
            if nargout > 2
                Dasr = sym(zeros(bvk.nbranches,bvk.nbranches));
                for i = 1:bvk.nshells; Dasr = Dasr + subs(subs(subs(D(:,:,i),kvec(1),0),kvec(2),0),kvec(3),0); end
                % simplify and remove TRUE = TRUE (first term) from the set of equations
                Dasr = unique(simplify(Dasr(:)==0)); Dasr = Dasr(2:end);
            end
        end

        function [bvk]        = get_bvk_force_constants(uc,md,bvk,pp,algo)
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

            import am_lib.* am_dft.*
            
            if nargin<5; algo=1; end

            % [cart] get displacements and forces
            u = matmul_( md.bas, mod_( md.tau-uc.tau +.5 )-.5 );
            f = matmul_( md.bas, md.force );

            % select a method (method 1 is the most accurate)
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

            import am_lib.* am_dft.*

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

            import am_lib.* am_dft.*

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

            import am_lib.* am_dft.*

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

        function [tb,H]       = get_tb_model(ip,pp,uc,spdf)
            % set oribtals per irreducible atom: spdf = {'d','p'};
            % may wish to do this to set it per species: x={'p','d'}; spdf={x{ic.species}};
            
            import am_lib.* am_dft.*

            % set sym digits
            digits(10);

            % initialize irreducible atom properties: for each irreducible atom,
            % set azimuthal quantum numbers J{:}, symmetries D{:}, and parity-transpose F{:}
            [J,D,F] = get_tb_symmetry_representation(spdf,pp.Q{1}(1:3,1:3,:));

            % primitive cell atoms define hamiltonian blocks dimensions and start/end sections
            p2i=uc.u2i(uc.p2u); for p=[1:pp.pc_natoms]; d(p)=sum(J{p2i(p)}*2+1); end; E=cumsum(d); S=E-d+1; nbands=E(end);

            % get form of force constants for irreducible prototypical bonds
            for p = 1:ip.nshells
                % get indicies
                x = ip.xy(1,p); i = uc.u2i(x); m = uc.u2p(x); dm = d(m);
                y = ip.xy(2,p); j = uc.u2i(y); n = pp.u2p(y); dn = d(n);

                % have not tested algo 2 yet!!!! try it with Si first...
                algo=2;
                switch algo
                    case 1
                        % original version
                        % use stabilzer group to determine crystallographic symmetry relations; A*B*C' equals kron(C,A)*B(:)
                        W = sum(kron_( D{j}(:,:,ip.s_ck(:,p)) , D{i}(:,:,ip.s_ck(:,p)) ) - eye(dm*dn),3);
                    case 2
                        % incoprorating flip symmetry (need to test still)
                        sym_list = find(ip.s_ck(:,i)); 
                        W = kron_( D{j}(:,:,sym_list) , D{i}(:,:,sym_list) );
                        if (i==j); for j = 1:numel(sym_list) % maybe F{i} should multiply the other right hand side of W?
                            if all(pp.Q{2}(:,sym_list(j))==[2;1]); W(:,:,j) = F{i}*W(:,:,j); end
                        end; end
                        W = sum( W - eye(dm*dn), 3);
                end

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
            H=sym(zeros(tb.nbands,tb.nbands,ip.nshells)); kvec=sym('k%d',[3,1],'real');
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
                H(mp,np,ir) = H(mp,np,ir) + vsk .* exp(sym(2i*pi * rij(:).','d') * kvec(:) );
            end
            end

            % simplify (speeds evaluation up significantly later)
            for i = 1:tb.nshells; H(:,:,i) = simplify(rewrite(H(:,:,i),'cos'),'steps',20); end

            % stupid matlab, doesn't allow for sum(H,3)
            Hsum=sym(zeros(tb.nbands,tb.nbands));
            for i = 1:tb.nshells; Hsum = Hsum + H(:,:,i); end

            % attach symbolic dynamical matrix to bvk
            tb.H = matlabFunction(Hsum);
        end

        function [tb]         = get_tb_matrix_elements(tb,dft,nskips)
            % nskips : number of dft bands to skip (e.g. 5)

            import am_lib.* am_dft.*

            % copy number of bands to skip
            tb.nskips = nskips;

            % fit neighbor parameter at high symmetry points using poor man's simulated anneal
            d4fc = repelem(tb.d,cellfun(@(x)size(x,2),tb.W)); nfcs=numel(d4fc); x=zeros(1,nfcs);
            d=unique(rnd_(d4fc)); d=conv([d,Inf],[1 1]/2,'valid'); nds = numel(d); r_best = Inf;

            % set simulated annealing temeprature and optimization options
            kT = 20; kT_decay_ = @(kT,i) kT .* exp(-i/5); rand_ = @(x) (0.5-rand(size(x))).*abs(x./max(x));
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

        function                plot_tb_vs_dft(tb,dft)

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
            %         UR = reshape(uc2ws_mex(UR,uc.bas,am_dft.tiny),size(UR));

            import am_lib.* am_dft.*

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

                % [pc-frac] define function to find the closest primitive cell vector
                % and shift reference atom to primitive cell 
                G_ = @(tau) round(tau - mod_(tau)); 
                tau = tau-G_(tau(:,:,:,1));
                
                % [uc-frac] record uc index of atoms
                %     relax matching criteria here by a factor of 10;
                %     solves a problem for systems with atoms are at 1/3
                %     position where one of the coordinates may be -0.6666
                %     and the other 0.3334. Applying mod takes -0.6666 to
                %     0.3334 causing a difference of 0.001.
                tau = mod_(matmul_( inv(uc.tau2pc), tau ));
                P1 = member_(tau(:,:,:,1)/10,uc.tau/10);
                P2 = member_(tau(:,:,:,2)/10,uc.tau/10);

                % create a unique pair label (it is very important to request 'stable' here
                % since the position of elements in the PM table also encode information).
                [V,~,V_p2i]=unique([P1(:),P2(:)],'rows','stable'); V=V.';

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
                        % => if this next line fails, check that the first occurance 
                        %    of species m in uc is in the primitive cell (as close to 
                        %    the origin as possible) 
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
            Q{1} = wdv_(Q{1});

            % save "primitive" pairs
            pp = pp_(uc,c_id,o_id,i_id,q_id,iq_id,Q);
        end

        function [it,pt]      = get_triplets(pc,uc,cutoff)

            import am_lib.* am_dft.*

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
                [V,~,V_p2i]=unique([P1(:),P2(:),P3(:)],'rows','stable'); V=V.';

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
            Q{1} = wdv_(Q{1});

            % save "primitive" pairs
            pt = pt_(uc,c_id,o_id,i_id,q_id,iq_id,Q);

        end

        function                print_pairs(uc,pp)

            import am_lib.* am_dft.*

            vec_ = @(xy) uc2ws(uc.bas*(uc.tau(:,xy(2,:))-uc.tau(:,xy(1,:))),uc.bas); Z=[];
            bar_ = @(x) repmat('-',[1,x]); fprintf('     %s primitive shells %s\n', bar_(30), bar_(30) );
            for m = 1:pp.pc_natoms
                Y=[]; ex_ = uniquemask_(pp.i{m});
                fprintf('     atom %i: %i shells\n', m, sum(ex_));
                fprintf('     %-10s    %-30s   %-4s   %-7s   %-7s   %-4s\n', 'd [cart]','bond [cart]','#','ic(i,j)','pc(m,n)','irr.');
                fprintf('     %-10s    %-30s   %-4s   %-7s   %-7s   %-4s\n', bar_(10),bar_(30),bar_(4),bar_(7),bar_(7),bar_(4));
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
                fprintf('     %10.5f  %10.5f %10.5f %10.5f   %4i   %-3i-%3i   %-3i-%3i   %4i\n', Y(:,rankc_(Y(1,:))) ); fprintf('\n');
                Z=[Z,Y];
            end
            w = accumarray(Z(end,:).',Z(5,:).',[],@sum); Z = Z(:,uniquemask_(Z(end,:).')); Z(5,:) = w;
            fprintf('     %s irreducible shells %s\n', bar_(29), bar_(29) );
            fprintf('     %-10s    %-30s   %-4s   %-7s   %-7s   %-4s\n', 'd [cart]','bond [cart]','#','ic(i,j)','pc(m,n)','irr.');
            fprintf('     %-10s    %-30s   %-4s   %-7s   %-7s   %-4s\n', bar_(10),bar_(30),bar_(4),bar_(7),bar_(7),bar_(4));
            fprintf('     %10.5f  %10.5f %10.5f %10.5f   %4i   %-3i-%3i   %-3i-%3i   %4i\n', Z(:,rankc_(Z(1,:))) );
            fprintf('     %s symmetry info %s\n', bar_(29), bar_(29) );

        end

        function                print_triplets(uc,pt)

            import am_lib.* am_dft.*

            vec_ = @(xy) uc2ws(pt.tau(:,xy(2,:)) - pt.tau(:,xy(1,:)),pt.bas); Z=[];
            bar_ = @(x) repmat('-',[1,x]); fprintf('     %s primitive shells %s\n', bar_(30), bar_(30) );
            for m = 1:pt.pc_natoms
                Y=[]; ex_ = uniquemask_(pt.i{m});
                fprintf('     atom %i: %i shells\n', m, sum(ex_));
                fprintf('       %-30s    %-30s   %-4s   %-11s   %-11s   %-4s\n','tau_1 [cart]','tau_2 [cart]','#','ic(i,j,k)','pc(m,n,o)','irr.');
                fprintf('       %-30s    %-30s   %-4s   %-11s   %-11s   %-4s\n',      bar_(30),      bar_(30),bar_(4),bar_(11),bar_(11),bar_(4));
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
                fprintf('     %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f   %4i   %3i-%3i-%3i   %3i-%3i-%3i   %4i\n', Y(2:end,rankc_(Y(1,:)))); fprintf('\n');
                Z=[Z,Y];
            end

            w = accumarray(Z(end,:).',Z(8,:).',[],@sum); Z = Z(:,uniquemask_(Z(end,:).')); Z(8,:) = w;
            fprintf('     %s irreducible shells %s\n', bar_(29), bar_(29) );
            fprintf('       %-30s    %-30s   %-4s   %-11s   %-11s   %-4s\n','tau_1 [cart]','tau_2 [cart]','#','ic(i,j,k)','pc(m,n,o)','irr.');
            fprintf('       %-30s    %-30s   %-4s   %-11s   %-11s   %-4s\n',      bar_(30),      bar_(30),bar_(4),bar_(11),bar_(11),bar_(4));
            fprintf('     %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f   %4i   %3i-%3i-%3i   %3i-%3i-%3i   %4i\n', Z(2:end,rankc_(Z(1,:))));
        end

    end


    % aux library

    methods (Static)

        % aux unit cells

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
            [symb,~,species]=unique(symb,'stable');

            % get space symmetries in conventional setting
            S = generate_sg(sg_code);

            % define function to apply symmetries to position vectors
            seitz_apply_ = @(S,tau) mod_(reshape(matmul_(S(1:3,1:3,:),tau),3,[],size(S,3)) + S(1:3,4,:));

            % apply symmetry operations to all atoms
            natoms = size(tau,2); nSs = size(S,3); i2u = repmat([1:natoms].',1,nSs); i2u=i2u(:).';
            tau = reshape(seitz_apply_(S,tau),3,[]); species = repmat(species(:),1,nSs); species=species(:).';

            % get unique species
            [~,ind] = uniquec_(tau); tau = tau(:,ind); species = species(ind); i2u = i2u(ind);

            % sort by species
            [~,ind] = sort(i2u);  tau = tau(:,ind); species = species(ind); i2u = i2u(ind);

            % define irreducible cell creation function and make structure
            uc_ = @(bas,tau,symb,species) struct('units','frac','bas',bas,...
                'symb',{symb},'mass',get_atomic_mass(get_atomic_number(symb)),...
                'nspecies',sum(unique(species).'==species,2).', ...
                'natoms',size(tau,2),'tau',tau,'species',species);
            uc = uc_(bas,tau,symb,species);
        end

        function [pc,p2u,u2p] = get_primitive_cell(uc, tol)
            % [pc,p2u,u2p] = get_primitive_cell(uc, tol)
            % NOTE: saves p2u entries which share a common closest
            % primitive lattice vector, not just the first primitive atoms
            % produced by the matrix A. When building shells, this property
            % is exploited.

            import am_lib.* am_dft.*

            % set tolerance if not specified
            if nargin<2; tol = am_dft.tiny; end

            % translate one atom to the origin if there isn't one already
            if all(sum(uc.tau,1)>am_dft.tiny); uc.tau = uc.tau-uc.tau(:,1); end

            % build permutation matrix for atoms related by translations
            T = get_symmetries(uc); nTs=size(T,2); PM=zeros(uc.natoms,nTs);
            for i = [1:nTs]; PM(:,i)=rankc_( [mod_(uc.tau(:,:,1)+T(1:3,i));uc.species] ); end

            % construct a sparse binary representation
            A=zeros(uc.natoms); A(sub2ind([1,1]*uc.natoms,repmat([1:uc.natoms].',nTs,1),PM(:)))=1; A=frref_(A); A=A(~all(A==0,2),:);

            % find the smallest primitive cell volume (the three smallest vectors which preserve periodic boundary conditions)
            inds=[0,0,0]; T_cart = uc.bas*T; fwd = rankc_(normc_(T_cart)); T_cart = T_cart(:,fwd); T = T(:,fwd);
            for j =           1:nTs; if sum(abs(T_cart(:,j)))                         >tol; inds(1)=j; break; end; end
            for j = (inds(1)+1):nTs; if sum(abs(cross(T_cart(:,inds(1)),T_cart(:,j))))>tol; inds(2)=j; break; end; end
            for j = (inds(2)+1):nTs; if abs(det(T_cart(:,[inds(1:2),j])+eye(3)*eps))  >tol; inds(3)=j; break; end; end
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
                    M=T_cart(:,ijk(:,i));M=M.'*M.*[1 2 2; 2 1 2; 2 2 1]; n(1,i) = numel(uniquetol(M(:), tol));
                    % for metric tensors with equal symmetries, get the one which has the
                    % most number of angles at 30,60,90,120 deg.
                    n(2,i) = numel(uniquetol([30,60,90,120,bas2abc(T_cart(:,ijk(:,i)))], tol));
                end
                inds_=rankc_(n); B = T(:,ijk(:,inds_(1)));
            end

            % set identifiers (see NOTE: cannot simply use p2u = findrow_(A)!)
            p2u = member_(mod_(B*mod_(B\uc.tau(:,findrow_(A),1))),mod_(uc.tau(:,:,1))).'; u2p = ([1:size(A,1)]*A);

            % define primitive cell creation function and make structure
            pc_ = @(uc,B,p2u) struct('units','frac','bas',uc.bas*B, ...
                'symb',{uc.symb},'mass',uc.mass,'nspecies',sum(unique(uc.species(p2u)).'==uc.species(p2u),2).', ...
                'natoms',numel(p2u),'tau',mod_(B\uc.tau(:,p2u)),'species',uc.species(p2u) );
            pc = pc_(uc,B,p2u);
            
            % set first atom to 0 position, important for symmetries
            pc.tau=mod_(pc.tau-pc.tau(:,1));
        end

        function [ic,i2p,p2i] = get_irreducible_cell(pc, tol)
            % idenitifes irreducible atoms

            import am_lib.* am_dft.*

            % set default numerical tolernece
            if nargin < 2; tol = am_dft.tiny; end
            
            % get seitz matrices
            [~,~,S] = get_symmetries(pc);

            % define function to apply symmetries to position vectors
            seitz_apply_ = @(S,tau) mod_(reshape(matmul_(S(1:3,1:3,:),tau),3,[],size(S,3)) + S(1:3,4,:), tol);

            % get permutation matrix and construct a sparse binary representation
            PM = member_(seitz_apply_(S,pc.tau),pc.tau, tol); A = get_connectivity(PM);

            % set identifiers
            i2p = round(findrow_(A)).'; p2i = round(([1:size(A,1)]*A));

            % define irreducible cell creation function and make structure
            ic_ = @(uc,i2p) struct('units','frac','bas',uc.bas, ...
                'symb',{uc.symb},'mass',uc.mass,'nspecies',sum(unique(uc.species(i2p)).'==uc.species(i2p),2).', ...
                'natoms',numel(i2p),'tau',uc.tau(1:3,i2p),'species',uc.species(i2p));
            ic = ic_(pc,i2p);
        end

        function [xc,x2i,i2x] = get_expanded_cell(ic, S)
            
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
            % define irreducible cell creation function and make structure
            xc_ = @(ic,tau,x2i) struct('units','frac','bas',ic.bas, ...
                'symb',{ic.symb},'mass',ic.mass,'nspecies',sum(unique(ic.species(x2i))==ic.species(x2i).',1), ...
                'natoms',numel(x2i),'tau',tau,'species',ic.species(x2i));
            xc = xc_(ic,tau,x2i);
        end

        function [uc,u2p,p2u] = get_supercell(pc, B)
            % [uc,u2p,p2u] = get_supercell(pc, B)
            % Example: to make a cell pc have the same shape as a reference cell ref, do this:
            % [~,ref] = get_cells('185_P63cm.poscar');
            % [~,pc]  = get_cells('194_P63mmc.poscar');
            % pc.bas = rotzd_(-30)*pc.bas; T = round(pc.bas\ref.bas);
            % sc = get_supercell(pc,T);
            %
            import am_lib.* am_dft.*

            % basic check
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

            % define irreducible cell creation function and make structure
            uc_ = @(uc,tau,B,s2u) struct('units','frac','bas',uc.bas*B,'bas2pc',inv(B),'tau2pc',B,...
                'symb',{{uc.symb{unique(uc.species(s2u))}}},'mass',uc.mass,'nspecies',sum(unique(uc.species(s2u)).'==uc.species(s2u),2).', ...
                'natoms',numel(s2u),'tau',tau,'species',uc.species(s2u));
            uc = uc_(pc,X(2:4,:),B,u2p);

            % % add maps
            % uc.u2p = u2p;
            % uc.p2u = p2u;
            % uc.u2i = pc.p2i(uc.u2p);
            % uc.i2u = uc.p2u(pc.i2p);
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
%                     [uc,pc,ic,cc]=get_cells(fposcar);
%                     [~,~,~,g] = get_symmetries(pc);
                    
                case 1
                    % find the transformation matrix mapping the pc point group to the one in standard setting
                    % find a transformation matrix which takes g to the standard setting h
                    [~,M] = find_pointgroup_transformation(R,generate_pg(identify_pointgroup(R)),1);
                    % get integer inverse matrices
                    C = inv(M); round_(C);
                    
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
                        T(:,i) = round_(T(:,i));
                    end

                    % sort rotation axes:
                    %   1) proper rotation axes first
                    %   1) higher order rotation axes first
                    %   2) vectors in positive octant first
                    cc_bas = pc.bas*T; 
                    fwd = rankc_([-sign(dt);-tr;-max(sign(cc_bas));sign(cc_bas)]);
                    filter_ = @(tr,dt,T,vc_bas,inds,ex_) deal(tr(ex_),dt(ex_),T(:,ex_),vc_bas(:,ex_),inds(ex_));
                    [tr,dt,T,cc_bas,inds] = filter_(tr,dt,T,cc_bas,inds, fwd );

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

        function formula      = get_formula(uc)
            import am_lib.*
            formula = ''; x = uc.nspecies./gcd_(uc.nspecies);
            for j = 1:numel(x)
                formula = [formula,strtrim(uc.symb{j})];
                if x(j)~=1
                    formula=[formula,strtrim(num2str(x(j)))];
                end
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

        function [J,D,F] = get_tb_symmetry_representation(spdf,R)
            % set symmetries D{:}, and parity-transpose F{:} for each
            % irreducible atom given a list of orbitals for each
            % irreducible atom, spdf = {'sp','d'}

            import am_lib.* am_dft.*

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

    end


    % atomic information

    methods (Static)

        function [Z]    = get_atomic_number(symb)
            symbol_database = {...
                 'h'  ,'he' ,'li' ,'be' ,'b'  ,'c'  ,'n'  ,'o'  ,'f'  ,'ne' ,'na' ,'mg' ,'al' ,'si' ,'p'  ,'s'  , ...
                 'cl' ,'ar' ,'k'  ,'ca' ,'sc' ,'ti' ,'v'  ,'cr' ,'mn' ,'fe' ,'co' ,'ni' ,'cu' ,'zn' ,'ga' ,'ge' , ...
                 'as' ,'se' ,'br' ,'kr' ,'rb' ,'sr' ,'y'  ,'zr' ,'nb' ,'mo' ,'tc' ,'ru' ,'rh' ,'pd' ,'ag' ,'cd' , ...
                 'in' ,'sn' ,'sb' ,'te' ,'i'  ,'xe' ,'cs' ,'ba' ,'la' ,'ce' ,'pr' ,'nd' ,'pm' ,'sm' ,'eu' ,'gd' , ...
                 'tb' ,'dy' ,'ho' ,'er' ,'tm' ,'yb' ,'lu' ,'hf' ,'ta' ,'w'  ,'re' ,'os' ,'ir' ,'pt' ,'au' ,'hg' , ...
                 'tl' ,'pb' ,'bi' ,'po' ,'at' ,'rn' ,'fr' ,'ra' ,'ac' ,'th' ,'pa' ,'u'  ,'np' ,'pu' ,'am' ,'cm' , ...
                 'bk' ,'cf' ,'es' ,'fm' ,'md' ,'no' ,'lr' ,'rf' ,'db' ,'sg' ,'bh' ,'hs' ,'mt' ,'ds' ,'rg' ,'uub', ...
                 'uut','uuq','uup','uuh'}; f_ = @(symb) find(strcmp(strtrim(lower(symb)),symbol_database));
            if iscell(symb)
                nZs = numel(symb); Z = zeros(1,nZs);
                for i = 1:nZs; Z(i) = f_(symb{i}); end
            else
                Z = f_(symb);
            end
        end

        function [symb] = get_atomic_symbol(Z)
            %
            symbol_database = {...
                 'h'  ,'he' ,'li' ,'be' ,'b'  ,'c'  ,'n'  ,'o'  ,'f'  ,'ne' ,'na' ,'mg' ,'al' ,'si' ,'p'  ,'s'  , ...
                 'cl' ,'ar' ,'k'  ,'ca' ,'sc' ,'ti' ,'v'  ,'cr' ,'mn' ,'fe' ,'co' ,'ni' ,'cu' ,'zn' ,'ga' ,'ge' , ...
                 'as' ,'se' ,'br' ,'kr' ,'rb' ,'sr' ,'y'  ,'zr' ,'nb' ,'mo' ,'tc' ,'ru' ,'rh' ,'pd' ,'ag' ,'cd' , ...
                 'in' ,'sn' ,'sb' ,'te' ,'i'  ,'xe' ,'cs' ,'ba' ,'la' ,'ce' ,'pr' ,'nd' ,'pm' ,'sm' ,'eu' ,'gd' , ...
                 'tb' ,'dy' ,'ho' ,'er' ,'tm' ,'yb' ,'lu' ,'hf' ,'ta' ,'w'  ,'re' ,'os' ,'ir' ,'pt' ,'au' ,'hg' , ...
                 'tl' ,'pb' ,'bi' ,'po' ,'at' ,'rn' ,'fr' ,'ra' ,'ac' ,'th' ,'pa' ,'u'  ,'np' ,'pu' ,'am' ,'cm' , ...
                 'bk' ,'cf' ,'es' ,'fm' ,'md' ,'no' ,'lr' ,'rf' ,'db' ,'sg' ,'bh' ,'hs' ,'mt' ,'ds' ,'rg' ,'uub', ...
                 'uut','uuq','uup','uuh'};
            nZs = numel(Z); symb = cell(1,nZs);
            for i = 1:nZs; symb(i) = symbol_database(Z(i)); end
        end

        function [mass] = get_atomic_mass(Z)
            mass_database = [...
                    1.007947000,     4.002602000,     6.941200000,     9.012182000,    10.811500000, ...
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
                  263.000000000,   262.000000000,   265.000000000,   266.000000000];
            nZs = numel(Z); mass = zeros(1,nZs);
            for i = 1:nZs; mass(i) = mass_database(Z(i)); end
        end

        function [r]    = get_atomic_radius(Z)
            % [r] = get_atomic_radius(Z)
            % Gives the radii of free atoms [nm].
            %
            % Refs: Vainshtein BK, Fridkin VM, Indenbom VL (1995) Structure of Crystals
            %         (3rd Edition). Springer Verlag, Berlin. 
            %       Clementi E, Raimondi DL, Reinhardt WP (1963). Journal of Chemical
            %         Physics 38:2686
            %    	M. De Graef and M. E. McHenry, Structure of Materials (Cambridge
            %    	  University Press, 2012), p 655.
            database = [...
                 0.053, 0.031, 0.167, 0.112, 0.087, 0.067, 0.056, 0.048, 0.042, 0.038, 0.190, 0.145, 0.118, 0.111, 0.098,  ...
                 0.088, 0.079, 0.071, 0.243, 0.194, 0.184, 0.176, 0.171, 0.166, 0.161, 0.156, 0.152, 0.149, 0.145, 0.142,  ...
                 0.136, 0.125, 0.114, 0.103, 0.094, 0.088, 0.265, 0.219, 0.212, 0.206, 0.198, 0.190, 0.183, 0.178, 0.173,  ...
                 0.169, 0.165, 0.161, 0.156, 0.145, 0.133, 0.123, 0.115, 0.108, 0.298, 0.253, 0.195, 0.185, 0.247, 0.206,  ...
                 0.205, 0.238, 0.231, 0.233, 0.225, 0.228, 0.226, 0.226, 0.222, 0.222, 0.217, 0.208, 0.200, 0.193, 0.188,  ...
                 0.185, 0.180, 0.177, 0.174, 0.171, 0.156, 0.154, 0.143, 0.135, 0.127, 0.120, 0.1  , 0.1  , 0.195, 0.180,  ...
                 0.180, 0.175, 0.175, 0.175, 0.175,     ];
            nZs = numel(Z); r = zeros(1,nZs);
            for i = 1:nZs; r(i) = database(Z(i)); end
        end
        
        function [r]    = get_ionic_radius(Z)
            % [r] = get_ionic_radius(Z)
            % Gives the ionic atomic radis [nm].
            %
            % Refs: Slater JC (1964) Journal of Chemical Physics 39:3199
            % 
            database = [...
                0.025, 0.031, 0.145, 0.105, 0.085, 0.070, 0.065, 0.060, 0.050, 0.038, 0.180, 0.150, 0.125,  ...
                0.110, 0.100, 0.100, 0.100, 0.071, 0.220, 0.180, 0.160, 0.140, 0.135, 0.140, 0.140, 0.140,  ...
                0.135, 0.135, 0.135, 0.135, 0.130, 0.125, 0.115, 0.115, 0.115, 0.088, 0.235, 0.200, 0.185,  ...
                0.155, 0.145, 0.145, 0.135, 0.130, 0.135, 0.140, 0.160, 0.155, 0.155, 0.145, 0.145, 0.140,  ...
                0.140, 0.108, 0.260, 0.215, 0.195, 0.185, 0.185, 0.185, 0.185, 0.185, 0.185, 0.180, 0.175,  ...
                0.175, 0.175, 0.175, 0.175, 0.175, 0.175, 0.155, 0.145, 0.135, 0.135, 0.130, 0.135, 0.135,  ...
                0.135, 0.150, 0.190, 0.180, 0.160, 0.190, 0.127, 0.120, 0.1  , 0.215, 0.195, 0.180, 0.180,  ...
                0.175, 0.175, 0.175, 0.175, 0.1  ];
            nZs = numel(Z); r = zeros(1,nZs);
            for i = 1:nZs; r(i) = database(Z(i)); end
        end
        
        function [r]    = get_crystal_radius(Z)
            % [r] = get_crystal_radius(Z)
            % Gives the radis [nm] of neutral atoms in a closely-packed crystal, since
            % increasing the coordination causes an increase in radius.
            %
            % Refs: M. De Graef and M. E. McHenry, Structure of Materials (Cambridge
            %    	  University Press, 2012), p 655.
            %       Shannon RD Prewitt CT (1969) Acta Crystallographica B25:925-946
            %       Shannon RD (1976) Acta Crystallographica A23:751-761
            database = [...
                0.010, 0.1  , 0.090, 0.041, 0.025, 0.029, 0.030, 0.121, 0.119, 0.1  , 0.116, 0.086, ... 
                0.053, 0.040, 0.031, 0.043, 0.167, 0.1  , 0.152, 0.114, 0.089, 0.075, 0.068, 0.076, ...
                0.081, 0.069, 0.054, 0.070, 0.071, 0.074, 0.076, 0.053, 0.072, 0.056, 0.182, 0.1  , ... 
                0.166, 0.132, 0.104, 0.086, 0.078, 0.079, 0.079, 0.082, 0.081, 0.078, 0.129, 0.092, ... 
                0.094, 0.069, 0.090, 0.111, 0.206, 0.062, 0.181, 0.149, 0.136, 0.115, 0.132, 0.130, ... 
                0.128, 0.110, 0.131, 0.108, 0.118, 0.105, 0.104, 0.103, 0.102, 0.113, 0.100, 0.085, ... 
                0.078, 0.074, 0.077, 0.077, 0.077, 0.074, 0.151, 0.083, 0.103, 0.149, 0.117, 0.108, ... 
                0.076, 0.1  , 0.194, 0.162, 0.126, 0.119, 0.109, 0.087, 0.1  , 0.100, 0.112, 0.111];
            nZs = numel(Z); r = zeros(1,nZs);
            for i = 1:nZs; r(i) = database(Z(i)); end
        end
        
        function [c]    = get_atomic_configuration(Z)
            %
            database = {...
                '1s1                                ', '1s2                                ', '[He] 2s1                           ', '[He] 2s2                           ', ...
                '[He] 2s2 2p1                       ', '[He] 2s2 2p2                       ', '[He] 2s2 2p3                       ', '[He] 2s2 2p4                       ', ...
                '[He] 2s2 2p5                       ', '[He] 2s2 2p6                       ', '[Ne] 3s1                           ', '[Ne] 3s2                           ', ...
                '[Ne] 3s2 3p1                       ', '[Ne] 3s2 3p2                       ', '[Ne] 3s2 3p3                       ', '[Ne] 3s2 3p4                       ', ...
                '[Ne] 3s2 3p5                       ', '[Ne] 3s2 3p6                       ', '[Ar] 4s1                           ', '[Ar] 4s2                           ', ...
                '[Ar] 3d1 4s2                       ', '[Ar] 3d2 4s2                       ', '[Ar] 3d3 4s2                       ', '[Ar] 3d5 4s1                       ', ...
                '[Ar] 3d5 4s2                       ', '[Ar] 3d6 4s2                       ', '[Ar] 3d7 4s2                       ', '[Ar] 3d8 4s2                       ', ...
                '[Ar] 3d10 4s1                      ', '[Ar] 3d10 4s2                      ', '[Ar] 3d10 4s2 4p1                  ', '[Ar] 3d10 4s2 4p2                  ', ...
                '[Ar] 3d10 4s2 4p3                  ', '[Ar] 3d10 4s2 4p4                  ', '[Ar] 3d10 4s2 4p5                  ', '[Ar] 3d10 4s2 4p6                  ', ...
                '[Kr] 5s1                           ', '[Kr] 5s2                           ', '[Kr] 4d1 5s2                       ', '[Kr] 4d2 5s2                       ', ...
                '[Kr] 4d4 5s1                       ', '[Kr] 4d5 5s1                       ', '[Kr] 4d5 5s2                       ', '[Kr] 4d7 5s1                       ', ...
                '[Kr] 4d8 5s1                       ', '[Kr] 4d10                          ', '[Kr] 4d10 5s1                      ', '[Kr] 4d10 5s2                      ', ...
                '[Kr] 4d10 5s2 5p1                  ', '[Kr] 4d10 5s2 5p2                  ', '[Kr] 4d10 5s2 5p3                  ', '[Kr] 4d10 5s2 5p4                  ', ...
                '[Kr] 4d10 5s2 5p5                  ', '[Kr] 4d10 5s2 5p6                  ', '[Xe] 6s1                           ', '[Xe] 6s2                           ', ...
                '[Xe] 5d1 6s2                       ', '[Xe] 4f1 5d1 6s2                   ', '[Xe] 4f3 6s2                       ', '[Xe] 4f4 6s2                       ', ...
                '[Xe] 4f5 6s2                       ', '[Xe] 4f6 6s2                       ', '[Xe] 4f7 6s2                       ', '[Xe] 4f7 5d1 6s2                   ', ...
                '[Xe] 4f9 6s2                       ', '[Xe] 4f10 6s2                      ', '[Xe] 4f11 6s2                      ', '[Xe] 4f12 6s2                      ', ...
                '[Xe] 4f13 6s2                      ', '[Xe] 4f14 6s2                      ', '[Xe] 4f14 5d1 6s2                  ', '[Xe] 4f14 5d2 6s2                  ', ...
                '[Xe] 4f14 5d3 6s2                  ', '[Xe] 4f14 5d4 6s2                  ', '[Xe] 4f14 5d5 6s2                  ', '[Xe] 4f14 5d6 6s2                  ', ...
                '[Xe] 4f14 5d7 6s2                  ', '[Xe] 4f14 5d9 6s1                  ', '[Xe] 4f14 5d10 6s1                 ', '[Xe] 4f14 5d10 6s2                 ', ...
                '[Xe] 4f14 5d10 6s2 6p1             ', '[Xe] 4f14 5d10 6s2 6p2             ', '[Xe] 4f14 5d10 6s2 6p3             ', '[Xe] 4f14 5d10 6s2 6p4             ', ...
                '[Xe] 4f14 5d10 6s2 6p5             ', '[Xe] 4f14 5d10 6s2 6p6             ', '[Rn] 7s1                           ', '[Rn] 7s22                          ', ...
                '[Rn] 6d1 7s2                       ', '[Rn] 6d2 7s2                       ', '[Rn] 5f2 6d1 7s2                   ', '[Rn] 5f3 6d1 7s2                   ', ...
                '[Rn] 5f4 6d1 7s2                   ', '[Rn] 5f6 7s2                       ', '[Rn] 5f7 7s2                       ', '[Rn] 5f7 6d1 7s2                   ', ...
                '[Rn] 5f9 7s2                       ', '[Rn] 5f10 7s2                      ', '[Rn] 5f11 7s2                      ', '[Rn] 5f12 7s2                      ', ...
                '[Rn] 5f13 7s2                      ', '[Rn] 5f14 7s2                      ', '[Rn] 5f14 7s2 7p1                  ', '[Rn] 5f14 6d2 7s2                  ', ...
                '[Rn] 5f14 6d3 7s2                  ', '[Rn] 5f14 6d4 7s2                  ', '[Rn] 5f14 6d5 7s2                  ', '[Rn] 5f14 6d6 7s2                  ', ...
                '[Rn] 5f14 6d7 7s2                  ', '[Rn] 5f14 6d8 7s2                  ', '[Rn] 5f14 6d9 7s2                  ', '[Rn] 5f14 6d10 7s2                 ', ...
                '[Rn] 5f14 6d10 7s2 7p1             ', '[Rn] 5f14 6d10 7s2 7p2             ', '[Rn] 5f14 6d10 7s2 7p3             ', '[Rn] 5f14 6d10 7s2 7p4             ', ...
                '[Rn] 5f14 6d10 7s2 7p5             ', '[Rn] 5f14 6d10 7s2 7p6             ', '[Og] 8s1                           ', '[Og] 8s2                           ', ...
                '[Og] 8s2 8p1                       ', '[Og] 7d1 8s2 8p1                   ', '[Og] 6f1 7d1 8s2 8p1               ', '[Og] 6f3 8s2 8p1                   ', ...
                '[Og] 5g1 6f3 8s2 8p1               ', '[Og] 5g2 6f2 7d1 8s2 8p1           ', '[Og] 5g3 6f2 8s2 8p2               ', '[Og] 5g4 6f2 8s2 8p2               ', ...
                '[Og] 5g5 6f2 8s2 8p2               ', '[Og] 5g6 6f2 8s2 8p2               ', '[Og] 5g7 6f2 8s2 8p2               ', '[Og] 5g8 6f2 8s2 8p2               ', ...
                '[Og] 5g8 6f3 8s2 8p2               ', '[Og] 5g8 6f4 8s2 8p2               ', '[Og] 5g9 6f4 8s2 8p2               ', '[Og] 5g10 6f4 8s2 8p2              ', ...
                '[Og] 5g11 6f3 7d1 8s2 8p2          ', '[Og] 5g12 6f3 7d1 8s2 8p2          ', '[Og] 5g13 6f2 7d2 8s2 8p2          ', '[Og] 5g14 6f3 7d1 8s2 8p2          ', ...
                '[Og] 5g15 6f2 7d2 8s2 8p2          ', '[Og] 5g16 6f2 7d2 8s2 8p2          ', '[Og] 5g17 6f2 7d2 8s2 8p2          ', '[Og] 5g18 6f1 7d3 8s2 8p2          ', ...
                '[Og] 5g18 6f3 7d2 8s2 8p2          ', '[Og] 5g18 6f4 7d2 8s2 8p2          ', '[Og] 5g18 6f5 7d2 8s2 8p2          ', '[Og] 5g18 6f6 7d2 8s2 8p2          ', ...
                '[Og] 5g18 6f6 7d3 8s2 8p2          ', '[Og] 5g18 6f6 7d4 8s2 8p2          ', '[Og] 5g18 6f8 7d3 8s2 8p2          ', '[Og] 5g18 6f9 7d3 8s2 8p2          ', ...
                '[Og] 5g18 6f11 7d2 8s2 8p2         ', '[Og] 5g18 6f12 7d2 8s2 8p2         ', '[Og] 5g18 6f13 7d2 8s2 8p2         ', '[Og] 5g18 6f14 7d2 8s2 8p2         ', ...
                '[Og] 5g18 6f14 7d3 8s2 8p2         ', '[Og] 5g18 6f14 7d4 8s2 8p2         ', '[Og] 5g18 6f14 7d4 8s2 8p2 9s1     ', '[Og] 5g18 6f14 7d5 8s2 8p2 9s1     ', ...
                '[Og] 5g18 6f14 7d6 8s2 8p2 9s1     ', '[Og] 5g18 6f14 7d8 8s2 8p2         ', '[Og] 5g18 6f14 7d9 8s2 8p2         ', '[Og] 5g18 6f14 7d10 8s2 8p2        ', ...
                '[Og] 5g18 6f14 7d10 8s2 8p2 9s1    ', '[Og] 5g18 6f14 7d10 8s2 8p2 9s2    ', '[Og] 5g18 6f14 7d10 8s2 8p2 9s2 9p1', '[Og] 5g18 6f14 7d10 8s2 8p2 9s2 9p2', ...
                '[Og] 5g18 6f14 7d10 8s2 8p3 9s2 9p2', '[Og] 5g18 6f14 7d10 8s2 8p4 9s2 9p2', '[Og] 5g18 6f14 7d10 8s2 8p5 9s2 9p2', '[Og] 5g18 6f14 7d10 8s2 8p6 9s2 9p2'};
            nZs = numel(Z); c = cell(1,nZs);
            for i = 1:nZs; c(i) = strtrim(database(Z(i))); end
        end
        
    end

    
    % unix functions and scripts

    methods (Static)

        % stand-alone

        function generate_scripts()

            import am_lib.* am_dft.*

            % initialize counter
            fprintf('Generating scripts:\n')

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

            fname='queue_vasp.sh'; fid=fopen([fname],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); fprintf(' ... %s (succeeded)\n',fname);
            %{
            #!/bin/bash
            # Antonio Mei Sep/2017
            usage_ () {
                echo "Queues unfinished calculations based on folders containing INCAR,POSCAR,KPOINTS,POTCAR,QSUB."
                echo ""
                echo "Usage: $0 [-h] [-l] [-v] [-n] [-t:target]" 
                echo ""
                echo ""
                echo "-h : prints this message"
                echo "-l : runs on local computer without queuing"
                echo "-t : target file to match"
                echo "-v : shows directories with an incomplete OUTCAR"
                echo "-n : (dry run) prints directors which would get queued"
                echo ""
                exit 1
            }

            main_ () {
                odir=$PWD
                dlist=$(find . -name ${TARGET} -exec sh -c '(cd `dirname {}` && pwd  )' ';' )
                for d in ${dlist}; do
                cd ${d}
                 if ! grep -sq "Total CPU time" OUTCAR; then
                  if [[ -e INCAR   ]]; then
                   if [[ -e POSCAR  ]]; then
                    if [[ -e KPOINTS ]]; then
                     if [[ -e POTCAR  ]]; then
                      echo $PWD
                      if ! ${ISDRY}; then
                          if hash qsub   2>/dev/null; then qsub   QSUB; fi
                          if hash sbatch 2>/dev/null; then sbatch QSUB; fi
                          if local; then vasp; fi
                      fi
                     fi
                    fi
                   fi
                  fi
                 fi
                cd ${odir}
                done
            }

            show_ () {
                grep -LR "Voluntary" --include="OUTCAR" . 
            }

            ISDRY=false; ISLOCAL=false; TARGET='QSUB'
            if (($# == 0)); then usage_; exit 1; fi
            while getopts "nlvht:" o; do
                case "${o}" in
                t)  TARGET=${OPTARG} ;;
                    n)  ISDRY=true ;;
                    l)  ISLOCAL=true ;;
                    v)  show_; exit 0 ;;
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

            import am_lib.* am_dft.*

            % define mex parameters
            FC        = 'ifort';
            FFLAGS    = '-O3 -parallel -fpp -fPIC -lmx -lmex -lmat -nofor_main -bundle -implicitnone -assume realloc_lhs';
            MPATH     = '/Applications/MATLAB_R2016b.app';
            LIBS      = ['-L',MPATH,'/bin/maci64 -I',MPATH,'/extern/include'];
            EXT       = '.mexmaci64';
            DEBUG     = '-debug';

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
                if system([FC,' ',FFLAGS,' ',DEBUG,' ',LIBS,' -c ',flib,'.f90'])
                    error(' ... %s (failed)\n',flib)
                else
                    fprintf(' ... %s (succeeded)\n',flib);
                end

            fprintf('Building mex interfaces:\n')
            for i = 1:numel(f)
                if system([FC,' ',FFLAGS,' ',DEBUG,' ',LIBS,' ',flib,'.o ',f{i},'.f90 -o ',f{i},EXT])
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
