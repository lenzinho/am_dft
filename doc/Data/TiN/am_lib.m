clear;clc;tic

flags='continue';

% get primitive and irreducible cells plus mappings
[uc,pc,ic] = get_cells('POSCAR.TiN.300K',flags);

% get brillouin zones
[fbz,ibz,bzp] = get_zones(pc,[61,61,61],flags);
% 
% % get Born-von Karman phonon model
cutoff=5.5; fname='infile.force_position.TiN.300K_2';
[bvk,ip] = get_bvk(cutoff,pc,uc,fname,flags);
plot_bvk_dispersion(bvk,bzp)
% 
% [ibz] = get_bvk_dispersion(bvk,ibz);

load('matlab.mat'); ninterps=20;
[bzp] = get_bz_path(pc,'fcc');
[bvk] = get_bvk_interpolation(bvk_TiN,bvk_VN,ninterps);

for i = 1:ninterps
    plot_bvk_dispersion(bvk(i),bzp); hold on;
end

% [ibz] = get_bvk_dispersion(bvk,ibz);

% % displace phonon mode 
% k=[0;0.5;.5]; amp=1; mode=6;
% [md] = get_bvk_displacement(bvk,ip,uc,k,amp,mode);

% % run md
% dt = 0.1; nsteps = 2000; Q = 5; T = 300;
% [md] = run_bvk_md(bvk,ip,uc,dt,nsteps,Q,T);
    



% get tight-binding
% cutoff=3; spdf={'d','p'}; nskips=5; fname='EIGENVAL.ibz'; Ef=0;
% [tb,ip] = get_tb(cutoff,pc,spdf,nskips,Ef,fname,flags);

% % plot electron band structure along path
% bzp = get_tb_dispersion(tb,bzp);
% 
% % define figure properties
% fig_ = @(h)       set(h,'color','white');
% axs_ = @(h,qt,ql) set(h,'Box','on','XTick',qt,'Xticklabel',ql);
% 
% figure(1); fig_(gcf); plot(bzp.x,sort(real(bzp.E)),'-k');
% axs_(gca,bzp.qt,bzp.ql); axis tight; ylabel('Energy E'); xlabel('Wavevector k');



fprintf('Done in %i seconds!\n',round(toc));



%%

% vasp

function        save_poscar(uc,fname)
    r_ = @(fid) fprintf(fid,'\n');
    fid=fopen(fname,'w');
    fprintf(fid,'%s ','POSCAR'); r_(fid);
    fprintf(fid,'%12.8f ',uc.latpar); r_(fid);
    fprintf(fid,'%12.8f ',uc.bas(:,1)/uc.latpar); r_(fid);
    fprintf(fid,'%12.8f ',uc.bas(:,2)/uc.latpar); r_(fid);
    fprintf(fid,'%12.8f ',uc.bas(:,3)/uc.latpar); r_(fid);
    fprintf(fid,' %s ',uc.symb{:}); r_(fid);
    fprintf(fid,' %i ',uc.nspecies); r_(fid);
    fprintf(fid,'Direct '); r_(fid);
    fprintf(fid,'%12.8f %12.8f %12.8f \n',uc.tau);
    fclose(fid);
end

function [uc] = load_poscar(fname)
    fid=fopen(fname,'r');                 % open file
    fgetl(fid); uc.units='frac';          % skip header but write units instead
    uc.latpar=sscanf(fgetl(fid),'%f');    % read lattice parameter
    a1=sscanf(fgetl(fid),'%f %f %f');     % first basis vector
    a2=sscanf(fgetl(fid),'%f %f %f');     % second basis vector
    a3=sscanf(fgetl(fid),'%f %f %f');     % third basis vector
    uc.bas=uc.latpar*[a1,a2,a3];          % construct the basis (column vectors)
    uc.recbas=inv(uc.bas);
    uc.vol=abs(det(uc.bas));
    uc.symb=regexp(fgetl(fid), '([^ \s][^\s]*)', 'match');
    uc.nspecies=sscanf(fgetl(fid),repmat('%f' ,1,length(uc.symb)))';
    uc.natoms=sum(uc.nspecies);
    coordtype=lower(strtrim(fgetl(fid)));
    l=0;
    for i=1:length(uc.nspecies)
        for j=1:uc.nspecies(i)
            l=l+1;
            uc.tau(:,l)=sscanf(fgetl(fid),'%f %f %f');
            uc.species(l)=i;
            uc.mass(l)=Z2mass(symb2Z(uc.symb{i}));
        end
    end
    if ~strcmp(coordtype(1),'d'); uc.tau=uc.recbas*uc.tau; end
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

function [dr,bz] = load_vasp_eigenval(fname)
    fprintf('loading dispersion from: %s \n',fname);
    fid=fopen(fname);
    % skip first five lines
    for i = 1:5; fgetl(fid); end
    buffer = strsplit(strtrim(fgetl(fid)));
    dr.nelecs = sscanf(buffer{1},'%i');
    bz.nks  = sscanf(buffer{2},'%i');
    dr.nbands = sscanf(buffer{3},'%i');
    fprintf(' ... electrons = %i \n',dr.nelecs);
    fprintf(' ... kpoints = %i \n',bz.nks);
    fprintf(' ... bands = %i \n',dr.nbands);
    for i = 1:bz.nks
        % skip line
        fgetl(fid);
        % get kpnts
        buffer = strsplit(strtrim(fgetl(fid)));
        bz.k(1,i) = sscanf(buffer{1},'%f');
        bz.k(2,i) = sscanf(buffer{2},'%f');
        bz.k(3,i) = sscanf(buffer{3},'%f');
        % loop over bands
        for j = 1:dr.nbands
            buffer = strsplit(strtrim(fgetl(fid)));
            dr.E(j,i)  = sscanf(buffer{2},'%f');
        end
        dr.E(:,i) = sort(dr.E(:,i));
    end
    fprintf(' ... done\n');
    fclose(fid);
end


% symmetry

function [P] = get_translations(tau,species)
    % define basic functions to enforce numeric precision
    tiny = 1E-12; mod_ = @(x) mod(x+tiny,1)-tiny; rnd_ = @(x) round(x,-log10(tiny));
    % column vector-based rank, sort, unique with numeric precision
    rankcol_ = @(x) sortrowsc(rnd_(x).',[1:size(x,1)]).'; sortcol_ = @(x) sortrows(rnd_(x).').';
    
    % define function to compare first two dimensons
    check3_ = @(A)  all(all(abs(A)<tiny,1),2);
    % get length of each column vector
    normc_ = @(A) sqrt(sum(A.^2,1));

    % sort atoms and species into a unique order (reference)
    X_ = @(tau,species) sortcol_([species;mod_(tau)]); X = X_(tau,species);

    % find vectors that preserve periodic boundary conditions
    N=1; P=mod_(tau(:,species==species(N))-tau(:,N)); nPs=size(P,2); P_ck=false(1,nPs);
    for j = 1:nPs; P_ck(j) = check3_( X_(tau(1:3,:)-P(:,j),species)-X ); end

    % sort P based on lengths
    P=[P(:,P_ck),eye(3)]; P=P(:,rankcol_(normc_(P))); 
end

function [S,R] = get_symmetries(uc)
    % define basic functions to enforce numeric precision
    tiny = 1E-12; mod_ = @(x) mod(x+tiny,1)-tiny; rnd_ = @(x) round(x,-log10(tiny));
    % column vector-based rank, sort, unique with numeric precision
    sortcol_ = @(x) sortrows(rnd_(x).').'; uniquecol_ = @(x) unique(rnd_(x).' ,'rows').'; 
    
    % define outer sum of two vector arrays
    sum_ = @(A,B) reshape(repmat(A,[1,size(B,2)]),size(A,1),[]) + reshape(repmat(B,[size(A,2),1]),size(B,1),[]);    
    % define function to compare first two dimensons
    check3_ = @(A) all(all(abs(A)<tiny,1),2);
    
    % define high dimensional matrix multiplication: A [(m,n),(a,b)] * B [(n,p),(b,c)] = C [(m,p),(a,c)]
    aug_ = @(x,i,y) [x(1:(i-1)),y,x(i:end)]; matmul_ = @(A,B) squeeze(sum(bsxfun(@times, reshape(A,aug_(size(A),3,1)), reshape(B,aug_(size(B),1,1))), 2));

    % define function to check for arithmetic holodries (symmetries for which R'*g*R = g; g = bas'*bas)
    N=9; Q=[-1:1]; nQs=numel(Q);[Y{N:-1:1}]=ndgrid(1:nQs); L=reshape(Q(reshape(cat(N+1,Y{:}),[],N)).',3,3,[]);
    get_rotations_frac_ = @(M) L(:,:,check3_(matmul_(matmul_(permute(L,[2,1,3]),M'*M),L)-M'*M));

    % get seitz operators which leave the atomic basis invariant
    X_= @(tau,species) sortcol_([species;mod_(tau)]); X = X_(uc.tau,uc.species); 
    R = get_rotations_frac_(uc.bas); T = uniquecol_( sum_(uc.tau,-uc.tau) ); 
    nRs = size(R,3); nTs = size(T,2); S = zeros(4,4,nRs*nTs); S(4,4,:)=1; k=0;
    for i = 1:nRs; for j = 1:nTs
        if check3_( X_(R(:,:,i)*uc.tau+T(:,j),uc.species) - X ); k=k+1; S(1:3,1:4,k)=[ R(:,:,i), T(:,j) ]; end
    end; end

    % trim unused space
    S=S(:,:,1:k); R=reshape(uniquecol_( reshape(S(1:3,1:3,:),[9,k]) ),3,3,[]);
end


% unit cells

function [uc,pc,ic]   = get_cells(fname,flags)
    % wrapper function
    % fname = 'infile.supercell' (poscar)
    
    % continue earlier calc?
    sfile = sprintf('%s','am_cells.mat');
    if and(strfind(flags,'continue'),exist(sfile,'file')); load(sfile); return; end 

    % load poscar
    uc = load_poscar(fname);

    % get primitive cell
    [pc,p2u,u2p] = get_primitive_cell(uc);  % write_poscar(get_supercell(pc,eye(3)*2),'POSCAR.2x2')

    % get irreducible cell
    [ic,i2p,p2i] = get_irreducible_cell(pc); 
    
    % complete mapping
    u2i = p2i(u2p); i2u = p2u(i2p);
    
    % save mapping to cells
    pc.p2i = p2i; pc.i2p = i2p;
    pc.p2u = p2u; pc.u2p = u2p;
    
    ic.i2p = i2p; ic.p2i = p2i;
    ic.i2u = i2u; ic.u2i = u2i;
    
    uc.u2p = u2p; uc.p2u = p2u;
    uc.u2i = u2i; uc.i2u = i2u;

    % save
    save('am_cells.mat','uc','pc','ic');
end

function [pc,p2u,u2p] = get_primitive_cell(uc)
    % define basic functions to enforce numeric precision
    tiny = 1E-12; mod_ = @(x) mod(x+tiny,1)-tiny; rnd_ = @(x) round(x,-log10(tiny));
    % column vector-based rank, sort, unique with numeric precision
    rankcol_ = @(x) sortrowsc(rnd_(x).',[1:size(x,1)]).'; 
    
    % define function to find the first nonzero value in each row of matrix A
    findrow_ = @(A) (sum(cumsum(A~=0,2)==0,2)+1).';

    % build permutation matrix for atoms related by the space symmetries
    T=get_translations(uc.tau,uc.species); nTs=size(T,2); PM=zeros(uc.natoms,nTs);
    for i = [1:nTs]; PM(:,i)=rankcol_( [mod_(uc.tau+T(1:3,i));uc.species] ); end

    % construct a sparse binary representation 
    A=zeros(uc.natoms); A(sub2ind([1,1]*uc.natoms,repmat([1:uc.natoms].',nTs,1),PM(:)))=1; A=frref(A); A=A(~all(A==0,2),:);
    
    % set identifiers
    p2u = findrow_(A); u2p = ([1:size(A,1)]*A);

    % set basis (the three smallest vectors which preserve periodic boundary conditions)
    inds=[0,0,0];
    for j = 1:nTs; if any(abs(T(:,j))>tiny); inds(1)=j; break; end; end
    for j = 1:nTs; if any(abs( cross(T(:,2),T(:,j)) )>tiny); inds(2)=j; break; end; end
    for j = 1:nTs; inds(3)=j; if abs(det(T(:,inds))+eye(3)*eps) > tiny; break; end; end
    B=T(:,inds); if det(B)<0; B=fliplr(B); end

    % define primitive cell creation function and make structure
    pc_ = @(uc,B,p2u) struct('units','frac','latpar',uc.latpar,'bas',uc.bas*B,'recbas',uc.recbas/B, ...
        'vol',det(uc.bas*B),'symb',{{uc.symb{uc.species(unique(p2u))}}},'nspecies',sum(unique(p2u).'==p2u,2).', ...
        'natoms',numel(p2u),'tau',mod_(B\uc.tau(:,p2u)),'species',uc.species(p2u),'mass',uc.mass(p2u) );
    pc = pc_(uc,B,p2u);
end

function [ic,i2p,p2i] = get_irreducible_cell(pc)
    % define basic functions to enforce numeric precision
    tiny = 1E-12; mod_ = @(x) mod(x+tiny,1)-tiny; rnd_ = @(x) round(x,-log10(tiny));
    % column vector-based rank, sort, unique with numeric precision
    rankcol_ = @(x) sortrowsc(rnd_(x).',[1:size(x,1)]).'; 
    
    % define function to find the first nonzero value in each row of matrix A
    findrow_ = @(A) (sum(cumsum(A~=0,2)==0,2)+1).';
   
    % get seitz matrices
    S=get_symmetries(pc); nSs=size(S,3);

    % build permutation matrix for atoms related by space symmetries
    PM=zeros(pc.natoms,nSs); for i=[1:nSs]; PM(:,i)=rankcol_( [mod_(S(1:3,1:3,i)*pc.tau+S(1:3,4,i));pc.species] ); end

    % construct a sparse binary representation 
    A=sparse(pc.natoms,pc.natoms); A(sub2ind([1,1]*pc.natoms,repmat([1:pc.natoms].',nSs,1),PM(:)))=1; 
    
    % reduce binary rep using fast sparse rref based on QR decomposition
    A=fsrref(A); A=full(A(any(A~=0,2),:));
    
    % set identifiers
    i2p = findrow_(A); p2i = ([1:size(A,1)]*A);
    
    % define irreducible cell creation function and make structure
    ic_ = @(uc,i2u) struct('units','frac','latpar',uc.latpar,'bas',uc.bas,'recbas',uc.recbas, ...
        'vol',uc.vol,'symb',{{uc.symb{uc.species(i2u)}}},'nspecies',sum(unique(i2u).'==i2u,2).', ...
        'natoms',numel(i2u),'tau',uc.tau(1:3,i2u),'species',uc.species(i2u),'mass',uc.mass(uc.species(i2u))); 
    ic = ic_(pc,i2p);
end

function [sc,s2u,u2s] = get_supercell(uc,B)
    % define basic functions to enforce numeric precision
    tiny = 1E-12; mod_ = @(x) mod(x+tiny,1)-tiny; rnd_ = @(x) round(x,-log10(tiny));
    % column vector-based rank, sort, unique with numeric precision
    uniquecol_ = @(x) unique(rnd_(x).' ,'rows').'; 

    % define outer sum of two arrays of vectors
    sum_ = @(A,B) reshape(repmat(A,[1,size(B,2)]),size(A,1),[]) + reshape(repmat(B,[size(A,2),1]),size(B,1),[]);
    
    % basic check
    if mod(det(B),1)~=0; error('determinant of B must be an integer'); end; 

    % generate primitive lattice vectors
    n=round(sum(abs(B),1)); [Y{3:-1:1}]=ndgrid(1:n(1),1:n(2),1:n(3)); nLs=prod(n); L=reshape(cat(3+1,Y{:})-1,[],3).'; 
    
    % expand atoms, coordinates supercell fractional, and reduce to primitive supercell
    X=uniquecol_([ reshape(repmat([1:uc.natoms],nLs,1),1,[]); mod_(inv(B)*sum_(L,uc.tau)) ]);
    
    % create mapping
    s2u = X(1,:); [~,u2s]=unique(s2u); u2s=u2s(:).';

    % define irreducible cell creation function and make structure
    sc_ = @(uc,tau,B,s2u) struct('units','frac','latpar',uc.latpar,'bas',uc.bas*B,'recbas',uc.recbas/B, ...
        'vol',det(uc.bas*B),'symb',{{uc.symb{uc.species(unique(s2u))}}},'nspecies',sum(unique(s2u).'==s2u,2).', ...
        'natoms',numel(s2u),'tau',tau,'species',uc.species(s2u),'mass',uc.mass(s2u));
    sc = sc_(uc,X(2:4,:),B,s2u);
end

function [ip] = get_irreducible_shells(pc,cutoff)
    % Cutoff needs to be small otherwise shells at large distances will be incomplete.

    % define basic functions to enforce numeric precision
    tiny = 1E-12; mod_ = @(x) mod(x+tiny,1)-tiny; rnd_ = @(x) round(x,-log10(tiny));
    % column vector-based rank, sort, unique with numeric precision
    rankcol_ = @(x) sortrowsc(rnd_(x).',[1:size(x,1)]).'; 
    sortcol_ = @(x) sortrows(rnd_(x).').'; uniquecol_ = @(x) unique(rnd_(x).' ,'rows').'; 
    
    % define high dimensional matrix multiplication: A [(m,n),(a,b)] * B [(n,p),(b,c)] = C [(m,p),(a,c)]
    aug_ = @(x,i,y) [x(1:(i-1)),y,x(i:end)]; matmul_ = @(A,B) squeeze(sum(bsxfun(@times, reshape(A,aug_(size(A),3,1)), reshape(B,aug_(size(B),1,1))), 2));
    % define function to find the first nonzero value in each row of matrix A
    findrow_ = @(A) (sum(cumsum(A~=0,2)==0,2)+1).';
    % get column normals
    normc_ = @(A) sqrt(sum(A.^2,1));
    % check second dimension
    check2_ = @(A) all(abs(A)<tiny,1);
    % define function to get check if rotation is a site-stabilizers (R*tau = tau)
    is_site_stabilizer_ = @(R,tau) all(all(abs(sortcol_(mod_(R*tau))-sortcol_(mod_(tau)))<tiny,1),2);
    % pad with zeros
    pad_ = @(A,n) [A,zeros(1,n(2)-size(A,2));zeros(n(1)-size(A,1),n(2))];

    % get point symmetries in [frac] and [cart]
    [~,R]=get_symmetries(pc); Rc=matmul_(matmul_(pc.bas,R),pc.recbas); nRs=size(R,3);
    
    % get supercell (TO DO: NEED AN INTELIGENT WAY OF SETTING B BASED ON CUTOFF)
    B=eye(3)*5; [sc,s2p,p2s] = get_supercell(pc,B);
    
    % determine primitive atoms shells
    Z=[]; Y_ = @(d,r,w,species,pc_id,ic_id,Rstb_ck,Rgen_ck) [d(:);r(:);w(:);species(:);pc_id(:);ic_id(:);pad_(Rstb_ck(:),[nRs,1]);pad_(Rgen_ck(:),[nRs,1])];
    for uc_id = p2s
        % get n-centered atomic basis
        T = sc.tau-sc.tau(:,uc_id); 

        % get site stabilizer checklist 
        R_stb_ck=false(1,nRs); for i=[1:nRs]; R_stb_ck(i)=is_site_stabilizer_(R(:,:,i),T); end; R_stb_id=find(R_stb_ck);

        % find which atoms are related by each symmetry operation (O_ck) by build and solving permutation matrix for orbits
        PM=zeros(sc.natoms,nRs); for j=[1:nRs]; PM(:,j)=rankcol_([mod_(R(:,:,j)*T);sc.species]); end
        A=sparse(sc.natoms,sc.natoms); A(sub2ind( size(A),repmat([1:sc.natoms]',nRs,1),PM(:) ))=1; A=fsrref(A); O_ck=full(A(any(A~=0,2),:));

        % identify reps and set number of orbits
        O_rep_id=findrow_(O_ck); norbits=numel(O_rep_id);

        % conver to [cart] and translate all orbit reps to wigner-seitz cell; then compute shell radii
        O_rep = uc2ws(sc.bas*T(:,O_rep_id),sc.bas); dc=normc_(O_rep);
        
        % goal of the loop over shells:
        % 	1) flip the bond to sort primitive atom indices
        %   2) save [s_ck] the point group symmetries which stabilize the bond 
        %   3) save [g_ck] the point group symmetries which generate the orbit
        Y=zeros(11+nRs*2,norbits);
        for i = 1:norbits; if dc(i)<cutoff

            % set the rep again, except more systematically this time (Z axis sorted preferentially)
            O_full = flipud(uniquecol_( flipud(matmul_(Rc(:,:,R_stb_id),O_rep(:,i))) )); w = size(O_full,2);
            
            % flip to order unit cell indicies
            if s2p(uc_id)>s2p(O_rep_id(i))
                mn=[O_rep_id(i), uc_id]; O_rep(:,i)=-O_full(:,1);
            else
                mn=[uc_id, O_rep_id(i)]; O_rep(:,i)=+O_full(:,w);
            end
            
            % now do it again with the new rep to get the ordering of unique correct 
            O_full = matmul_(Rc(:,:,R_stb_id),O_rep(:,i));

            % create orbital generator checklist
            [~,u]=unique( rnd_(flipud(O_full)).' ,'rows'); Rgen_ck=zeros(1,nRs); Rgen_ck(R_stb_id(u))=1;

            % save results
            Y(:,i) = Y_(dc(i),O_rep(:,i),w,sc.species(mn),s2p(mn),pc.p2i(s2p(mn)),check2_(matmul_(Rc,O_rep(:,i))-O_rep(:,i)),Rgen_ck);
        end; end; 
    
        % append Z
        Z = [Z,Y(:,dc<cutoff)];
    end
    
    % reduce to unique rows
    [~,ip_id]=unique( rnd_(Z([1:7,10:11],:)).','rows');
    
    % save results Y([1,[2,3,4],5,[6,7],[8,9],[10,11]],nshells)
    pair_ = @(pc,R,Z) struct('units','cart','bas',pc.bas,'recbas',pc.recbas, ...
        'symb',{{pc.symb{:}}},'species',pc.species,'mass',pc.mass,...
        'nshells',size(Z,2),'norbits',Z(5,:),'v',Z([2,3,4],:),'m',Z(8,:),'n',Z(9,:),'i',Z(10,:),'j',Z(11,:), ...
        'nRs',size(R,3),'R',R,'s_ck',logical(Z(11+[1:nRs],:)),'g_ck',logical(Z(11+nRs+[1:nRs],:)));
    ip = pair_(pc,Rc, Z(:,ip_id));
    
    % print table
    bar_ = @(x) repmat('-',[1,x]);
    fprintf('%-10s\n', 'irreducible shells');
    fprintf('%-10s    %-30s    %-4s    %-7s    %-7s    %-7s\n', 'd [cart]','bond [cart]','#','species','pc(m,n)','ic(i,j)'); 
    fprintf('%-10s    %-30s    %-4s    %-7s    %-7s    %-7s\n', bar_(10),bar_(30),bar_(4),bar_(7),bar_(7),bar_(7));
    fprintf('%10.5f  %10.5f %10.5f %10.5f    %4i    %-3i-%3i    %-3i-%3i    %-3i-%3i\n', Z(1:11,ip_id) );
end


% brillouin zones

function [fbz,ibz,bzp] = get_zones(pc,n,flags)

    % continue earlier calc?
    sfile = sprintf('%s','am_zones.mat');
    if and(strfind(flags,'continue'),exist(sfile,'file')); load(sfile); return; end 

    % get full brillouin zone
    [fbz] = get_fbz(pc,n);

    % get irreducible zone
    [ibz,i2f,f2i] = get_ibz(fbz,pc);

    % get zone path
    [bzp] = get_bz_path(pc,'fcc');

    % save mapping to zones
    fbz.f2i = f2i; fbz.i2f = i2f;
    ibz.i2f = i2f; ibz.f2i = f2i;
    
    % save
    save(sfile,'fbz','ibz','bzp','i2f','f2i');
end

function [bzp] = get_bz_path(pc,brav)

    switch strtrim(lower(brav))
        case 'fcc'
            % define kpoint path
            G=[0;0;0]; X1=[0;1;1]/2; X2=[2;1;1]/2; L=[1;1;1]/2; K=[6;3;3]/8;
            N = 50; ql={'G','X','K','G','L'}; qs=[G,X2,K,G]; qe=[X1,K,G,L]; nqs=size(qs,2); 
        otherwise 
            error('invalid bravais lattice');
    end

    % get path: convert to [cart-recp] to get x spacings right then convert back to [frac-recp]
    [k,x,qt] = get_path(pc.recbas*qs,pc.recbas*qe,nqs,N); k=pc.bas*k;
    
    % create path object
    bzp_ = @(pc,ql,qt,nks,x,k) struct('units','frac-recp','latpar',pc.latpar,'bas',pc.bas,'recbas',pc.recbas,'vol',pc.vol,'ql',{{ql{:}}},'qt',qt,'nks',nks,'x',x,'k',k);
    bzp = bzp_(pc,ql,qt,size(k,2),x,k);

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

function [fbz] = get_fbz(pc,n)

    % check
    if any(mod(n,1)~=0); error('n must have integers'); end
    
    % generate primitive lattice vectors
    Q_ = @(i) [0:(n(i)-1)]/n(i); [Y{1:3}]=ndgrid(Q_(1),Q_(2),Q_(3)); k=reshape(cat(3+1,Y{:}),[],3).';
    
    % define irreducible cell creation function and make structure
    fbz_ = @(uc,n,k) struct('units','frac-recp','latpar',uc.latpar,'bas',uc.bas,'recbas',uc.recbas,'vol',uc.vol,...
        'n',n,'nks',size(k,2),'k',k,'w',ones([1,size(k,2)]));
    fbz = fbz_(pc,n,k);

end

function [ibz,i2f,f2i] = get_ibz(fbz,pc)
    % define basic functions to enforce numeric precision
    tiny = 1E-12; mod_ = @(x) mod(x+tiny,1)-tiny; rnd_ = @(x) round(x,-log10(tiny));
    % column-based rank, sort, find, unique
    rankcol_ = @(x) sortrowsc(rnd_(x).',[1:size(x,1)]).'; 
    % matmul: A [(m,n),(a,b)] * B [(n,p),(b,c)] = C [(m,p),(a,c)]
    aug_ = @(x,i,y) [x(1:(i-1)),y,x(i:end)]; matmul_ = @(A,B) squeeze(sum(bsxfun(@times, reshape(A,aug_(size(A),3,1)), reshape(B,aug_(size(B),1,1))), 2));
    
    % get point symmetries [real-frac --> rec-frac] by applying basis transformation twice
    [~,R] = get_symmetries(pc); nRs = size(R,3); R=matmul_(pc.bas^2,matmul_(R,pc.recbas^2)); 

    % build permutation matrix for kpoints related by point symmetries
    PM=zeros(fbz.nks,nRs); for i=[1:nRs]; PM(:,i)=rankcol_(mod_(R(:,:,i)*fbz.k)); end

    % get irreducible rows (requires complete symmetry group and complete k-point basis!)
    [A,i2f,f2i]=unique(sort(PM,2),'rows'); i2f=i2f(:).'; f2i=f2i(:).'; A([false(size(A(:,1))),diff(A,1,2)==0])=0; w=sum(A~=0,2).';
    if abs(sum(w)-prod(fbz.n))>tiny; error('mismatch: kpoint mesh and point group symmetry'); end 
    
    % get stabilizer and star generator checklist
    % s_ck = [PM(i2f,:)==i2f.'].'; % g_ck  = first occurance of each number in each row... to do later...
    
    % get irreducible tetrahedra
    [tet,~,tet_f2i] = unique(sort(f2i(get_tetrahedra(fbz.recbas,fbz.n))).','rows'); tet=tet.'; tetw = hist(tet_f2i,[1:size(tet,2)].'-.5);

    % define irreducible cell creation function and make structure
    ibz_ = @(fbz,uc,i2f,w,R,tet,tetw) struct('units','frac-recp','latpar',uc.latpar,'bas',uc.bas,'recbas',uc.recbas,'vol',uc.vol,...
    	'n',fbz.n,'nks',numel(i2f),'k',fbz.k(:,i2f),'w',w,'ntets',size(tet,2),'tet',tet,'tetw',tetw,...
        'nRs',size(R,3),'R',R);
    ibz = ibz_(fbz,pc,i2f,w,R,tet,tetw);
end


% phonons

function [bvk,ip] = get_bvk(cutoff,pc,uc,fname,flags)
    % for paper:
    % cutoff = 5; % Angstroms 
    % fname = 'infile.force_position.4.00-300'
    
    % continue earlier calc?
    sfile = sprintf('am_bvk_%s.mat',fname);
    if and(strfind(flags,'continue'),exist(sfile,'file')); load(sfile); return; end 

    % get irreducible shells
    ip  = get_irreducible_shells(pc,cutoff);

    % force force constant model
    bvk = get_bvk_model(ip);

    % get force constants
    bvk = get_bvk_force_constants(bvk,ip,uc,fname);

    % save everything generated
    save(sfile,'bvk','ip')
end

function [bvk] = get_bvk_model(ip)
        %
        % NOTE #1: sum rules are enforced again in the extraction of
        % symmetry-adapted force constants and the implementation here is
        % redundant. Still, keep in mind that, it will add more parameters to
        % the evaluation of the dynamical matrix.
        %
        % NOTE #2: kpoints are in fractional reciprocal units but force
        % constants and atomic positions are in cartesian units.
        %

        % define basic functions to enforce numeric precision
        tiny = 1E-12; rnd_ = @(x) round(x,-log10(tiny)); normc_ = @(A) sqrt(sum(A.^2,1)); findrow_ = @(A) (sum(cumsum(A~=0,2)==0,2)+1).';

        % define high dimensional kronecker self-product: C(:,:,k) = kron(A(:,:,i),B(:,:,j))
        kron_ = @(A,B) reshape(bsxfun(@times, permute(A,[4 1 5 2 3]), permute(B,[1 4 2 5 3])), size(A,1)*size(B,1),[],size(A,3));

        % initialize dynamical matrix and kvector
        nbands=3*numel(unique(ip.m)); D=sym(zeros(nbands)); kvec=sym('k%d',[3,1],'real'); mass=sym('m%d',[1,numel(unique([ip.j,ip.j]))]);
        [~,~,shell]=unique(rnd_(normc_(ip.v))); shell=shell(:).'; alloc_ = @(cel) struct('W',cel,'c',cel,'shell',cel); sav = alloc_({cell(1,ip.nshells)});
        for p = 1:ip.nshells
            % use stabilzer group to determine crystallographic symmetry relations; A*B*C' equals kron(C,A)*B(:)
            W = sym(sum(  kron_(ip.R(:,:,ip.s_ck(:,p)),ip.R(:,:,ip.s_ck(:,p)))  -eye(9),3));

            % enforce intrinsic symmetry (immaterial order of differentiation: c == c.')
            F=zeros(9,9); F(sub2ind([9,9],[1:9],[1,4,7,2,5,8,3,6,9]))=1; W = W + F - eye(9);

            % get nullspace and normalize to first nonzero element
            W = null(W); for i = 1:size(W,2); W(:,i) = W(:,i)/W(find(W(:,i),1),i); end; W = rref(W.').';

            % define parameters and save it for ASR later if it's an self-force 
            c = sym(sprintf('c%02i_%%d%%d',p),[3,3],'real'); c = c(findrow_(double(W).'));

            % print symmetry adapted force constants
            fprintf('v = [%10f,%10f,%10f] \n', ip.v(:,p)); phi = reshape( W*c.',[3,3]); disp( phi )

            % get block of dynamical matrix (see NOTE #2)
            B = zeros(3,3); 
            for i = find(ip.g_ck(:,p)).'
                B = B + ip.R(:,:,i) * phi * ip.R(:,:,i)' ...
                    .* exp( sym(2i*pi) * dot(rnd_(ip.recbas*ip.R(:,:,i)*ip.v(:,p)),kvec) ); 
            end

            % update block matrix with masses
            B = B/sqrt( mass(ip.i(p))*mass(ip.j(p)) );

            % augment dynamical matrix
            mp = 3*(ip.m(p)-1)+[1:3]; np = 3*(ip.n(p)-1)+[1:3];
            D( mp,np ) = D( mp,np ) + B;  if (ip.m(p)~=ip.n(p)); D( np,mp ) = D( np,mp ) + B'; end

            % save important stuff (sort W to be in line with c, matlabFunction sorts D variables)
            [sav.c{p},n] = sort(c); sav.W{p} = W(:,n); sav.shell{p} = ones(size(c))*shell(p);
        end

        % simplify the dynamical matrix
        D = simplify(D,'Seconds',500);

        % get acoustic sum rules (see Wooten p 646, very illuminating H2O example)
        D_gamma = subs(subs(D,kvec,[0;0;0]), mass,ones(1,numel(mass))); nprims = numel(unique([ip.m,ip.n])); asr = sym(zeros(3*nprims,3));
        for m = 1:nprims; for n = 1:nprims
            mp = 3*(m-1)+[1:3]; np = 3*(n-1)+[1:3]; asr(mp,1:3) = asr(mp,1:3) + D_gamma(mp,np);
        end; end; asr = unique(simplify(asr==0)); asr = asr(asr~=true);

        % impose acoustic sum rules on D (see NOTE #1)
        c = [sav.c{:}];
        for j = 1:size(asr,1)
            for i = find([sav.shell{:}]==1)
                sol = solve(asr(j),c(i));
                if ~isempty(sol); D=subs(D,c(i),sol); end
            end
        end

        % create bvk structure
        bvk_ = @(ip,sav,D) struct('units','cart,frac-recp','bas',ip.bas,'recbas',ip.recbas,'natoms',numel(unique(ip.m)),'mass',ip.mass, ...
            'nshells',size(sav.W,2),'W',{sav.W},'shell',{sav.shell},'nbands',nbands,'D',matlabFunction(D));
    bvk = bvk_(ip,sav,D);
end

function [bvk] = get_bvk_force_constants(bvk,ip,uc,fname_force_position)
    % Extracts symmetry adapted force constants.
    % get_bvk_force_constants(bvk,ip,load_poscar('infile.supercell'),'infile.force_position.4.04-300')
    % reference supercell : [sc] = load_poscar('infile.supercell');
    %
    % NOTE #1: matching to an irreducible atom may not work if this
    % instance of the irreducible atom has an orbit which is rotated, i.e.
    % this irreducible atom is related through a space symmetry which is
    % not pure translational. I will fix it when this case arises.

    % define basic functions to enforce numeric precision
    tiny = 1E-12; mod_ = @(x) mod(x+tiny,1)-tiny;

    % define high dimensional matrix multiplication: A [(m,n),(a,b)] * B [(n,p),(b,c)] = C [(m,p),(a,c)]
    aug_ = @(x,i,y) [x(1:(i-1)),y,x(i:end)]; matmul_ = @(A,B) squeeze(sum(bsxfun(@times, reshape(A,aug_(size(A),3,1)), reshape(B,aug_(size(B),1,1))), 2));

    % get central atom c_id and orbit o_id supercell indices
    [c_id,o_id,k2ijor] = get_bvk_lookup_tables(ip,uc);

    % load force positions md( 1:3=positions, 4:6=forces , natoms , nsteps )
    fprintf(' ... loading force vs positions\n');
    tic; [md] = load_force_position(uc.natoms,fname_force_position); nsteps = size(md,3); toc

    % convert atomic coordinates [cart,ang] into displacements [cart,ang] by first converting to frac [unitless] and taking modulo
    md(1:3,:,:) = matmul_( uc.bas, mod_(matmul_(uc.recbas,md(1:3,:,:))-uc.tau+0.5)-0.5 );

    % loop over irreducible types
    for i = unique(ip.i)
        % get number of central atoms (stars with orbits around them)
        nstars = numel(c_id{i});

        % get forces : f = [ (x,y,z), (1:natoms)*nsteps ] 
        % get displacements : u [ (x,y,z)*orbits, (1:natoms)*nsteps ]
        %  ... and solve for the generalized force constants: FC = f / u 
        fc = - reshape(md(4:6,c_id{i},:),3,nstars*nsteps) / reshape(md(1:3,o_id{i},:),3*size(o_id{i},1),nstars*nsteps);

        % rotate to proper irreducible bond orientation, but before enforce ASR (correct the self force-constant matrix)
        fc3x3 = reshape(fc,3,3,[]); fc3x3(:,:,1) = -sum(fc3x3(:,:,2:end),3);
        for k = 1:size(fc3x3,3)
            fc3x3(:,:,k) = ip.R(:,:,k2ijor{i}(k,4))' * fc3x3(:,:,k) * ip.R(:,:,k2ijor{i}(k,4));
        end

        % solve for the symmetry adapted force constants : A x = B
        for j = 1:ip.nshells
            mask = k2ijor{i}(:,2)==j; sum_mask = sum(mask);
            if sum_mask>0
                A = repmat(double(bvk.W{j}),sum_mask,1);
                B = reshape(fc3x3(:,:,mask),[],1);
                % get force constants as row vectors
                bvk.fc{j} = reshape( A \ B , 1, []); 
            end
        end

        % print results
        textwidth=64;
        bar_    = @(x) repmat('-',[1,x]);
        strpad_ = @(A,n) [A,repmat(' ',[1,n-length(A)])]; 
        fprintc_= @(t,w) fprintf([strjust(strpad_(t,w),'center'),'\n']);

        % header
        fprintc_(sprintf('irreducible atom %i',i),textwidth);
        fprintc_(bar_(textwidth),textwidth);

        % print table : sum of generalized force constants (should be close to zero)
        fprintf('%-32s\n', ' generalized force constant sum'); 
        fprintf('%-32s\n', bar_(32));
        fprintf('%10.5f %10.5f %10.5f\n', sum(reshape(fc,3,3,[]),3).' );

        % print table : rotated generalized and symmetry adapted force constants
        fprintf('%-32s        %-32s    (%5s, %5s, %5s)\n','         generalized fc','          sym-adapted fc','atom','shell','orbit');
        for k = 1:size(fc3x3,3)
            fprintf('%-32s        %-32s    (%5i, %5i, %5i)\n',bar_(32),bar_(32),k2ijor{i}(k,1:3)); 
            fprintf('%10.5f %10.5f %10.5f        %10.5f %10.5f %10.5f\n', [fc3x3(:,:,k), reshape(  double(bvk.W{k2ijor{i}(k,2)}) * bvk.fc{k2ijor{i}(k,2)}.'  ,3,3) ].');
        end
    end

    function [md] = load_force_position(natoms,fname)
        %
        % First, preprocess the outcar using the bash function below. Then,
        % call load_force_position on the file produced.
        %
        % #!/bin/bash
        % # Preprocess outcar to remove the last run, which may not have finished
        % # Antonio Mei Nov/2014
        % # Antonio Mei Jan/2017
        %
        % usage_ () {
        %     echo "Creates infile.force_position based on supplied outcar files."
        %     echo ""
        %     echo "Usage: $0 [-h] [-t] [-f] -o <outcar_list> -n <natoms> [-c <compress_name>]"
        %     echo ""
        %     echo "Example: $0 -f -t -o \"\$(find . -name OUTCAR | grep 4.00-300)\" -n 250"
        %     echo ""
        %     echo "-h : prints this message"
        %     echo "-n : [REQUIRED] number of atoms in the simulation cell"
        %     echo "-o : [REQUIRED] list of outcar files to parse"
        %     echo "-t : trims the last md run (useful for removing runs which have not completed)"
        %     echo "-f : overwrites existing infile.force_position"
        %     echo "-c : compresses infile.force_position to a tar.gz file"
        %     echo ""
        %     echo "infile.force_position file contents:"
        %     echo "   x position   y position   z position     x force      y force      z force"
        %     exit 1
        % }
        %
        % main_ () {
        %     # trim the last md run which may not have completed
        %     trim_ () { tac $1 | awk '!found && /POSITION/{found=1;next}1' | tac ; }
        %     # get position and forces
        %     get_  () { cat $2 | grep -h -A $(($1+1)) POSITION  ; }
        %     # cut header lines
        %     cut_  () { cat $1 | sed '/^--$/d' | sed '/--/d' | sed '/POSITION/d' ; }
        %     # compress produced infile.force_position
        %     compress_ () { tar -zcvf infile.force_position.tar.gz infile.force_position ; }
        %     #
        %     if ${ISFORCE}; then
        %         if [ -f "./infile.force_position" ]; then
        %             rm ./infile.force_position
        %             printf " ... ./infile.force_position overwritten\n"
        %         fi
        %     fi
        %     #
        %     if ${ISTRIM}; then
        %         printf " ... trim:\n"
        %         for F in "${FLIST}"; do
        %             printf " ...     %-100s\n" "${F}"
        %             trim_ ${F} | get_ ${NATOMS} | cut_ >> infile.force_position
        %         done
        %     else
        %         printf " ... batch parsing without trim\n"
        %         get_ ${NATOMS} "${FLIST}" | cut_ >> infile.force_position
        %     fi
        %     #
        %     printf " ... infile.force_position created\n"
        %     #
        %     if ${ISCOMPRESS}; then
        %         printf " ... infile.force_position.tar.gz compressed\n"
        %         compress_ 
        %     fi
        % }
        %
        % ISCOMPRESS=false; ISTRIM=false; ISFORCE=false;
        % if (($# == 0)); then usage_; exit 1; fi
        % while getopts "n:o:htfc" o; do
        %     case "${o}" in
        %         o)  FLIST=${OPTARG} ;;
        %         n)  NATOMS=${OPTARG} ;;
        %         c)  ISCOMPRESS=true ;;
        %         t)  ISTRIM=true ;;
        %         f)  ISFORCE=true ;;
        %         h)  usage_; exit 0 ;;
        %         *)  usage_; exit 1 ;;
        %     esac
        % done
        %
        % main_
        %

        % count number of lines in file and check that all runs completed properly
        nlines = count_lines(fname); if mod(nlines,natoms)~=0; error('lines appear to be missing.'); end;

        % open file and parse
        fid = fopen(fname); md = reshape(fscanf(fid,'%f'),6,natoms,nlines/natoms); fclose(fid);

        function [nlines] = count_lines(fname)
            if ispc
                [~,a] = system(sprintf('type %s | find /c /v ""',fname));
            elseif or(ismac,isunix)
                [~,a] = system(sprintf('wc -l %s',fname));
            end
            nlines = sscanf(a,'%i');
        end
    end
end

function [bvk] = get_bvk_interpolation(bvk_1,bvk_2,n)
    % interpolates froce constants and masses from bvk_1 and bvk_2 on n points (includes end points)

    bvk_ = @(bvk,mass,fc) struct('units','cart,frac-recp','bas',bvk.bas,'recbas',bvk.recbas,'natoms',bvk.natoms,'mass',mass, ...
        'nshells',bvk.nshells,'W',{bvk.W},'shell',{bvk.shell},'nbands',bvk.nbands,'D',bvk.D,'fc',{fc});

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

function [bz]  = get_bvk_dispersion(bvk,bz)
    % set units [eV]
    units = 0.06465555; % sqrt( eV / angstrom * angstrom / amu) * hbar

    % load force constants (ignoring self-forces)
    fc = [bvk.fc{:}]; fc = fc([bvk.shell{:}]~=1);

    % get eigenvalues
    bz.hw = zeros(bvk.nbands,bz.nks); bz.U = zeros(bvk.nbands,bvk.nbands,bz.nks); 
    for i = 1:bz.nks
        % define input ...
        input = num2cell([fc,bz.k(:,i).',bvk.mass]);
        % ... and evaluate (U are column vectors)
        [bz.U(:,:,i),bz.hw(:,i)] = eig(bvk.D(input{:}),'vector'); bz.hw(:,i) = sqrt(real(bz.hw(:,i)))*units;
    end
end

function [md]  = run_bvk_md(bvk,ip,uc,dt,nsteps,Q,T)    
    % set time step [ps ~ 0.1], number of MDs steps, Nose-Hoover "mass" Q, and temperature T [K]
    % dt = 0.1; nsteps = 10000; Q = 1; T = 300;
    %
    % Uses a verlet algorithm. at present the velocities increase rapidly
    % after a couple time steps. I tried reducing the time step size to see
    % if that was related to a diff.eq. stability issue, but it didn't seem
    % to help. I am suspecting the need to implement a thermostat to keep
    % the temperature in check, i.e. Nose-Hoover.
    %
    % NOTE #1: matching to an irreducible atom may not work if this
    % instance of the irreducible atom has an orbit which is rotated, i.e.
    % this irreducible atom is related through a space symmetry which is
    % not pure translational. I will fix it when this case arises.
    %
    
    % define basic functions to enforce numeric precision
    tiny = 1E-12; mod_ = @(x) mod(x+tiny,1)-tiny; normc_ = @(A) sqrt(sum(A.^2,1));

    % loop over irreducible shells to get central atom and orbit supercell
    % indices (but organize it by the irreducible type of the center atom)
    [c_id,o_id,k2ijoy] = get_bvk_lookup_tables(ip,uc);
    
    % build force constants
    for i = unique(ip.i)
        nshells = size(k2ijoy{i},1); fc3x3 = zeros(3,3,nshells);
        for k = 1:nshells
            fc3x3(:,:,k) = ip.R(:,:,k2ijoy{i}(k,4)) * reshape(double(bvk.W{k2ijoy{i}(k,2)})*bvk.fc{k2ijoy{i}(k,2)}(:),3,3) * ip.R(:,:,k2ijoy{i}(k,4))';
        end
        % reinforce ASR
        fc3x3(:,:,1) = -sum(fc3x3(:,:,2:end),3);
        % reshape and flip sign of forces
        fc{i} = -reshape(fc3x3,3,[]);
    end

    % allocate position and velocity arrays
    tau = zeros(3,uc.natoms,nsteps); vel = zeros(3,uc.natoms,nsteps); force = zeros(3,uc.natoms);

    % set initial value conditions: random displacement between +/- 0.001 [frac]
    tau(:,:,1) = uc.tau + (.5-rand(3,uc.natoms))*0.001;
    vel(:,:,1) = zeros(3,uc.natoms); 

    % run md using verlet algorithm
    KE=zeros(1,nsteps); PE=zeros(1,nsteps); k_boltz = 0.000086173303; % [eV/K]
    for j = 2:nsteps
        % 1) get displacements [frac->cart:Ang]
        u = uc.bas * mod_( tau(:,:,j-1)-uc.tau +.5)-.5;

        % 2) compute force [eV/Ang]
        for i = 1:numel(c_id); force(:,c_id{i}) = fc{i} * reshape(u(:,o_id{i}), size(o_id{i}).*[3,1]); end

        % 3) compute potential energy [eV/Ang * Ang -> eV]
        PE(j) = - u(:).'*force(:);

        % 4) compute kinetic energy [ need a change of units here? amu * (Ang/ps)^2 -> 1.036382E-4 eV ]
        KE(j) = sum(uc.mass.*normc_(uc.bas*vel(:,:,j-1))/2); 

        % 5) compute Nose-Hoover drag: p_eta = KE - TE; degrees of freedom = 3 (1/2 PE + 1/2 KE per direction)
        nosehoover = vel(:,:,j-1)/Q * ( KE(j) - 3*uc.natoms*k_boltz*T );

        % 6) get acceleration in [frac]
        accel = uc.recbas * force ./ repmat(uc.mass,3,1);

        % update md [frac]: x' = x + v * dt; v' = v + a * dt; Nose-Hoover dv/dt becomes a - p_eta / Q * v;
        tau(:,:,j) = mod_( tau(:,:,j-1) + dt * vel(:,:,j-1) );
        vel(:,:,j) =       vel(:,:,j-1) + dt * ( accel  - nosehoover );
        
        if  j==2
            fprintf('%10s   %10s     %10s   %10s   %10s \n','step','temp','PE','KE','PE+KE');
        elseif mod(j,50)==0
            fprintf('%10i   %10f K   %10f   %10f   %10f \n',j, KE(j)/(3*uc.natoms)/k_boltz,PE(j),KE(j),PE(j)+KE(j));
            %
            figure(1); set(gcf,'color','white'); 
            plot([1:nsteps],[KE;PE;KE+PE].');legend('KE','PE','KE+PE')
            xlabel('step'); ylabel('energy [eV]'); xlim([0,nsteps]); drawnow;
        end
    end

    % define md creation function
    md_ = @(uc,tau,vel,dt,KE,PE) struct('header','md simulation','latpar',uc.latpar,...
        'bas',uc.bas,'recbas',uc.recbas,'vol',uc.vol,'symb',{{uc.symb{:}}},...
        'nspecies',uc.nspecies,'natoms',uc.natoms,'tau',tau,'vel',vel, ...
        'species',uc.species,'mass',uc.mass,'dt',dt,'nsteps',size(tau,3),'KE',KE,'PE',PE);
    md = md_(uc,tau,vel,dt,KE,PE);


    % 
    % plot positions, position and velocity histograms
    hist_v = zeros(uc.natoms,nsteps); hist_u = zeros(uc.natoms,nsteps);
    for i = 1:nsteps; hist_v(:,i) = normc_(md.vel(:,:,i)); hist_u(:,i) = normc_(mod_(md.tau(:,:,i)-uc.tau+.5)-.5); end
    for i = 1:50:nsteps 
        subplot(2,1,1); plot3(tau(1,:,i), tau(2,:,i), tau(3,:,i),'.', 'markersize', 20); view([0 0 1]); daspect([1 1 1]); drawnow; title(num2str(i));
        subplot(2,2,3); histogram(hist_v(:,i),linspace(0,max(hist_v(:)),100)); ylim([0 uc.natoms/10]); drawnow; 
        subplot(2,2,4); histogram(hist_u(:,i),linspace(0,max(hist_u(:)),100)); ylim([0 uc.natoms/10]); drawnow; 
    end

    % units = 4.1356; % [1/ps = THz -> meV]
    % 
    % psd_ = @(x) abs(fft(x)).^2;
    % 
    % psd = 0;
    % for j = 1:3
    % for i = 1:uc.natoms
    %     psd = psd + psd_(squeeze(vel(j,i,3000:end)));
    % end
    % end
    % 
    % semilogy( psd(1:end/2) ); axis tight
end

function [md]  = get_bvk_displacement(bvk,ip,uc,k,amp,mode)
    % k = [0;0;0.5];
    % amp = 1; % displacement amplitude
    % mode = 1
    % define basic functions to enforce numeric precision
    tiny = 1E-12; mod_ = @(x) mod(x+tiny,1)-tiny; normc_ = @(A) sqrt(sum(A.^2,1));
    
    % set units [eV]
    units = 0.06465555; % sqrt( eV / angstrom * angstrom / amu) * hbar

    % load force constants (ignoring self-forces)
    fc = [bvk.fc{:}]; fc = fc([bvk.shell{:}]~=1);
    % define input ...
    input = num2cell([fc,k(:).',bvk.mass]);
    % ... and evaluate (U are column vectors)
    [U,hw] = eig(bvk.D(input{:}),'vector'); hw=sqrt(real(hw))*units;
    % sort eigenvalues (and eigenvectors)
    [hw,inds] = sort(hw); U=U(:,inds);

    % print
    fprintf('Energies [meV]\n');fprintf('%5.2f \n',hw*1E3);
    fprintf('Mode selected: %i \n',mode);
    
    % get central atom and orbit supercell indices
    [c_id,o_id,k2ijoy] = get_bvk_lookup_tables(ip,uc);

    % build force constants
    fc = get_fc_matrix(bvk,ip,k2ijoy);

    % set time steps and amplitude
    nsteps=1000; t=linspace(0,1,nsteps); dt=t(2)-t(1); vel=zeros(3,uc.natoms,nsteps);
    b=zeros(1,nsteps); tau=zeros(3,uc.natoms,nsteps); KEr=zeros(1,nsteps); PEr=zeros(1,nsteps); F=zeros(3,uc.natoms);
    for i = 1:nsteps
        % get phonon coordinate [cart] (Ziman Eq. 1.6.16)
        b(i) = amp * exp(2i*pi*t(i));
        % get displacements [cart]
        u = reshape(U([1:3].'+3*(uc.u2p-1),mode),3,uc.natoms).*exp(2i*pi*k(:).'*uc.tau)./sqrt(uc.mass) .* real(b(i));
        % save displacement [frac]
        tau(:,:,i) = mod_( uc.tau + uc.recbas*real(u) );
        % get velocity [cart]
        v = reshape(U([1:3].'+3*(uc.u2p-1),mode),3,uc.natoms).*exp(2i*pi*k(:).'*uc.tau)./sqrt(uc.mass) .* imag(b(i));
        % save velocities [frac]
        vel(:,:,i) = uc.recbas * imag(v);
        
        % get references:
        % 1) potential energy [ fc [eV/Ang] * u [Ang] = [eV] ]
        for j = 1:numel(c_id); F(:,c_id{j}) = fc{j} * reshape(real(u(:,o_id{j})), size(o_id{j}).*[3,1]); end
        % 3) compute potential energy [eV]
        PEr(i) = - real(u(:)).'*F(:);
        % 4) compute kinetic energy [eV]
        KEr(i) = sum(uc.mass.*sum((uc.recbas*vel(:,:,i)).^2,1))/2;
    end
    
    % get potential energy (Ziman Eq. 1.6.17)
    PE = abs(hw(mode)).^2 * real(b).^2 / 2;
    
    % get kinetic energy (Ziman Eq. 1.6.17; this seems wrong...)
    KE = abs(hw(mode)).^2 * imag(b).^2 / 2;

    md_ = @(uc,tau,vel,dt,KE,PE) struct('units','frac','latpar',uc.latpar,...
        'bas',uc.bas,'recbas',uc.recbas,'vol',uc.vol,'symb',{{uc.symb{:}}},...
        'nspecies',uc.nspecies,'natoms',uc.natoms,'tau',tau,'vel',vel, ...
        'species',uc.species,'mass',uc.mass,'dt',dt,'nsteps',size(tau,3),'KE',KE,'PE',PE);
    md = md_(uc,tau,vel,dt,KE,PE);
    
    figure(1); plotyy([1:nsteps],KEr,[1:nsteps],KE); legend('REF','KE');
    figure(2); plotyy([1:nsteps],PEr,[1:nsteps],PE); legend('REF','PE');
    
    function fc = get_fc_matrix(bvk,ip,k2ijoy)
        % build force constants
        fc(unique(ip.i))=cell(1);
        for ii = unique(ip.i)
            nshells = size(k2ijoy{ii},1); fc3x3 = zeros(3,3,nshells);
            for ki = 1:nshells
                fc3x3(:,:,ki) = ip.R(:,:,k2ijoy{ii}(ki,4)) * reshape(double(bvk.W{k2ijoy{ii}(ki,2)})*bvk.fc{k2ijoy{ii}(ki,2)}(:),3,3) * ip.R(:,:,k2ijoy{ii}(ki,4))';
            end
            % enforce ASR
            fc3x3(:,:,1) = -sum(fc3x3(:,:,2:end),3);
            % reshape and flip sign of forces
            fc{ii} = -reshape(fc3x3,3,[]);
        end
    end
end

function         plot_bvk_dispersion(bvk,bzp)
    % get phonon band structure along path
    bzp = get_bvk_dispersion(bvk,bzp);

    % and plot the results
    fig_ = @(h)       set(h,'color','white');
    axs_ = @(h,qt,ql) set(h,'Box','on','XTick',qt,'Xticklabel',ql);

    figure(1); fig_(gcf);
    plot(bzp.x,sort(real(bzp.hw)*1E3),'-k',bzp.x,-sort(abs(imag(bzp.hw))),':r');
    axs_(gca,bzp.qt,bzp.ql); axis tight; ylabel('Energy [meV]'); xlabel('Wavevector k');
end


% electrons

function [tb,ip] = get_tb(cutoff,pc,spdf,nskips,Ef,fname,flags)
    % for paper:
    % cutoff = 3; % Angstroms 
    % fname = 'EIGENVAL'
    
    % continue earlier calc?
    sfile = sprintf('am_tb.mat');
    if and(strfind(flags,'continue'),exist(sfile,'file')); load(sfile); return; end 

    % get shells
    ip = get_irreducible_shells(pc,cutoff);

    % get tb model
    tb = get_tb_model(ip,spdf);

    % get force constants
    tb = get_tb_matrix_elements(tb,nskips,Ef,fname);

    % save everything generated
    save(sfile,'tb','ip');
end

function [tb] = get_tb_model(ip,spdf)
    % set oribtals per irreducible atom: spdf = {'d','p'};
    %
    % NOTE #1: kpoints are in fractional reciprocal units but force
    % constants are in cartesian units.

    % define basic functions to enforce numeric precision
    tiny = 1E-12; rnd_ = @(x) round(x,-log10(tiny)); modc_ = @(x) rnd_(mod(x+tiny+.5,1)-tiny-.5); 
    normc_ = @(A) sqrt(sum(A.^2,1)); findrow_ = @(A) (sum(cumsum(A~=0,2)==0,2)+1).';

    % initialize atoms: J = azimuthals, D = symmetry reps, F = flip super-operator
    [J,D,F] = initialize_atoms(spdf,ip.R);

    % primitive cell atoms define hamiltonian blocks dimensions and start/end sections
    p2i(ip.m)=ip.i; p2i(ip.n)=ip.j; pc_atoms=unique([ip.m,ip.n]); d(pc_atoms)=0; 
    for p=pc_atoms; d(p)=sum(J{p2i(p)}*2+1); end; E=cumsum(d); S=E-(d)+1; nbands=E(end);

    % initialize dynamical matrix and kvector
    H = sym(zeros(nbands)); kvec=sym('k%d',[3,1],'real'); [~,~,shell]=unique(rnd_(normc_(ip.v))); shell=shell(:).';
    alloc_ = @(cel) struct('W',cel,'shell',cel); sav = alloc_({cell(1,ip.nshells)});
    for p = 1:ip.nshells
        % set abbreviations
        dij=d(ip.i(p))*d(ip.j(p));

        % use stabilzer group to determine crystallographic symmetry relations; A*B*C' equals kron(C,A)*B(:)
        W = sym(zeros(dij)); for i = find(ip.s_ck(:,p)).'; W = W + kron( D{ip.j(p)}(:,:,i),D{ip.i(p)}(:,:,i)) - eye(dij); end

        % partity transpose 
        if (ip.i(p)==ip.j(p)); W = W + F{ip.i(p)} - eye(dij); end

        % get nullspace and normalize to first nonzero element
        W = null(W); for i = 1:size(W,2); W(:,i) = W(:,i)/W(find(W(:,i),1),i); end; W = rref(W.').';

        % define parameters
        c = sym(sprintf('c%02i_%%d%%d',p),[d(ip.i(p)),d(ip.j(p))],'real'); c = c(findrow_(double(W).'));

        % print symmetry adapted matrix elements
        fprintf('v = [%10f,%10f,%10f] \n', ip.v(:,p)); vsk = reshape(W*c(:),d([ip.i(p),ip.j(p)])); disp( vsk )

        % get block of dynamical matrix (see NOTE #1)
        B = zeros( d( ip.i(p) ) , d( ip.j(p) ) ); 
        for i = find(ip.g_ck(:,p)).'
            B = B + D{ip.i(p)}(:,:,i) * vsk * D{ip.j(p)}(:,:,i)' ....
                .* exp( sym(2i*pi) * dot(rnd_(ip.recbas*ip.R(:,:,i)*ip.v(:,p)),kvec) ); 
        end

        % augment hamiltonian
        mp = S(ip.m(p)):E(ip.m(p)); np = S(ip.n(p)):E(ip.n(p));
        H( mp,np ) = H( mp,np ) + B;  if (ip.m(p)~=ip.n(p)); H( np,mp ) = H( np,mp ) + B'; end

        % save important stuff
        [~,n] = sort(c); sav.W{p} = W(:,n); sav.shell{p} = ones(size(c))*shell(p); 
    end

    % simplify the dynamical matrix
    H = simplify(H,'steps',500);

    % build tb structure
    tb_ = @(ip,sav,H,spdf,nbands) struct('units','cart,frac-recp','bas',ip.bas,'recbas',ip.recbas,'natoms',numel(unique(ip.m)),'spdf',{spdf}, ...
        'nshells',size(sav.W,2),'W',{sav.W},'shell',{sav.shell},'nbands',nbands,'H',matlabFunction(H));
    tb = tb_(ip,sav,H,spdf,nbands);
end

function [tb] = get_tb_matrix_elements(tb,nskips,Ef,fname)
    % Ef     : fermi energy (read from OUTCAR)
    % nskips : number of dft bands to skip (e.g. 5)
    % fname  : eigenval file (e.g. 'EIGENVAL')

    % load dispersion [frac-recp] and shift Fermi energy to zero
    [dft,bz]=load_vasp_eigenval(fname); dft.E = dft.E - Ef; 

    % define randomizer
    rand_ = @(x) (0.5-rand(size(x))).*abs(x./max(x));

    % set array parameter size calculator for selected shells
    nxs_ = @(selector_shell) sum(any([tb.shell{:}]==selector_shell(:),1));

    % fit gamma point band energies to zero-th neighbor shell parameters
    selector_shell = [1]; selector_kpoint = 1; x=zeros(1,nxs_(selector_shell));
    x_best = optimize_and_plot(tb,bz,dft,x,selector_kpoint,selector_shell,nskips);

    % fit neighbor parameter at high symmetry points using poor man's simulated anneal
    shells=unique([tb.shell{:}]);
    for j = shells(2:end)

        % select shells and kpoints
        selector_shell=[selector_shell,j]; selector_kpoint = round(linspace(1,bz.nks,10)); kT = 20;

        for i = 1:20
            if i == 1 % initialize
                x = zeros(1,nxs_(selector_shell)); x(1:numel(x_best)) = x_best; r_best = Inf;
            else % modify x values with a damped temperature function
                x = x_best + kT * rand_(x_best) * exp(-i/25);
            end

            % optimize
            [x,r] = optimize_and_plot(tb,bz,dft,x,selector_kpoint,selector_shell,nskips);

            % save r_best parameter
            if r < r_best; r_best = r; x_best = x; end
        end
    end

    % five final last pass with all parameters and all kpoints
    selector_shell=shells; selector_kpoint=[1:bz.nks]; kT = 5;
    [x,r] = optimize_and_plot(tb,bz,dft,x,selector_kpoint,selector_shell,nskips);

    % save r_best parameter
    if r < r_best; r_best = r; x_best = x; end

    % save refined matrix elements and conform to bvk
    for i = [1:tb.nshells]; d(i)=size(tb.W{i},2); end; Evsk=cumsum(d); Svsk=Evsk-d+1;
    for i = [1:tb.nshells]; tb.vsk{i} = x(Svsk(i):Evsk(i)); end;
    
    %%

    function H = compute_tb_hamiltonian(tb,x,k)
        % define padding
        pad_ = @(A,n) [A,zeros(1,n(2)-size(A,2));zeros(n(1)-size(A,1),n(2))];
        % pad tight binding matrix elements
        nargs = nargin(tb.H); padded_x = pad_(x,[1,nargs-3]);
        % get hamiltonians
        nks = size(k,2); H = zeros(tb.nbands,tb.nbands,nks);
        for m = 1:nks
            % build input
            input = num2cell([padded_x,k(:,m).']);
            % evaluate H
            H(:,:,m) = tb.H(input{:});
        end
    end

    function [R] = compute_tb_residual(tb,bz,dft,x,selector_kpoint,nskips)
        % generate input
        H = compute_tb_hamiltonian(tb,x,bz.k(:,selector_kpoint));
        % allocate space for residual vector
        R = zeros(1,numel(selector_kpoint)*tb.nbands);
        % determine residual vector
        inds = [1:tb.nbands]+nskips;
        for n = 1:size(H,3)
            % build residual vector
            R([1:tb.nbands]+(n-1)*tb.nbands) = dft.E(inds,selector_kpoint(n)) - sort(real(eig( H(:,:,n) )));
        end
    end

    function [x,r] = optimize_and_plot(tb,bz,dft,x,selector_kpoint,selector_shell,nskips)
        % trim x according to selector_shell NOTE shells must be in
        % monotonically increasing order and not skip any shells!
        trim = any([tb.shell{:}]==selector_shell(:),1); x = x(trim(1:numel(x)));

        % define cost functional
        cost_ = @(x) compute_tb_residual(tb,bz,dft,x,selector_kpoint,nskips); 

        % define optimization parameters
        opts = optimoptions('lsqnonlin','Display','none','MaxIter',7);

        % perform optimization
        [x,r] = lsqnonlin(cost_,x,[],[],opts);

        % plot band structure (quick and dirty)
        H = compute_tb_hamiltonian(tb,x,bz.k); E = zeros([size(H,1),size(H,2)]);
        for m = 1:bz.nks; E(:,m)=sort(real(eig( H(:,:,m) ))); end
        figure(1); plot([1:bz.nks],real(E),'-k', [1:bz.nks],dft.E([1:(end-nskips)]+nskips,:),':r');
        set(gca,'XTick',[]); axis tight; grid on; ylabel('Energy E'); xlabel('Wavevector k'); drawnow;
    end

end

function [bz] = get_tb_dispersion(tb,bz)
    % get eigenvalues
    bz.E = zeros(tb.nbands,bz.nks); bz.V = zeros(tb.nbands,tb.nbands,bz.nks);
    for i = 1:bz.nks
        % define input ...
        input = num2cell([tb.vsk{:},bz.k(:,i).']);
        % ... and evaluate (V are column vectors)
        [bz.V(:,:,i),bz.E(:,i)] = eig(tb.H(input{:}),'vector'); 
    end
end

function [ibz] = get_nesting(tb,ibz,bzp,Ef,degauss)
    % Ef, fermi energy
    % degauss, degauss = 0.04 61x61x61 kpoint mesh

    % define input ...
    input = num2cell([tb.vsk{:},bzp.k(:,ik).']);
    % ... and evaluate
    bzp.E(:,ik) = sqrt(eig(tb.H(input{:})));
    
    % define gaussian function 
    lorentz_ = @(x) 1./(pi*(x.^2+1)); % gauss_ = @(x) exp(-abs(x).^2)./sqrt(pi); 

    % evaluate eigenvalues on ibz and save to fbz mesh
    nbands = 8; E = zeros(nbands,numel(f2i));
    for i=[1:ibz.nks]; E(:,f2i==i)=repmat(eig(getH(v,(ibz.recbas*ibz.latpar)*ibz.k(:,i)),'vector'),1,ibz.w(i)); end

    % define nonzero, reshape array into tensor, flatten
    expand_ = @(x) reshape(x,ibz.n); flatten_ = @(x) x(:);

    % define fermi list
    Ef_list = [-1:0.05:0.6];

    % define figure plots
    colormap(get_colormap('magma',256)); [xx,yy]=meshgrid(bzp.x,Ef_list);
    fig_ = @(h) set(h,'color','white'); axs_ = @(h,qt,ql) set(h,'Box','on','XTick',qt,'Xticklabel',ql);

    % loop over fermi levels
    nEfs = numel(Ef_list); nk_vs_ef = zeros(bzp.nks,nEfs); 
    for j = 1:nEfs
        Ef = Ef_list(j);

        % compute spectral function A on the full mesh
        A = expand_( sum(lorentz_((E-Ef)/degauss)/degauss,1) );

        % compute nesting on ibz and save on fbz mesh
        n=zeros(ibz.n); for i=[1:ibz.nks]; [q{1:3}]=ind2sub(ibz.n,i2f(i)); n(f2i==i)=(A(:).'*flatten_(circshift(A,[q{:}]-1)))/prod(ibz.n); end

        % interpolate n on [frac] path in strides (memory overload otherwise)
        nk=zeros(bzp.nks,1); stride=100; for i=[1:stride]'+[0:stride:(bzp.nks-stride)]; nk(i) = real(fftinterp(n,(bzp.bas)*bzp.k(:,i))); end

        nk_vs_ef(:,j) = nk;

        % plot figure
        figure(1); fig_(gcf); surf(xx.',yy.',nk_vs_ef); 
        view(2); axis tight; axs_(gca,bzp.qt,bzp.ql); drawnow;
    end
    
    
    
    
    
    
    % define gaussian function (N=61,degauss=0.04 is good)
    lorentz_ = @(x) 1./(pi*(x.^2+1)); % gauss_ = @(x) exp(-abs(x).^2)./sqrt(pi); 

    % evaluate eigenvalues on ibz and save to on fbz mesh
    nbands = 8; E = zeros(nbands,numel(f2i));
    for i = 1:ibz.nks; E(:,f2i==i) = repmat(eig(getH(v,(ibz.recbas*ibz.latpar)*ibz.k(:,i)),'vector'),1,ibz.w(i)); end

    % define nonzero, reshape array into tensor, flatten
    expand_ = @(x) reshape(x,ibz.n); flatten_ = @(x) x(:);

    % compute spectral function A on the full mesh
    A = expand_( sum(lorentz_((E-Ef)/degauss)/degauss,1) );

    % compute nesting on ibz and save on fbz mesh
    n = zeros(ibz.n); m_ = @(i,j) mod(i-1,ibz.n(j))+1; fbz_nks=numel(f2i);
    for i = 1:ibz.nks; [q{1:3}]=ind2sub(ibz.n,i2f(i)); n(f2i==i) = (A(:).'*flatten_(circshift(A,[q{:}]-1)))/fbz_nks; end

    % define kpoint path (cart)
    G=[0;0;0]; X1=[0;1;1]/2; X2=[2;1;1]/2; L=[1;1;1]/2; K=[6;3;3]/8; iM = ones(3)-2*eye(3); M=inv(iM);
    Np = 300; ql={'G','X','K','G','L'}; qs=iM*[G,X2,K,G]; qe=iM*[X1,K,G,L]; nqs = size(qs,2);

    % build path (transform by iM to correct distances ~ 1/nm)
    [k,x,qt] = get_path(qs,qe,nqs,Np); nks = numel(x); k = M*k;

    % interpolate n on the path in strides (memory overload otherwise)
    nk = zeros(nks,1); stride=100; for i = [1:stride]'+[0:stride:(nks-stride)]; nk(i) = fftinterp(n, k(:,i)); end

    % define figure properties
    fig_ = @(h)       set(h,'color','white');
    axs_ = @(h,qt,ql) set(h,'Box','on','XTick',qt,'Xticklabel',ql);

    % plot band structure 
    figure(1); fig_(gcf);
    semilogy(x,real(nk),'-'); axs_(gca,qt,ql); axis tight; 
    ylabel('Nesting [a.u.]'); set(gca,'YTickLabel',[]); xlabel('Wavevector k'); 

end

function get_nesting_vs_ef()
    % define basic functions to enforce numeric precision
    tiny = 1E-12; mod_ = @(x) mod(x+tiny,1)-tiny; rnd_ = @(x) round(x,-log10(tiny));

    % define tight binding matrix elements
    v = [-0.4212,1.1975,-4.1841,-1.0193,-1.0322,-0.0565,0.1132,-0.5218,-0.1680,0.0635,-0.0546,-0.1051,0.4189,0.3061];

    % define gaussian function (N=61,degauss=0.04 is good)
    degauss = 0.04; lorentz_ = @(x) 1./(pi*(x.^2+1)); % gauss_ = @(x) exp(-abs(x).^2)./sqrt(pi); 

    % evaluate eigenvalues on ibz and save to fbz mesh
    nbands = 8; E = zeros(nbands,numel(f2i));
    for i=[1:ibz.nks]; E(:,f2i==i)=repmat(eig(getH(v,(ibz.recbas*ibz.latpar)*ibz.k(:,i)),'vector'),1,ibz.w(i)); end

    % define nonzero, reshape array into tensor, flatten
    expand_ = @(x) reshape(x,ibz.n); flatten_ = @(x) x(:);

    % define fermi list
    Ef_list = [-1:0.05:0.6];

    % define figure plots
    colormap(get_colormap('magma',256)); [xx,yy]=meshgrid(bzp.x,Ef_list);
    fig_ = @(h) set(h,'color','white'); axs_ = @(h,qt,ql) set(h,'Box','on','XTick',qt,'Xticklabel',ql);

    % loop over fermi levels
    nEfs = numel(Ef_list); nk_vs_ef = zeros(bzp.nks,nEfs); 
    for j = 1:nEfs
        Ef = Ef_list(j);

        % compute spectral function A on the full mesh
        A = expand_( sum(lorentz_((E-Ef)/degauss)/degauss,1) );

        % compute nesting on ibz and save on fbz mesh
        n=zeros(ibz.n); for i=[1:ibz.nks]; [q{1:3}]=ind2sub(ibz.n,i2f(i)); n(f2i==i)=(A(:).'*flatten_(circshift(A,[q{:}]-1)))/prod(ibz.n); end

        % interpolate n on [frac] path in strides (memory overload otherwise)
        nk=zeros(bzp.nks,1); stride=100; for i=[1:stride]'+[0:stride:(bzp.nks-stride)]; nk(i) = real(fftinterp(n,(bzp.bas)*bzp.k(:,i))); end

        nk_vs_ef(:,j) = nk;

        % plot figure
        figure(1); fig_(gcf); surf(xx.',yy.',nk_vs_ef); 
        view(2); axis tight; axs_(gca,bzp.qt,bzp.ql); drawnow;
    end

    % define figure properties
    colormap((get_colormap('magma',256)));
    figure(1); fig_(gcf); surf(xx.',yy.',log(nk_vs_ef./max(nk_vs_ef,[],1))); shading interp;
    view(2); axis tight; axs_(gca,bzp.qt,bzp.ql); ylabel('Energy [eV]'); xlabel('wavevector k'); 
    % caxis(log([0.14 0.7])); ylim([-1 0.4]); set(gca,'Ytick',[-1:0.2:0.4]);

    set(gca,'LooseInset',get(gca,'TightInset')); set(gcf,'PaperSize',[10 10]);
    set(gcf,'PaperPosition',[0,0,1.6180,1]*3); print(gcf,'-djpeg','-r600','nesting_vs_Ef.jpeg');
end

function [J,D,F] = initialize_atoms(spdf_list,R)

    % get symmetries
    nRs=size(R,3);

    % transform symmetries to the tight binding representation (wiger functions)
    W=cell(1,3); for j=[1:3]; W{j} = get_wigner(j,R); end

    % set orbitals J{:}, symmetries D{:}, and parity-transpose T{:} for each irreducible atom
    natoms=numel(spdf_list); F=cell(1,natoms);  D=cell(1,natoms);
    for i = 1:natoms
        % set orbitals
        J{i} = spdf2J(spdf_list{i});

        % set start and end points for J
        E=cumsum(J{i}*2+1); S=E-(J{i}*2+1)+1;

        % construct D matrix and lay the ground work construction of parity super-operator
        d = max(E); P = zeros(1,d); D{i} = sym(zeros(d,d,nRs));
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

    function [J] = spdf2J(spdf)
        o = lower(strtrim(spdf)); norbitals=length(o); J=zeros(1,norbitals);
        for l = 1:norbitals
            J(l) = find(strcmp(o(l),{'s','p','d','f'}))-1;
        end
    end
end


% aux brillouin zones

function tet = get_tetrahedra(recbas,n)
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
end

function box = grid2box(n)
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

function tetrahedron = box2tetrahedron(recbas)
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

% aux phonons

function [c_id,o_id,k2ijor] = get_bvk_lookup_tables(ip,uc)
    % define basic functions to enforce numeric precision
    tiny = 1E-12; mod_ = @(x) mod(x+tiny,1)-tiny; 
    
    % define comparison function
    match2_ = @(q,k) find(all(abs(q-k)<tiny,1)); % index of q(3,1) in k(3,:) or [] otherwise

    % define high dimensional matrix multiplication: A [(m,n),(a,b)] * B [(n,p),(b,c)] = C [(m,p),(a,c)]
    aug_ = @(x,i,y) [x(1:(i-1)),y,x(i:end)]; matmul_ = @(A,B) squeeze(sum(bsxfun(@times, reshape(A,aug_(size(A),3,1)), reshape(B,aug_(size(B),1,1))), 2));

    % count number of orbits involving primitive cell atom i
    norbits_ = @(i) sum(ip.norbits( or(ip.i==i,ip.j==i)) );

    % loop over irreducible shells to get central atom and orbit supercell indices 
    % organize it by the irreducible type of the center atom
    for i = unique(ip.i)
        c_id{i} = find(uc.u2i==i); nstars = numel(c_id{i}); k = 0;
        o_id{i} = zeros(norbits_(i),nstars);
        k2ijor{i} = zeros(norbits_(i),4);
        for j = 1:ip.nshells; if or(ip.i(j)==i,ip.j(j)==i)        
            y = find(ip.g_ck(:,j)).';
            for o = 1:ip.norbits(j); k = k + 1;
                % generate orbit [cart -> frac]
                orbit = mod_( uc.recbas * matmul_(ip.R(:,:,ip.g_ck(:,j)),ip.v(:,j)) );
                % indices: irreducible atom i, irreducible shell j, orbit o, symmetry y
                k2ijor{i}(k,:)=[i,j,o,y(o)];
                % identify indices of orbit points for each atom in the supercell atoms of irreducible type m (SEE NOTE #1)
                for s = 1:nstars; o_id{i}(k,s) = match2_( orbit(:,o), mod_(uc.tau-uc.tau(:,c_id{i}(s)) ) ); end
            end
        end; end
    end
end

% aux electrons

function [Wtes] = get_wigner(j,R)

    % define tiny, kronecker delta, and heavside
    tiny = 1E-8; d_ = @(x,y) abs(x-y)<tiny; t_ = @(x,y) logical(x>y);

    % matrix indices
    [m,mp]=meshgrid([j:-1:-j]);

    % define angular momentum operators (Jp raising, Jm lowering, ...)
    Jm = d_(m,mp+1).*sqrt((j+m).*(j-m+1)); Jp = d_(m,mp-1).*sqrt((j-m).*(j+m+1));
    Jx = (Jp+Jm)/2; Jy = (Jp-Jm)/2i; Jz = d_(m,mp).*m; J = cat(3,Jx,Jy,Jz);

    % define basis change: spherical to tesseral harmonics (real basis)
    T = d_(0,m) .* d_(mp,m) + ...
        t_(m,0) .* - sqrt(-1/2) .* ( (-1).^m.*d_(m,-mp) - d_(m, mp) ) + ...
        t_(0,m) .* - sqrt( 1/2) .* ( (-1).^m.*d_(m, mp) + d_(m,-mp) );

    % batch convert to wigner
    nRs = size(R,3); Wtes = sym(zeros(2*j+1,2*j+1,nRs));
    for i = [1:nRs]; Wtes(:,:,i) = get_wigner_engine(J,T,j,sym(R(:,:,i))); end
end

function [Wtes] = get_wigner_engine(J,T,j,R)
    % define wigner function for spherical (complex) and tesseral (real) harmonics
    % Note: for l = 1, Wtes_(R) = R
    
    % define conversion of proper rotations to axis & angle representation
    ang_ = @(R) acos((trace(R)-1)/2); ax_ = @(R) circshift( get_ax(double(R)), 1);

    % define spin-vector dot products [Eq. 1.37 Blundell]
    dotV_ = @(S,V) S(:,:,1)*V(1) + S(:,:,2)*V(2) + S(:,:,3)*V(3);

    if isa(R, 'sym') % faster for symbolic and squares matrices with dimensions > 9
        Wsph = expm( -sqrt(-1) * dotV_( sym(J), ax_(det(R)*R) ) * sym( ang_(det(R)*R) ) ) * det(R)^j;
        Wtes = simplify( sym(T') * Wsph * sym(T) );
    else % faster for numerical square matrices with dimensions < 9
        [V,D] = eig( -sqrt(-1) * dotV_( J, ax_(det(R)*R) ) * ang_(det(R)*R) ); 
        Wsph = V*diag(exp(diag(D)))/V * det(R)^j;
        Wtes = T' * Wsph * T ;
    end

end

function [ax] = get_ax(R)
    % define basic parameters and functions
    tiny = 1E-8; normalize_ = @(v) v/norm(v); 

    % check for identity
    if abs(trace(R)-3)< tiny; ax = [0;0;1]; return; end

    % get rotation axis
    ax = null(R-eye(3));

    % get random point on plane perpendicular to the rotation axis
    v1 = rand(3,1); v1 = normalize_(v1 - dot(v1,ax)*ax);

    % rotate point on the perpendicular plane
    v2 = R*v1;

    % get cross product (add tiny cross component to deal with 180 deg rotation)
    c = normalize_(cross(v1,v2+cross(v1,ax)*tiny));

    % adjust sign
    ax = sign(dot(c,ax))*ax;
end

% aux functions


function [K0] = uc2ws(K,M)
    % uc2ws uses M real (reciprocal) lattice vectors to reduces K(1:3,:) vectors 
    % in cartesian (reciprocal) coordinates to the definiging Wigner-Seitz cell.
    % Note: K0 = cell2mat(arrayfun(@(j) uc2ws_engine(K(:,j),G,G2),[1:size(K,2)],'unif',0)); 
    %       is slower than looping.
    
    % generate mesh
    [X{1:3}] = meshgrid([-1:1]); G=M*[X{1}(:),X{2}(:),X{3}(:)].'; G2=sum(G.^2,1);

    % call engine
    m=size(K,2); K0 = zeros(3,m); for j = 1:m; K0(:,j) = uc2ws_engine(K(:,j),G,G2); end
    
    function [K] = uc2ws_engine(K,G,G2)
        % define tiny
        tiny = 1E-12;
        
        go = true;
        while go; go=false;
            for i = 1:27
                P = 2*(K).'*G(:,i)-G2(i);
                if P > + tiny
                    K = K - G(:,i); go=true;
                end
            end
        end
    end
end

function [fq] = fftinterp(f,q)
    % fourier interpolate f(k) at points q; v must be periodic over [0,1)
    %
    % generate f using like this:
    %
    % mpgrid_ = @(N) [0:(N-1)]/N;
    % [xk,yk,zk]=meshgrid(mpgrid_(n(1)),mpgrid_(n(2)),mpgrid_(n(3)));
    % k = [xk(:),yk(:),zk(:)];
    % f = cos(2*pi*xk)+cos(2*pi*yk)+cos(2*pi*zk);
    %

    % define flatten
    flatten_ = @(x) x(:);
    
    % mesh dimensions
    n = size(f);

    % generate Fourier mesh
    fftmesh_ = @(N) [0:(N-1)]-floor(N/2);
    [x,y,z]=meshgrid(fftmesh_(n(1)),fftmesh_(n(2)),fftmesh_(n(3)));
    r = [x(:),y(:),z(:)].';

    % construct inverse Fourier transform kernel
    Ki = ones(size(q,2),size(r,2)); i2pi = sqrt(-1)*2*pi;
    for i = 1:3; Ki = Ki.*exp(i2pi*q(i,:).'*r(i,:))./sqrt(n(i)); end

    % perform interpolation 
    fq = Ki * flatten_(fftshift(fftn( f ))) / sqrt(prod(n));
end

function [A] = frref(A)
%FRREF   Fast reduced row echelon form.
%   R = FRREF(A) produces the reduced row echelon form of A.
% 
%   Description: 
%   For full matrices, the algorithm is based on the vectorization of MATLAB's
%   RREF function. A typical speed-up range is about 2-4 times of 
%   the MATLAB's RREF function. However, the actual speed-up depends on the 
%   size of A. The speed-up is quite considerable if the number of columns in
%   A is considerably larger than the number of its rows or when A is not dense.
%
%   For sparse matrices, the algorithm ignores the tol value and uses sparse
%   QR to compute the frref form, improving the speed by a few orders of 
%   magnitude.
%
%   Authors: Armin Ataei-Esfahani (2008)
%            Ashish Myles (2012)
%
%   Revisions:
%   25-Sep-2008   Created Function
%   21-Nov-2012   Added faster algorithm for sparse matrices

    % set size and default tolerences
    [m,n] = size(A); tol=max(m,n)*eps(class(A))*norm(A,'inf');

    % Loop over the entire matrix.
    i = 1; j = 1;
    % t1 = clock;
    while (i <= m) && (j <= n)
       % Find value and index of largest element in the remainder of column j.
       [p,k] = max(abs(A(i:m,j))); k = k+i-1;
       if (p <= tol)
          % The column is negligible, zero it out.
          A(i:m,j) = 0; %(faster for sparse) %zeros(m-i+1,1);
          j = j + 1;
       else
          % Swap i-th and k-th rows.
          A([i k],j:n) = A([k i],j:n);
          % Divide the pivot row by the pivot element.
          Ai = A(i,j:n)/A(i,j);    
          % Subtract multiples of the pivot row from all the other rows.
          A(:,j:n) = A(:,j:n) - A(:,j)*Ai;
          A(i,j:n) = Ai;
          i = i + 1;
          j = j + 1;
       end
    end
end

function [F] = fsrref(A)
    % fast sparse rref based on QR decomposition
    [m,n] = size(A);

    % Non-pivoted Q-less QR decomposition computed by Matlab actually
    % produces the right structure (similar to rref) to identify independent
    % columns.
    R = qr(A);
    % i_dep = pivot columns = dependent variables
    %       = left-most non-zero column (if any) in each row
    % indep_rows (binary vector) = non-zero rows of R
    [indep_rows, i_dep] = max(R ~= 0, [], 2);
    indep_rows = full(indep_rows); % probably more efficient
    i_dep = i_dep(indep_rows);
    i_indep = setdiff(1:n, i_dep);

    % solve R(indep_rows, i_dep) x = R(indep_rows, i_indep)
    %   to eliminate all the i_dep columns
    %   (i.e. we want F(indep_rows, i_dep) = Identity)
    F = sparse([],[],[], m, n);
    F(indep_rows, i_indep) = R(indep_rows, i_dep) \ R(indep_rows, i_indep);
    F(indep_rows, i_dep) = speye(length(i_dep));

    % result
    A = F;
end

function [y] = nlinspace(d1, d2, n)
    %LINSPACENDIM Linearly spaced multidimensional matrix.
    
    if nargin == 2
        n = 100;
    end
    n  = double(n);
    d1 = squeeze(d1); d2 = squeeze(d2);

    if ndims(d1)~= ndims(d2) || any(size(d1)~= size(d2))
        error('d1 and d2 must have the same number of dimension and the same size'),
    end

    NDim = ndims(d1);
    %%%%%%%% To know if the two first dimensions are singleton dimensions
    if NDim==2 && any(size(d1)==1)
        NDim = NDim-1;
        if all(size(d1)==1)
            NDim = 0;
        end
    end

    pp      = (0:n-2)./(floor(n)-1);

    Sum1    = TensorProduct(d1, ones(1,n-1));
    Sum2    = TensorProduct((d2-d1), pp);
    y = cat(NDim+1, Sum1  + Sum2, shiftdim(d2, size(d1, 1)==1 ));

    %%%%% An old function that I wrote to replace the built in Matlab function:
    %%%%% KRON
    function Z = TensorProduct(X,Y)
        %   Z = TensorProduct(X,Y) returns the REAL Kronecker tensor product of X and Y. 
        %   The result is a multidimensional array formed by taking all possible products
        %   between the elements of X and those of Y. 
        %

        sX=size(X);sY=size(Y);

        ndim1=ndims(X);ndim2=ndims(Y);

        indperm=[ndim2+1:ndim1+ndim2,1:ndim2];

        % to remove all singleton dimensions 
        Z=squeeze(repmat(X,[ones(1,ndims(X)),sY]).*permute(repmat(Y,[ones(1,ndims(Y)),sX]),indperm));
    end
end

% aux aesthetic

function [cmap] =  get_colormap(palette,n)

% color palettes
switch palette
case 'spectral' % sns.color_palette("Spectral", 10))
cmap = [ ...
 0.81414841553744144, 0.21968473703143937, 0.30480585554066825;
 0.93302576331531295, 0.39131103777417953, 0.27197233193060932;
 0.98177624099394856, 0.60738179087638855, 0.34579008992980509;
 0.99469434864380779, 0.80922723167082844, 0.48696657138712268;
 0.99823144954793597, 0.94517493598601399, 0.65705499929540301;
 0.95578623869839840, 0.98231449547935934, 0.68004615517223588;
 0.82029989537070780, 0.92756632496328917, 0.61268745099796973;
 0.59100347582031698, 0.83552480795804196, 0.64429067864137535;
 0.36001538412243711, 0.71618609919267540, 0.66551328406614418;
 0.21299500558890549, 0.51141871132102668, 0.73079586379668293];
case 'RdGy' % sns.color_palette("RdGy", 10)
cmap = [ ... 
 0.66920416904430757, 0.08489042841920666, 0.16401384522517523;
 0.81153403660830326, 0.32110727271612954, 0.27581700390460445;
 0.92226067360709696, 0.56747406545807333, 0.44867361117811766;
 0.97970011655022116, 0.78408305785235233, 0.68489044554093303;
 0.99646289909587182, 0.93633218372569371, 0.90096117468441239;
 0.94517493598601399, 0.94517493598601399, 0.94517493598601399;
 0.82583622722064742, 0.82583622722064742, 0.82583622722064742;
 0.67058825492858887, 0.67058825492858887, 0.67058825492858887;
 0.48481355812035354, 0.48481355812035354, 0.48481355812035354;
 0.28235295195789900, 0.28235295195789900, 0.28235295195789900];
case 'OrRd' % sns.color_palette("OrRd", 10)
cmap = [ ... 
 0.99717031927669753, 0.92618224200080423, 0.82362169448067157;
 0.99434063855339494, 0.87504806588677797, 0.71132643012439500;
 0.99215686321258545, 0.81522492450826312, 0.60281432586557726;
 0.99215686321258545, 0.74140716650906735, 0.52604385754641370;
 0.98965013518052947, 0.61802385975332819, 0.40985776314548417;
 0.96984237011741192, 0.49634757789911010, 0.32496733069419859;
 0.92950404111076801, 0.37896194370353925, 0.26911189196740881;
 0.85863899343154015, 0.22246828552554634, 0.14805075201918097;
 0.76452135198256554, 0.08341407314235088, 0.05387158794145958;
 0.64518263491929750, 0.00000000000000000, 0.00000000000000000];
case 'GnBu' % sns.color_palette("GnBu", 10)
cmap = [ ... 
 0.90354479621438422, 0.96276816830915568, 0.88175317890503824;
 0.84367551873711977, 0.93903883485233086, 0.82059209066278793;
 0.77674741815118231, 0.91252595677095305, 0.76221454704509062;
 0.67044984663234042, 0.87118801229140341, 0.71497118192560527;
 0.54602076925483400, 0.82405229806900027, 0.74740485934650192;
 0.41868512595401092, 0.76462900287964763, 0.78985007019603959;
 0.29457901909070855, 0.68936564711963433, 0.82066898416070377;
 0.19123414667213665, 0.57420994464088881, 0.75866206744137932;
 0.09219531267881393, 0.47040370585871677, 0.70579009827445538;
 0.03137255087494850, 0.36416763593168822, 0.62755865349489104];
case 'YlGn' % sns.color_palette("YlGn", 10)
cmap = [ ... 
 0.97736255421357998, 0.99151095783009247, 0.77353326993830063;
 0.91649366126340981, 0.96738177818410542, 0.68725876527674057;
 0.82256056420943313, 0.92890427182702462, 0.62565169474657845;
 0.69264131013084862, 0.87280277574763576, 0.56364477802725399;
 0.54557478696692230, 0.80901192987666415, 0.50422146378778943;
 0.39277201715637655, 0.73826991179410151, 0.43489427531466762;
 0.24521339342874640, 0.65799309996997613, 0.35630912044469049;
 0.15663207261001361, 0.54283739749123061, 0.27953865212552687;
 0.06082276345468031, 0.45650136143553488, 0.23653979979309381;
 0.00000000000000000, 0.36962707428371205, 0.20039984916939454];
case 'YlGnBu' % In [11]: sns.color_palette("YlGnBu", 10)
cmap = [ ... 
 0.94906574698055490, 0.98019223493688246, 0.73779317210702333;
 0.86337563290315522, 0.94648212545058308, 0.69933104444952576;
 0.73388697750428145, 0.89564014462863695, 0.71040370814940512;
 0.52129181202720187, 0.81296425566953767, 0.73107268038917994;
 0.34262207781567294, 0.74626683557734774, 0.75589390151640945;
 0.20396771430969238, 0.66137641387827251, 0.76296810332466569;
 0.11534025493790122, 0.55215688698432031, 0.74519032660652607;
 0.13010381214758929, 0.40156863822656519, 0.67432527892729821;
 0.13988466630963720, 0.27690888619890397, 0.61514804293127623;
 0.11343329991487897, 0.17880815288015439, 0.51487891335113378];
case 'GnBu_d' % sns.color_palette("GnBu_d")
cmap = [ ...
 0.21697808798621680, 0.32733564601225013, 0.36941176807179171;
 0.23442778952760630, 0.45820839330261826, 0.54352941859002213;
 0.25140587751382310, 0.58554403931486831, 0.71294118666181383;
 0.32480841754308700, 0.68493145540648814, 0.78994746862673293;
 0.45066770474895150, 0.75099834881576832, 0.77038576275694604;
 0.58002308326608990, 0.81890043370863974, 0.75028067616855398];
case 'cubehelix_greenblue' % sns.cubehelix_palette(8, start=.5, rot=-.75)
cmap = [ ...
 0.84232988177938480, 0.87374044279641840, 0.75249540307310370;
 0.68251876357072430, 0.81069190728320800, 0.63524701801182060;
 0.51093657786460060, 0.73671986965753190, 0.57304087941263320;
 0.37208664465749840, 0.63786334195260290, 0.55503689058379240;
 0.28846276635927040, 0.51638144597481420, 0.54342177164221150;
 0.24670155725826660, 0.37340824813905654, 0.49725690696587516;
 0.22179626547231540, 0.23841378594571613, 0.39797674055755683;
 0.17250549177124480, 0.11951843162770594, 0.24320155229883056];
case 'cubehelix_purple' % sns.cubehelix_palette(8)
cmap = [ ... 
 0.93126922233253720, 0.82019217960821180, 0.79714809746635920;
 0.88228981687371890, 0.69582086670574200, 0.70654571194854310;
 0.81353802547006760, 0.57050551823578220, 0.63928085946815500;
 0.71958007083491190, 0.45537982893127477, 0.58610629958109260;
 0.60469068026344690, 0.35739308184976665, 0.53374078536924060;
 0.46496993672552040, 0.26868986121314253, 0.46365277636406470;
 0.32101947432593470, 0.19303051265196464, 0.37078816777247920;
 0.17508656489522050, 0.11840023306916837, 0.24215989137836502];
case 'red2blue' % sns.color_palette("RdBu_r", 7)
cmap = [ ... 
 0.16339870177063293, 0.44498270983789490, 0.6975009791991290;
 0.42068437209316328, 0.67643216077019186, 0.8186851319144753;
 0.76147636946509856, 0.86851211856393251, 0.9245674785445717;
 0.96908881383783674, 0.96647443490869855, 0.9649365649503820;
 0.98246828247519102, 0.80069205340217142, 0.7061130509657018;
 0.89457901435739851, 0.50380624217145586, 0.3997693394913390;
 0.72848905885920801, 0.15501730406985564, 0.1973856272650700];
case 'colorbrewer' % http://colorbrewer2.org/#type=diverging&scheme=Spectral&n=11
cmap = [ ...
 0.22656250000000000, 0.00390625000000000, 0.2578125000000000;
 0.83203125000000000, 0.24218750000000000, 0.3085937500000000;
 0.95312500000000000, 0.42578125000000000, 0.2617187500000000;
 0.98828125000000000, 0.67968750000000000, 0.3789062500000000;
 0.99218750000000000, 0.87500000000000000, 0.5429687500000000;
 0.99609375000000000, 0.99609375000000000, 0.7460937500000000;
 0.89843750000000000, 0.95703125000000000, 0.5937500000000000;
 0.66796875000000000, 0.86328125000000000, 0.6406250000000000;
 0.39843750000000000, 0.75781250000000000, 0.6445312500000000;
 0.19531250000000000, 0.53125000000000000, 0.7382812500000000;
 0.36718750000000000, 0.30859375000000000, 0.6328125000000000];
case 'magma'
cmap = [ ...
 1.46159096e-03,   4.66127766e-04,   1.38655200e-02;
 2.25764007e-03,   1.29495431e-03,   1.83311461e-02;
 3.27943222e-03,   2.30452991e-03,   2.37083291e-02;
 4.51230222e-03,   3.49037666e-03,   2.99647059e-02;
 5.94976987e-03,   4.84285000e-03,   3.71296695e-02;
 7.58798550e-03,   6.35613622e-03,   4.49730774e-02;
 9.42604390e-03,   8.02185006e-03,   5.28443561e-02;
 1.14654337e-02,   9.82831486e-03,   6.07496380e-02;
 1.37075706e-02,   1.17705913e-02,   6.86665843e-02;
 1.61557566e-02,   1.38404966e-02,   7.66026660e-02;
 1.88153670e-02,   1.60262753e-02,   8.45844897e-02;
 2.16919340e-02,   1.83201254e-02,   9.26101050e-02;
 2.47917814e-02,   2.07147875e-02,   1.00675555e-01;
 2.81228154e-02,   2.32009284e-02,   1.08786954e-01;
 3.16955304e-02,   2.57651161e-02,   1.16964722e-01;
 3.55204468e-02,   2.83974570e-02,   1.25209396e-01;
 3.96084872e-02,   3.10895652e-02,   1.33515085e-01;
 4.38295350e-02,   3.38299885e-02,   1.41886249e-01;
 4.80616391e-02,   3.66066101e-02,   1.50326989e-01;
 5.23204388e-02,   3.94066020e-02,   1.58841025e-01;
 5.66148978e-02,   4.21598925e-02,   1.67445592e-01;
 6.09493930e-02,   4.47944924e-02,   1.76128834e-01;
 6.53301801e-02,   4.73177796e-02,   1.84891506e-01;
 6.97637296e-02,   4.97264666e-02,   1.93735088e-01;
 7.42565152e-02,   5.20167766e-02,   2.02660374e-01;
 7.88150034e-02,   5.41844801e-02,   2.11667355e-01;
 8.34456313e-02,   5.62249365e-02,   2.20755099e-01;
 8.81547730e-02,   5.81331465e-02,   2.29921611e-01;
 9.29486914e-02,   5.99038167e-02,   2.39163669e-01;
 9.78334770e-02,   6.15314414e-02,   2.48476662e-01;
 1.02814972e-01,   6.30104053e-02,   2.57854400e-01;
 1.07898679e-01,   6.43351102e-02,   2.67288933e-01;
 1.13094451e-01,   6.54920358e-02,   2.76783978e-01;
 1.18405035e-01,   6.64791593e-02,   2.86320656e-01;
 1.23832651e-01,   6.72946449e-02,   2.95879431e-01;
 1.29380192e-01,   6.79349264e-02,   3.05442931e-01;
 1.35053322e-01,   6.83912798e-02,   3.14999890e-01;
 1.40857952e-01,   6.86540710e-02,   3.24537640e-01;
 1.46785234e-01,   6.87382323e-02,   3.34011109e-01;
 1.52839217e-01,   6.86368599e-02,   3.43404450e-01;
 1.59017511e-01,   6.83540225e-02,   3.52688028e-01;
 1.65308131e-01,   6.79108689e-02,   3.61816426e-01;
 1.71713033e-01,   6.73053260e-02,   3.70770827e-01;
 1.78211730e-01,   6.65758073e-02,   3.79497161e-01;
 1.84800877e-01,   6.57324381e-02,   3.87972507e-01;
 1.91459745e-01,   6.48183312e-02,   3.96151969e-01;
 1.98176877e-01,   6.38624166e-02,   4.04008953e-01;
 2.04934882e-01,   6.29066192e-02,   4.11514273e-01;
 2.11718061e-01,   6.19917876e-02,   4.18646741e-01;
 2.18511590e-01,   6.11584918e-02,   4.25391816e-01;
 2.25302032e-01,   6.04451843e-02,   4.31741767e-01;
 2.32076515e-01,   5.98886855e-02,   4.37694665e-01;
 2.38825991e-01,   5.95170384e-02,   4.43255999e-01;
 2.45543175e-01,   5.93524384e-02,   4.48435938e-01;
 2.52220252e-01,   5.94147119e-02,   4.53247729e-01;
 2.58857304e-01,   5.97055998e-02,   4.57709924e-01;
 2.65446744e-01,   6.02368754e-02,   4.61840297e-01;
 2.71994089e-01,   6.09935552e-02,   4.65660375e-01;
 2.78493300e-01,   6.19778136e-02,   4.69190328e-01;
 2.84951097e-01,   6.31676261e-02,   4.72450879e-01;
 2.91365817e-01,   6.45534486e-02,   4.75462193e-01;
 2.97740413e-01,   6.61170432e-02,   4.78243482e-01;
 3.04080941e-01,   6.78353452e-02,   4.80811572e-01;
 3.10382027e-01,   6.97024767e-02,   4.83186340e-01;
 3.16654235e-01,   7.16895272e-02,   4.85380429e-01;
 3.22899126e-01,   7.37819504e-02,   4.87408399e-01;
 3.29114038e-01,   7.59715081e-02,   4.89286796e-01;
 3.35307503e-01,   7.82361045e-02,   4.91024144e-01;
 3.41481725e-01,   8.05635079e-02,   4.92631321e-01;
 3.47635742e-01,   8.29463512e-02,   4.94120923e-01;
 3.53773161e-01,   8.53726329e-02,   4.95501096e-01;
 3.59897941e-01,   8.78311772e-02,   4.96778331e-01;
 3.66011928e-01,   9.03143031e-02,   4.97959963e-01;
 3.72116205e-01,   9.28159917e-02,   4.99053326e-01;
 3.78210547e-01,   9.53322947e-02,   5.00066568e-01;
 3.84299445e-01,   9.78549106e-02,   5.01001964e-01;
 3.90384361e-01,   1.00379466e-01,   5.01864236e-01;
 3.96466670e-01,   1.02902194e-01,   5.02657590e-01;
 4.02547663e-01,   1.05419865e-01,   5.03385761e-01;
 4.08628505e-01,   1.07929771e-01,   5.04052118e-01;
 4.14708664e-01,   1.10431177e-01,   5.04661843e-01;
 4.20791157e-01,   1.12920210e-01,   5.05214935e-01;
 4.26876965e-01,   1.15395258e-01,   5.05713602e-01;
 4.32967001e-01,   1.17854987e-01,   5.06159754e-01;
 4.39062114e-01,   1.20298314e-01,   5.06555026e-01;
 4.45163096e-01,   1.22724371e-01,   5.06900806e-01;
 4.51270678e-01,   1.25132484e-01,   5.07198258e-01;
 4.57385535e-01,   1.27522145e-01,   5.07448336e-01;
 4.63508291e-01,   1.29892998e-01,   5.07651812e-01;
 4.69639514e-01,   1.32244819e-01,   5.07809282e-01;
 4.75779723e-01,   1.34577500e-01,   5.07921193e-01;
 4.81928997e-01,   1.36891390e-01,   5.07988509e-01;
 4.88088169e-01,   1.39186217e-01,   5.08010737e-01;
 4.94257673e-01,   1.41462106e-01,   5.07987836e-01;
 5.00437834e-01,   1.43719323e-01,   5.07919772e-01;
 5.06628929e-01,   1.45958202e-01,   5.07806420e-01;
 5.12831195e-01,   1.48179144e-01,   5.07647570e-01;
 5.19044825e-01,   1.50382611e-01,   5.07442938e-01;
 5.25269968e-01,   1.52569121e-01,   5.07192172e-01;
 5.31506735e-01,   1.54739247e-01,   5.06894860e-01;
 5.37755194e-01,   1.56893613e-01,   5.06550538e-01;
 5.44015371e-01,   1.59032895e-01,   5.06158696e-01;
 5.50287252e-01,   1.61157816e-01,   5.05718782e-01;
 5.56570783e-01,   1.63269149e-01,   5.05230210e-01;
 5.62865867e-01,   1.65367714e-01,   5.04692365e-01;
 5.69172368e-01,   1.67454379e-01,   5.04104606e-01;
 5.75490107e-01,   1.69530062e-01,   5.03466273e-01;
 5.81818864e-01,   1.71595728e-01,   5.02776690e-01;
 5.88158375e-01,   1.73652392e-01,   5.02035167e-01;
 5.94508337e-01,   1.75701122e-01,   5.01241011e-01;
 6.00868399e-01,   1.77743036e-01,   5.00393522e-01;
 6.07238169e-01,   1.79779309e-01,   4.99491999e-01;
 6.13617209e-01,   1.81811170e-01,   4.98535746e-01;
 6.20005032e-01,   1.83839907e-01,   4.97524075e-01;
 6.26401108e-01,   1.85866869e-01,   4.96456304e-01;
 6.32804854e-01,   1.87893468e-01,   4.95331769e-01;
 6.39215638e-01,   1.89921182e-01,   4.94149821e-01;
 6.45632778e-01,   1.91951556e-01,   4.92909832e-01;
 6.52055535e-01,   1.93986210e-01,   4.91611196e-01;
 6.58483116e-01,   1.96026835e-01,   4.90253338e-01;
 6.64914668e-01,   1.98075202e-01,   4.88835712e-01;
 6.71349279e-01,   2.00133166e-01,   4.87357807e-01;
 6.77785975e-01,   2.02202663e-01,   4.85819154e-01;
 6.84223712e-01,   2.04285721e-01,   4.84219325e-01;
 6.90661380e-01,   2.06384461e-01,   4.82557941e-01;
 6.97097796e-01,   2.08501100e-01,   4.80834678e-01;
 7.03531700e-01,   2.10637956e-01,   4.79049270e-01;
 7.09961888e-01,   2.12797337e-01,   4.77201121e-01;
 7.16387038e-01,   2.14981693e-01,   4.75289780e-01;
 7.22805451e-01,   2.17193831e-01,   4.73315708e-01;
 7.29215521e-01,   2.19436516e-01,   4.71278924e-01;
 7.35615545e-01,   2.21712634e-01,   4.69179541e-01;
 7.42003713e-01,   2.24025196e-01,   4.67017774e-01;
 7.48378107e-01,   2.26377345e-01,   4.64793954e-01;
 7.54736692e-01,   2.28772352e-01,   4.62508534e-01;
 7.61077312e-01,   2.31213625e-01,   4.60162106e-01;
 7.67397681e-01,   2.33704708e-01,   4.57755411e-01;
 7.73695380e-01,   2.36249283e-01,   4.55289354e-01;
 7.79967847e-01,   2.38851170e-01,   4.52765022e-01;
 7.86212372e-01,   2.41514325e-01,   4.50183695e-01;
 7.92426972e-01,   2.44242250e-01,   4.47543155e-01;
 7.98607760e-01,   2.47039798e-01,   4.44848441e-01;
 8.04751511e-01,   2.49911350e-01,   4.42101615e-01;
 8.10854841e-01,   2.52861399e-01,   4.39304963e-01;
 8.16914186e-01,   2.55894550e-01,   4.36461074e-01;
 8.22925797e-01,   2.59015505e-01,   4.33572874e-01;
 8.28885740e-01,   2.62229049e-01,   4.30643647e-01;
 8.34790818e-01,   2.65539703e-01,   4.27671352e-01;
 8.40635680e-01,   2.68952874e-01,   4.24665620e-01;
 8.46415804e-01,   2.72473491e-01,   4.21631064e-01;
 8.52126490e-01,   2.76106469e-01,   4.18572767e-01;
 8.57762870e-01,   2.79856666e-01,   4.15496319e-01;
 8.63320397e-01,   2.83729003e-01,   4.12402889e-01;
 8.68793368e-01,   2.87728205e-01,   4.09303002e-01;
 8.74176342e-01,   2.91858679e-01,   4.06205397e-01;
 8.79463944e-01,   2.96124596e-01,   4.03118034e-01;
 8.84650824e-01,   3.00530090e-01,   4.00047060e-01;
 8.89731418e-01,   3.05078817e-01,   3.97001559e-01;
 8.94700194e-01,   3.09773445e-01,   3.93994634e-01;
 8.99551884e-01,   3.14616425e-01,   3.91036674e-01;
 9.04281297e-01,   3.19609981e-01,   3.88136889e-01;
 9.08883524e-01,   3.24755126e-01,   3.85308008e-01;
 9.13354091e-01,   3.30051947e-01,   3.82563414e-01;
 9.17688852e-01,   3.35500068e-01,   3.79915138e-01;
 9.21884187e-01,   3.41098112e-01,   3.77375977e-01;
 9.25937102e-01,   3.46843685e-01,   3.74959077e-01;
 9.29845090e-01,   3.52733817e-01,   3.72676513e-01;
 9.33606454e-01,   3.58764377e-01,   3.70540883e-01;
 9.37220874e-01,   3.64929312e-01,   3.68566525e-01;
 9.40687443e-01,   3.71224168e-01,   3.66761699e-01;
 9.44006448e-01,   3.77642889e-01,   3.65136328e-01;
 9.47179528e-01,   3.84177874e-01,   3.63701130e-01;
 9.50210150e-01,   3.90819546e-01,   3.62467694e-01;
 9.53099077e-01,   3.97562894e-01,   3.61438431e-01;
 9.55849237e-01,   4.04400213e-01,   3.60619076e-01;
 9.58464079e-01,   4.11323666e-01,   3.60014232e-01;
 9.60949221e-01,   4.18323245e-01,   3.59629789e-01;
 9.63310281e-01,   4.25389724e-01,   3.59469020e-01;
 9.65549351e-01,   4.32518707e-01,   3.59529151e-01;
 9.67671128e-01,   4.39702976e-01,   3.59810172e-01;
 9.69680441e-01,   4.46935635e-01,   3.60311120e-01;
 9.71582181e-01,   4.54210170e-01,   3.61030156e-01;
 9.73381238e-01,   4.61520484e-01,   3.61964652e-01;
 9.75082439e-01,   4.68860936e-01,   3.63111292e-01;
 9.76690494e-01,   4.76226350e-01,   3.64466162e-01;
 9.78209957e-01,   4.83612031e-01,   3.66024854e-01;
 9.79645181e-01,   4.91013764e-01,   3.67782559e-01;
 9.81000291e-01,   4.98427800e-01,   3.69734157e-01;
 9.82279159e-01,   5.05850848e-01,   3.71874301e-01;
 9.83485387e-01,   5.13280054e-01,   3.74197501e-01;
 9.84622298e-01,   5.20712972e-01,   3.76698186e-01;
 9.85692925e-01,   5.28147545e-01,   3.79370774e-01;
 9.86700017e-01,   5.35582070e-01,   3.82209724e-01;
 9.87646038e-01,   5.43015173e-01,   3.85209578e-01;
 9.88533173e-01,   5.50445778e-01,   3.88365009e-01;
 9.89363341e-01,   5.57873075e-01,   3.91670846e-01;
 9.90138201e-01,   5.65296495e-01,   3.95122099e-01;
 9.90871208e-01,   5.72706259e-01,   3.98713971e-01;
 9.91558165e-01,   5.80106828e-01,   4.02441058e-01;
 9.92195728e-01,   5.87501706e-01,   4.06298792e-01;
 9.92784669e-01,   5.94891088e-01,   4.10282976e-01;
 9.93325561e-01,   6.02275297e-01,   4.14389658e-01;
 9.93834412e-01,   6.09643540e-01,   4.18613221e-01;
 9.94308514e-01,   6.16998953e-01,   4.22949672e-01;
 9.94737698e-01,   6.24349657e-01,   4.27396771e-01;
 9.95121854e-01,   6.31696376e-01,   4.31951492e-01;
 9.95480469e-01,   6.39026596e-01,   4.36607159e-01;
 9.95809924e-01,   6.46343897e-01,   4.41360951e-01;
 9.96095703e-01,   6.53658756e-01,   4.46213021e-01;
 9.96341406e-01,   6.60969379e-01,   4.51160201e-01;
 9.96579803e-01,   6.68255621e-01,   4.56191814e-01;
 9.96774784e-01,   6.75541484e-01,   4.61314158e-01;
 9.96925427e-01,   6.82827953e-01,   4.66525689e-01;
 9.97077185e-01,   6.90087897e-01,   4.71811461e-01;
 9.97186253e-01,   6.97348991e-01,   4.77181727e-01;
 9.97253982e-01,   7.04610791e-01,   4.82634651e-01;
 9.97325180e-01,   7.11847714e-01,   4.88154375e-01;
 9.97350983e-01,   7.19089119e-01,   4.93754665e-01;
 9.97350583e-01,   7.26324415e-01,   4.99427972e-01;
 9.97341259e-01,   7.33544671e-01,   5.05166839e-01;
 9.97284689e-01,   7.40771893e-01,   5.10983331e-01;
 9.97228367e-01,   7.47980563e-01,   5.16859378e-01;
 9.97138480e-01,   7.55189852e-01,   5.22805996e-01;
 9.97019342e-01,   7.62397883e-01,   5.28820775e-01;
 9.96898254e-01,   7.69590975e-01,   5.34892341e-01;
 9.96726862e-01,   7.76794860e-01,   5.41038571e-01;
 9.96570645e-01,   7.83976508e-01,   5.47232992e-01;
 9.96369065e-01,   7.91167346e-01,   5.53498939e-01;
 9.96162309e-01,   7.98347709e-01,   5.59819643e-01;
 9.95932448e-01,   8.05527126e-01,   5.66201824e-01;
 9.95680107e-01,   8.12705773e-01,   5.72644795e-01;
 9.95423973e-01,   8.19875302e-01,   5.79140130e-01;
 9.95131288e-01,   8.27051773e-01,   5.85701463e-01;
 9.94851089e-01,   8.34212826e-01,   5.92307093e-01;
 9.94523666e-01,   8.41386618e-01,   5.98982818e-01;
 9.94221900e-01,   8.48540474e-01,   6.05695903e-01;
 9.93865767e-01,   8.55711038e-01,   6.12481798e-01;
 9.93545285e-01,   8.62858846e-01,   6.19299300e-01;
 9.93169558e-01,   8.70024467e-01,   6.26189463e-01;
 9.92830963e-01,   8.77168404e-01,   6.33109148e-01;
 9.92439881e-01,   8.84329694e-01,   6.40099465e-01;
 9.92089454e-01,   8.91469549e-01,   6.47116021e-01;
 9.91687744e-01,   8.98627050e-01,   6.54201544e-01;
 9.91331929e-01,   9.05762748e-01,   6.61308839e-01;
 9.90929685e-01,   9.12915010e-01,   6.68481201e-01;
 9.90569914e-01,   9.20048699e-01,   6.75674592e-01;
 9.90174637e-01,   9.27195612e-01,   6.82925602e-01;
 9.89814839e-01,   9.34328540e-01,   6.90198194e-01;
 9.89433736e-01,   9.41470354e-01,   6.97518628e-01;
 9.89077438e-01,   9.48604077e-01,   7.04862519e-01;
 9.88717064e-01,   9.55741520e-01,   7.12242232e-01;
 9.88367028e-01,   9.62878026e-01,   7.19648627e-01;
 9.88032885e-01,   9.70012413e-01,   7.27076773e-01;
 9.87690702e-01,   9.77154231e-01,   7.34536205e-01;
 9.87386827e-01,   9.84287561e-01,   7.42001547e-01;
 9.87052509e-01,   9.91437853e-01,   7.49504188e-01];
end

% interpolating function
map_ = @(n,cmap) interp1([0:(size(cmap,1)-1)]./(size(cmap,1)-1),cmap,linspace(0,1,n));

% interpolate and return
cmap = map_(n,cmap);

end





