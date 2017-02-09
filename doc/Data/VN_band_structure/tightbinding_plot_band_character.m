clear;clc;

% set tight binding matrix elements
v = [-0.4212,1.1975,-4.1841,-1.0193,-1.0322,-0.0565,0.1132,-0.5218,-0.1680,0.0635,-0.0546,-0.1051,0.4189,0.3061];

% define kpoint path (cart, units of 2pi/a)
iM=ones(3)-2*eye(3); G=[0;0;0]; X1=[0;1;1]/2; X2=[2;1;1]/2; L=[1;1;1]/2; K=[6;3;3]/8;
N = 100; ql={'G','X','K','G','L'}; qs=iM*[G,X2,K,G]; qe=iM*[X1,K,G,L]; nqs = size(qs,2);  

% build path
[k,x,qt] = get_path(qs,qe,nqs,N); nks = numel(x); 

% normalize columns of matrix
normc_ = @(m) ones(size(m,1),1)*sqrt(1./sum(m.*m)).*m;

% get eigenvalues and band character weights
nbands = 8; degen = rand(3,1)*1E-5; c = zeros(nbands,nks); E = zeros(nbands,nks);
for i = 1:nks
    % compute energies
    H = getH(v,k(1:3,i)+degen); [V,S]=eig(H,'vector'); [E(:,i),ind]=sort(real(S));
    
    % compute weights for band character
    Vp(1,:) = sum(abs(V([1,2,3],:)),1); % p
    Vp(2,:) = sum(abs(V([4,5,7],:)),1); % d_t2g
    Vp(3,:) = sum(abs(V([6,8]  ,:)),1); % d_eg
    Vp = normc_(Vp); c(:,i)=assign_color( Vp(:,ind) );
end

% define figure properties
fig_ = @(h)       set(h,'color','white');
axs_ = @(h,qt,ql) set(h,'Box','on','XTick',qt,'Xticklabel',ql,'ylim',[-10 6]);

% plot band structure  (colorshift allows rotation of cyclic colormap without affecting data)
figure(2); clf; fig_(gcf); colorshift = 25;
hold on; for j = 1:8; plotc(x,E(j,:),mod(c(j,:)+colorshift/100,1)); end; hold off; axis tight;
axs_(gca,qt,ql); ylabel('Energy E'); xlabel('Wavevector k');

% plot fermi level
line([0,x(end)],[0,0],'linewidth',2,'color',[1,1,1]*0.5,'linestyle',':');

% colar pallette borrowed from python's seaborn library (added last row for periodicity)
% import seaborn as sns; sns.color_palette("Spectral", 10)
spectral=[ ...
 0.81414841553744144, 0.21968473703143937, 0.30480585554066825;
 0.93302576331531295, 0.39131103777417953, 0.27197233193060932;
 0.98177624099394856, 0.60738179087638855, 0.34579008992980509;
 0.99469434864380779, 0.80922723167082844, 0.48696657138712268;
 0.99823144954793597, 0.94517493598601399, 0.65705499929540301;
 0.95578623869839840, 0.98231449547935934, 0.68004615517223588;
 0.82029989537070780, 0.92756632496328917, 0.61268745099796973;
 0.59100347582031698, 0.83552480795804196, 0.64429067864137535;
 0.36001538412243711, 0.71618609919267540, 0.66551328406614418;
 0.21299500558890549, 0.51141871132102668, 0.73079586379668293;
 0.81414841553744144, 0.21968473703143937, 0.30480585554066825];
map_ = @(n,cmap) interp1([0:(size(cmap,1)-1)]./(size(cmap,1)-1),[cmap],linspace(0,1,n));

% apply color map and label axes
colormap(flipud(circshift(map_(100,spectral),-colorshift+10)).^(1.25)); 
h = colorbar; caxis([0,1]); set(h,'Ticks',sort(mod([0.33,0.66,1]-10/100+colorshift/100,1)),'TickLabels',{'p','d_{t2g}','d_{eg}'});

% -------------------------------------------------------------------------
% plot DFT results (requires EIGENVAL)
[dft,bz]=load_vasp_eigenval('EIGENVAL.nscf'); 

% set Fermi energy manually (see OUTCAR)
EF = 9.50006808+3.2; dft.E = dft.E-EF;

% define DFT bands to skip;
skip_nbands = 5;

% rebuild path with number of kpoints matching KPOINTS
N = 40; [k,x,qt] = get_path(qs,qe,nqs,N); nks = numel(x); 

% plot dft bands
inds = [1:5:nks];
hold on; h_dft = plot(x(inds),dft.E(skip_nbands+[1:size(E,1)],inds),'.','color',[.95,.95,1]*0.82,'markersize',7); hold off;
uistack(h_dft,'bottom')
% -------------------------------------------------------------------------
% save figure
set(gca,'LooseInset',get(gca,'TightInset')); set(gcf,'PaperSize',[10 10]);
set(gcf,'PaperPosition',[0,0,1.6180,1]*3); set(gca,'Ytick',[-10:2:6]); 
print(gcf,'-dpdf','band');
%%

function H = getH(tb,k)
% parse wavevector and 14 tight binding parameters
k=num2cell(k);[k1,k2,k3]=deal(k{:}); tb=num2cell(tb);[a11,a12,b11,c11,c12,d11,d12,d13,d14,d15,d16,e11,e12,e13]=deal(tb{:});

% trigonometric abbreviations
sx = sin(pi*k1); sy = sin(pi*k2); sz = sin(pi*k3); sxy=sx*sy; sxz=sx*sz; syz=sy*sz;
cx = cos(pi*k1); cy = cos(pi*k2); cz = cos(pi*k3); cxy=cx*cy; cxz=cx*cz; cyz=cy*cz;

% build dynamical matrix
H = reshape([b11+cxy.*e13.*4.0+cxz.*e13.*4.0+cyz.*e11.*4.0,e12.*sxy.*-4.0,e12.*sxz.*-4.0,c11.*sz.*2.0i,c11.*sy.*-2.0i,c12.*sx.*2.0i,0.0,sqrt(3.0).*c12.*sx.*-2.0i,e12.*sxy.*-4.0,b11+cxy.*e13.*4.0+cxz.*e11.*4.0+cyz.*e13.*4.0,e12.*syz.*-4.0,0.0,c11.*sx.*-2.0i,c12.*sy.*-4.0i,c11.*sz.*-2.0i,0.0,e12.*sxz.*-4.0,e12.*syz.*-4.0,b11+cxy.*e11.*4.0+cxz.*e13.*4.0+cyz.*e13.*4.0,c11.*sx.*2.0i,0.0,c12.*sz.*2.0i,c11.*sy.*-2.0i,sqrt(3.0).*c12.*sz.*2.0i,c11.*sz.*-2.0i,0.0,c11.*sx.*-2.0i,a11+cxy.*d12.*4.0+cxz.*d13.*4.0+cyz.*d12.*4.0,d11.*syz.*-4.0,d15.*sxz.*-8.0,d11.*sxy.*-4.0,0.0,c11.*sy.*2.0i,c11.*sx.*2.0i,0.0,d11.*syz.*-4.0,a11+cxy.*d13.*4.0+cxz.*d12.*4.0+cyz.*d12.*4.0,d15.*sxy.*-4.0,d11.*sxz.*4.0,sqrt(3.0).*d15.*sxy.*-4.0,c12.*sx.*-2.0i,c12.*sy.*4.0i,c12.*sz.*-2.0i,d15.*sxz.*-8.0,d15.*sxy.*-4.0,a12+cxy.*d14.*4.0+cxy.*d16.*4.0-cxz.*d14.*2.0+cxz.*d16.*4.0+cyz.*d14.*4.0+cyz.*d16.*4.0,d15.*syz.*-4.0,sqrt(3.0).*cy.*d14.*(cx-cz).*-2.0,0.0,c11.*sz.*2.0i,c11.*sy.*2.0i,d11.*sxy.*-4.0,d11.*sxz.*4.0,d15.*syz.*-4.0,a11+cxy.*d12.*4.0+cxz.*d12.*4.0+cyz.*d13.*4.0,sqrt(3.0).*d15.*syz.*4.0,sqrt(3.0).*c12.*sx.*2.0i,0.0,sqrt(3.0).*c12.*sz.*-2.0i,0.0,sqrt(3.0).*d15.*sxy.*-4.0,sqrt(3.0).*cy.*d14.*(cx-cz).*-2.0,sqrt(3.0).*d15.*syz.*4.0,a12+cxy.*d16.*4.0+cxz.*d14.*6.0+cxz.*d16.*4.0+cyz.*d16.*4.0],[8,8]);
end

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

function [h] = plotc(x,y,c)
    x = x(:).'; y=y(:).'; c=c(:).'; z=zeros(size(x));
    % col = x;  % This is the color, vary with x in this case.
    surface([x;x],[y;y],[z;z],[c;c],'facecol','no','edgecol','interp','linew',1);
end

function [th] = assign_color(V)
% assigns a number betwee [0,1] based on how close the vectors V are to the identity vectors.

% set number of points
n = size(V,1);

% define points around circle on complex plane
p = exp(2*pi*1i*[1:n]/n);

% get complex points
c = (p*abs(V)).';

% get linear mapping between [0,1]
th = (atan2d(imag(c),real(c))/180+1)/2;
end

function [dr,bz] = load_vasp_eigenval(fname)
    fprintf('loading dispersion from: %s \n',fname);
    fid=fopen(fname);
    % skip first five lines
    for i = 1:5
    fgetl(fid);
    end
    buffer = strsplit(strtrim(fgetl(fid)));
    dr.nelecs = sscanf(buffer{1},'%i');
    bz.nkpts  = sscanf(buffer{2},'%i');
    dr.nbands = sscanf(buffer{3},'%i');
    fprintf(' ... electrons = %i \n',dr.nelecs);
    fprintf(' ... kpoints = %i \n',bz.nkpts);
    fprintf(' ... bands = %i \n',dr.nbands);
    for i = 1:bz.nkpts
        % skip line
        fgetl(fid);
        % get kpnts
        buffer = strsplit(strtrim(fgetl(fid)));
        bz.kpt(1,i) = sscanf(buffer{1},'%f');
        bz.kpt(2,i) = sscanf(buffer{2},'%f');
        bz.kpt(3,i) = sscanf(buffer{3},'%f');
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

