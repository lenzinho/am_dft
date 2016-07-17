% fresh start
clear;clc;
tiny=1.0E-6;
% load symmetry rep
dumpdir = '../../../debug/tbpg/';
load_dumps
% load pc
[pc]  = load_poscar('../../uc/outfile.primitive');
% load optimizd matrix elements
load VN_pd_optimized_tb_parameters.mat

%
bas    = (ones(3)-eye(3))/2;
rot    = vrrotvec2mat([0 0 1 pi/2]);
recbas = inv(bas);
m = 50; n = 50; o = 50;
hh     = linspace(-1,1,m);
kk     = linspace(-1,1,n);
ll     = linspace(-1,1,o);
[h,k,l]= meshgrid(hh,kk,ll);
z = 0;
for i = 1:m
for ii = 1:n
for iii = 1:o
    z = z+1;
    bz.kpt(:,z) = ([h(i,ii,iii);k(i,ii,iii);l(i,ii,iii)])/4.0712;
    bz.i(z) = i;
    bz.j(z) = ii;
    bz.k(z) = iii;
end
end
end
bz.nkpts = z;

D(1:m,1:n,1:o,8) = 0;
for z = 1:bz.nkpts
    H = get_H_numeric_cart(tbpg, xmin, bz.kpt(:,z), [1:100]);
    D(bz.i(z),bz.j(z),bz.k(z),:) = sort(real(eig(H)));
end

%% determine what is inside and outside of FBZ
% make FBZ just a little bigger for nice effects
for z = 1:bz.nkpts
    bool(bz.i(z),bz.j(z),bz.k(z))=~isinfbz(bas*4.0712*1.001,bz.kpt(:,z));
end
for z = 1:bz.nkpts
    if (bool(bz.i(z),bz.j(z),bz.k(z))) 
        boolnan(bz.i(z),bz.j(z),bz.k(z))=1;
    else
        boolnan(bz.i(z),bz.j(z),bz.k(z))=NaN;
    end
end

%% draw BZ boundary
clc;clf;
figure(1); set(gcf,'color','w'); hold on;
view(3); daspect([1 1 1]); axis tight; box on; grid on
fv = isosurface(X,Y,Z,bool);
p2 = patch(fv,'FaceColor','blue','facealpha',1,'edgecolor','k');
reducepatch(p2,0.1)
% isonormals(V,p2)

%% get V "green's function"
for z = 1:bz.nkpts
    V(bz.i(z),bz.j(z),bz.k(z)) = sum( 1./(D(bz.i(z),bz.j(z),bz.k(z),:).^2+(0.1).^2));
end
%%
clc
hsp = surf(linspace(-1,1,100),linspace(-1,1,100),zeros(100));
% rotate(hsp,[0,1,0],-45)
rotate(hsp,[0,1,0],-45)
xd = hsp.XData;
yd = hsp.YData;
zd = hsp.ZData;
handle = slice(h,k,l,V.*boolnan,xd,yd,zd);
% handle = slice(h,k,l,V,h(:,:,10),h(:,:,10),l(:,:,10));
daspect([1 1 1])
view([0 0 1]);axis tight;
handle.EdgeColor='none';

%% PLOT FERMI SURFACE MOVIE
% clc;clf; 
[h,k,l]= meshgrid(hh,kk,ll);

subplot(1,2,1); 
set(gcf,'color','w'); hold on;
view(3); daspect([1 1 1]); box on; grid on; axis off; axis tight;
arrow([0,0,0],[ 1, 1, 1],'linewidth',2)
arrow([0,0,0],[ 1,-1,-1],'linewidth',2)
arrow([0,0,0],[-1, 1,-1],'linewidth',2)
set(gca,'cameraviewanglemode','manual');

% arrow FIXLIMITS
% draw BZ boundary
fv = isosurface(h,k,l,bool);
p2 = patch(fv,'FaceColor',[1,1,1]*0,'facealpha',0.05,'edgecolor','none');
isonormals(V,p2)
% loop over fermi surfaces
Eflist = linspace(min(D(:))+0.1,max(D(:))-0.1,3000);
for j = 1%1:length(Eflist);
    Ef = 0;%Eflist(j);
    s = 0;
    for i = 1:size(D,4)
        V = D(:,:,:,i).*boolnan;
        if min(V(:))<Ef
        if max(V(:))>Ef
        s=s+1;
        [faces,verts,colors] = isosurface(h,k,l,V,Ef,h.^2+k.^2+l.^2); 
        handle(s) = patch('Vertices', verts, 'Faces', faces, ... 
              'FaceVertexCData', colors, ... 
              'FaceColor','interp', ... 
              'edgecolor','interp');
        end
        end
    end
    caxis([0 1])
    camorbit(0.3, 0)
    camorbit(0.1, 0.1)
    
    % plot EF lin
    subplot(1,2,2); 
    s=s+1;
    handle(s) = line([0,0.0093],[Ef,Ef]);
    set(handle(s),'linewidth',1.5)
    title(sprintf('E = %5.2f eV',Ef));
    subplot(1,2,1); 
    
    drawnow;
    pause(0.1)
    F(j) = getframe(gcf);
    delete(handle(1:s))
end

%%
fig = figure;
movie(fig,F,1)

%%

v = VideoWriter('myFile.superfast','MPEG-4');
v.FrameRate = 120;
open(v)
writeVideo(v,F)
close(v)

