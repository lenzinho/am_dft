clear;clc;

read_ = @(fname) fscanf(fopen(fname),'%f');
th = read_('th.itx');
phi = read_('phi.itx');
data = read_('data.itx');

data = reshape(data,length(phi),length(th));
% surf(data,'edgecolor','none');

E_photon = 137.85;
E_binding = 142;
E_inner = 10;
hbar = 1; m = 1;

[phi,th,E] = meshgrid(th,phi,E_photon);

x = 1/hbar*sqrt(2*m*E).*sind(phi).*cosd(th);
y = 1/hbar*sqrt(2*m*E).*sind(phi).*sind(th);
z = 1/hbar*sqrt(2*m*(E.*cosd(phi).^2+E_inner));  

surf(x,y,data,'edgecolor','none');
axis tight; box on; grid off; daspect([1 1 1]);
view([0 0 1])

caxis([6e-6 1.5e-05])

% symmeterize
hold on;
phi = phi + 180;
x = 1/hbar*sqrt(2*m*E).*sind(phi).*cosd(th);
y = 1/hbar*sqrt(2*m*E).*sind(phi).*sind(th);
z = 1/hbar*sqrt(2*m*(E.*cosd(phi).^2+E_inner));  
surf(x,y,data,'edgecolor','none');
hold off;

