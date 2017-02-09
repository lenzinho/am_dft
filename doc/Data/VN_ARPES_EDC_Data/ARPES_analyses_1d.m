clear;clc;

% define figure properties
fig = @(h) set(h,'color','white');
axs = @(h) set(h,'Box','on','ylim',[4 0.2]);

% set primitive basis [nm]
bas =  [0.000000  0.500000  0.500000; 
        0.500000 -0.000000  0.500000;
        0.500000  0.500000 -0.000000]*4.071153;
recbas = inv(bas);

% define in-plane and out-of-plane directions
c = (ones(3)-eye(3))/2; % hkl => frac
para = (c)*[+1,+0,+0]';
perp = (c)*[+0,+1,+1]';
para=para./norm(perp);
perp=perp./norm(perp);

% define fermi energy photon energy, relative electron mass (m/m_e),
% analyzer work function [eV], kinetic energy [eV]; Note: fermi level
% defined above. 
hv = 146; Ef = 142; WF_analyzer = 4.2; m = 1; V0 = 12; 

% set tight binding matrix elements
v = [  1.1584087
      -0.4766591
      -4.0853000
       1.0456333
       1.0382811
      -0.1302027
      -0.2203503
      -0.0743411
      -0.5334793
       0.0554776
       0.1127403
      -0.0910807
       0.4326968
       0.3172856];

% load ARPES data
be = importdata('energy.txt');
th = importdata('angle.txt');  
I  = importdata('spectra.txt');
I  = I.data(:,2:end);
th = th(:,2); 

% down sample data "dn" times
dn = 20;
I  = I(1:dn:end,1:dn:end);
th = th(1:dn:end);
be = be(1:dn:end)-Ef;
ke = hv-(Ef+be)-WF_analyzer;

% exclude negative values
% ex_ = (ke>0);
% be = be(ex_);
% ke = ke(ex_);
% I = I(ex_,:);

[  ~,be2] = meshgrid(th,be);
[th2,ke2] = meshgrid(th,ke);

% % plot ARPES
% figure(1); fig(gcf);
% subplot(2,1,1);
% surf(th2,be2,I-1E8,'edgecolor','none'); 
% view([0 0 1]); axis tight;

% dimensions of data
nkes = size(ke2,1);
nths = size(th2,2);
nbands = 3+5;

% load point group matrices
tb  = xml_read('../save.tb');
sym = reshape(str2num(tb.group.sym.value),tb.group.sym.shape);

% Define final wavevector calculation based on free-electron final state
% In-plane electron wavevector is conserved; out-of-plane wavevector isn't
% [Eq. 7.2, 7.3 and p 535, Hunfer]
a = 3.622627; % sqrt(eV*m_e)/hbar/nm
kf_para_ = @(ke,th,m)    a*sqrt(2*m*ke).*sind(th);
kf_perp_ = @(ke,th,m,V0) a*sqrt(2*m*ke .*cosd(th).^2+V0);

% determine wavevector [nm => frac]
nks = length(th);

% define gaussian function
degauss = 1;
gauss_ = @(x) exp(-abs(x).^2)./sqrt(pi);

% get wavevectors
k = zeros(3,nths);
for i = 1:3
    k(i,:) = k(i,:)+perp(i)*kf_perp_(ke(1),th,m,V0)';
    k(i,:) = k(i,:)+para(i)*kf_para_(ke(1),th,m)';
end
k = recbas*k;

% get eigenvalues
E = zeros(nbands,nths);
for ith = 1:nths
    H = get_H_numeric_frac(sym, v, k(1:3,ith), [1:100]);
    [V(:,:,ith),E(:,ith)] = eig(H,'vector');
end

%  get amplitude
for ith = 1:nths
    x = repmat(E(:,ith),1,nkes) - repmat(permute(be,[2,1]),nbands,1);
%     w = repmat( abs(V(3,:,ith)).^2, nkes, 1)';
% 	A(:,ith) = sum( w .*  gauss_(x/degauss)/degauss, 1);
	A(:,ith) = sum( gauss_(x/degauss)/degauss, 1);
end

% plot ARPES
figure(1); fig(gcf);
subplot(2,1,1);
surf(th2,be2,I-1E8,'edgecolor','none'); 
hold on; plot(th,sort(real(E),1),'color',[1 1 1],'linewidth',1.5); hold off; box on;
view([0 0 1]); axis tight; ylim([-4 0]);

subplot(2,1,2);
surf(th2,be2,A,'edgecolor','none'); 
view([0 0 1]); axis tight; box on;






