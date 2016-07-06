% plot

clear;clc;
tiny=1.0E-6;
dumpdir = '/Volumes/Lenzinho/MyLinux/calc.vasp.5.3.3/development/am_lib/examples/custom/VN_pd';
load_dumps

clf; hold on;
clist = get(gca,'defaultAxesColorOrder');
plot(tbdr.E_input','-','color',clist(1,:));
plot(tbdr.E_optimized','-','color',clist(2,:));
axis tight; ylim = get(gca,'ylim'); set(gca,'ylim',ylim)
plot(tbdr.E_dft','-','color',clist(3,:));
% 
% width = 1;
% clist = get(gca,'defaultAxesColorOrder');
% set(h1,'color',clist(1,:),'linewidth',width)
% set(h2,'color',clist(2,:),'linewidth',width)
% set(h3,'color',clist(3,:),'linewidth',width)
% 
% legend([h1(1),h2(1),h3(1)],{'starting','optimizd','DFT'})