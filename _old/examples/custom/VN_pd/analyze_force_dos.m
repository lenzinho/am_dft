%% plot DOS
clear;clc;
s = xml_read('save.dos');
%
x = str2num(s.E.value);
y = reshape(str2num(s.pD.value),s.pD.shape);
z = str2num(s.iD.value);
%
figure(1);clf;set(gcf,'color','white');hold on;
xp=x;
yp=              z; h(1)=plot(x,yp,'c.');
yp=cumsum(sum(y(  :,:),1))*(xp(2)-xp(1)); h(1)=plot(x,yp,'k');
yp=sum(y(6:8,:),1); h(2)=plot(x,yp,'b');
yp=sum(y(1:5,:),1); h(3)=plot(x,yp,'g');
%
xlabel('E [eV]');
ylabel('g [states/eV]');
axis tight; grid on; box on;
line([-10 10],[0 0],'color','k')
%
% ylim([-1 2]*6)
% xlim([-10 2])
legend(h,{'tot','p','d'})
%
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 7.62 7.62/1.61803398875]*1.5);
set(gcf,'PaperSize',[7.62 7.62/1.61803398875]*1.5);
print(gcf,'-dpdf','dos.pdf')
%% plot force-projected DOS
clear;clc;
s = xml_read('save.dos_force_7');
%
x = str2num(s.E.value);
y = reshape(str2num(s.pD.value),s.pD.shape);
%
figure(1);clf;set(gcf,'color','white');hold on;
xp=x; ytot =(sum(sum(y,1).*(xp<0))*(xp(2)-xp(1)));
yp=sum(y(  :,:),1)./ytot; h(1)=plot(x,yp,'k');
yp=sum(y(6:8,:),1)./ytot; h(2)=plot(x,yp,'b');
yp=sum(y(1:5,:),1)./ytot; h(3)=plot(x,yp,'g');
%
xlabel('E [eV]');
ylabel('force-projected dos [states/eV * eV/eta]');
axis tight; grid on; box on;
line([-10 10],[0 0],'color','k')
%
title(ytot)
ylim([-0.2 0.4])
xlim([-10 2])
legend(h,{'tot','p','d'})
%
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 7.62 7.62/1.61803398875]*1.5);
set(gcf,'PaperSize',[7.62 7.62/1.61803398875]*1.5);
print(gcf,'-dpdf','dos_force_def_7.pdf')

%% plot dispersion
clear;clc
s = xml_read('save.bzp');
d = xml_read('save.dr');
%
t = str2num(s.tick.value); % ticks
l = strsplit(s.kpt_symb.value,' '); % labels
x = str2num(s.x.value); % x
y = reshape(str2num(d.E.value),d.E.shape); % energies
%
figure(1);clf;set(gcf,'color','white');hold on;
plot(x,y,'k')
set(gca,'XTick',t)
set(gca,'XTickLabel',l)
axis tight; grid on; box on;

%% plot force-projected dos
s = xml_read('save.bzp_tbforce');
d = xml_read('save.dr_tbforce');
%
t = str2num(s.tick.value); % ticks
l = strsplit(s.kpt_symb.value,' '); % labels
x = str2num(s.x.value); % x
y = reshape(str2num(d.E.value),d.E.shape); % energies
c = reshape(str2num(d.weight.value),d.weight.shape); % color
%
figure(1);clf;set(gcf,'color','white');hold on;
for i = 1:size(y,1)
    xp=x(:)'; yp=y(i,:); cp=c(i,:); zp=zeros(size(x));
    surface([xp;xp],[yp;yp],[zp;zp],[cp;cp],'facecol','no','edgecol','interp','linew',1);
end
%
% clist = haxby_colormap(100);
clist = lightbertlein_colormap(100,'redblue');
colormap(clist(:,:).^(4/2))
cbh = colorbar('southoutside');
caxis([-3 3]);
cbh.Position(4)=cbh.Position(4)/3;
cbh.Position(2)=cbh.Position(2)/4;
%
ylabel('E [eV]');
xlabel('k');
set(gca,'XTick',t);
set(gca,'XTickLabel',l);
axis tight; grid on; box on;
ylim([-10 6])
%
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 7.62 7.62/1.61803398875]*1.5);
set(gcf,'PaperSize',[7.62 7.62/1.61803398875]*1.5);
print(gcf,'-dpdf','tbforce_def_0.pdf')
