%% step 1: fresh start
clear;clc;
tiny=1.0E-6;
%% step 2: oad symmetry rep
dumpdir = '/Volumes/Lenzinho/MyLinux/calc.vasp.5.3.3/development/am_lib/examples/custom/VN_pd/out/tb/debug/tbpg';
load_dumps
%% step 3: set Fermi energy, load dispersion
EF = 9.50006808;
[dr,bz]=load_vasp_eigenval('../../../EIGENVAL');
% dr.E=dr.E-EF;
%% step 4: load primitive cell, convert k-points to cart 
[pc]  = load_poscar('../../uc/outfile.primitive');
bz.kpt = pc.recbas*bz.kpt;
%% step 5: set # of bands to skip
% skip the s band
band_skip = 1;
%% step 6: fit gamma
% choose kpt_mask to only select gamma
selector_kpoint = 1;
% choose shell_mask to only select zero-th neighbors
selector_shell = [1,2];
% set size of x to be the mininmum value that works
x = rand(3,1);
% perform fit
fun = @(x)tb_model_fitter_compute_residual(bz,dr,tbpg,x,selector_kpoint,selector_shell,band_skip);
opts= optimoptions('lsqnonlin','Display','iter','MaxIter',5);
[x,res] = lsqnonlin(fun,x,[],[],opts);
% plot results
for j = 1:bz.nkpts
    H = get_H_numeric_cart(tbpg, x, bz.kpt(:,j), selector_shell);
    tbdr.E(:,j) = sort(real(eig(H)));
end
figure(1); set(gcf,'color','white'); clf;
plot(1:bz.nkpts,dr.E,'-')
hold on;
plot(1:bz.nkpts,tbdr.E,'.')
hold off; axis tight; grid on;
ylim([-10 15]);
% save x
xmin = x;
%% step 7: fit high symmetry points to 1st nearest neighbors
% add more points to mask
selector_kpoint = round(linspace(1,bz.nkpts,10));
% set irreducible shell to include first neighbors
selector_shell = [1,2,3];
% simulated annealing temperature 
kT = 5;
% simulated annealing-type minimization
for i = 1:100
    if i == 1
        % initialize
        x = zeros(5,1);
        x(1:length(xmin)) = xmin;
        last_res = Inf;
    else
        % modify x values on second pass
        x = xmin + kT.*rand(size(xmin)).*abs(xmin)./max(abs(xmin));
    end
    % perform fit
    fun = @(x)tb_model_fitter_compute_residual(bz,dr,tbpg,x,selector_kpoint,selector_shell,band_skip);
    opts= optimoptions('lsqnonlin','Display','iter','MaxIter',5,'Display','off');
    [x,res] = lsqnonlin(fun,x,[],[],opts);
    %
    fprintf('%i %f %f\n',i,res,last_res)
    % save best parameter
    if res < last_res
        last_res = res;
        xmin = x;
        % plot result
        figure(1); set(gcf,'color','white'); clf;
        plot(1:bz.nkpts,dr.E,'-')
        hold on;
        plot(1:bz.nkpts,tbdr.E,'.')
        hold off; axis tight; box on; grid on; ylim([-10 15]); drawnow;
        for j = 1:bz.nkpts
            H = get_H_numeric_cart(tbpg, x, bz.kpt(:,j), selector_shell);
            tbdr.E(:,j) = sort(real(eig(H)));
        end
    end
end
%% step 8: repeat with second nearest neighbors
% add more points to mask
selector_kpoint = round(linspace(1,bz.nkpts,10));
% set irreducible shell to include first neighbors
selector_shell = [1,2,3,4,5];
% simulated annealing temperature 
kT = 5;
% simulated annealing-type minimization
for i = 1:100
    if i == 1
        % initialize
        x = zeros(14,1);
        x(1:length(xmin)) = xmin;
        last_res = Inf; % alread set last step
    else
        % modify x values on second pass
        x = xmin + kT.*(rand(size(xmin))-0.5).*abs(xmin)./max(abs(xmin));
    end
    % perform fit
    fun = @(x)tb_model_fitter_compute_residual(bz,dr,tbpg,x,selector_kpoint,selector_shell,band_skip);
    opts= optimoptions('lsqnonlin','Display','iter','MaxIter',5,'Display','off');
    [x,res] = lsqnonlin(fun,x,[],[],opts);
    %
    fprintf('%i %f %f\n',i,res,last_res)
    % save best parameter
    if res < last_res
        last_res = res;
        xmin = x;
        % plot result
        figure(1); set(gcf,'color','white'); clf;
        plot(1:bz.nkpts,dr.E,'-')
        hold on;
        plot(1:bz.nkpts,tbdr.E,'.')
        hold off; axis tight; box on; grid on; ylim([-10 15]); drawnow;
        for j = 1:bz.nkpts
            H = get_H_numeric_cart(tbpg, x, bz.kpt(:,j), selector_shell);
            tbdr.E(:,j) = sort(real(eig(H)));
        end
    end
end
%% generate plot
% band_skip = 1;
% % load xmin
load VN_pd_optimized_tb_parameters.mat
% % load dumps
% dumpdir = '/Volumes/Lenzinho/MyLinux/calc.vasp.5.3.3/development/am_lib/examples/custom/VN_pd/out/tb/debug/tbpg';
% load_dumps
% set shells
selector_shell = [1,2,3,4,5];
% set kpoint
kstart=[
    0.000000   0.000000   0.000000
    0.000000   0.500000   0.500000
    0.000000   0.000000   0.000000
    ]';
kend=[
    1.000000   0.500000   0.500000
    0.000000   0.000000   0.000000
    0.500000   0.500000   0.500000
    ]';
klabel={'G','X','G','L'};

dumpdir = '/Volumes/Lenzinho/MyLinux/calc.vasp.5.3.3/development/am_lib/examples/custom/VN_fit/out/tb/debug/tbpg';
load_dumps
[pc] = load_poscar('../../uc/outfile.primitive');
[bz] = get_kpoint_path(pc.bas,kstart,kend,40);

k=0;
for i = 1:bz.npaths
for j = 1:bz.ndivs
    k=k+1;
%     H(:,:,k) = get_H_numeric_frac(tbpg, xmin, bz.path(i).kpt_frac(:,j),selector_shell);
    H(:,:,k) = get_H_numeric_cart(tbpg, xmin, bz.path(i).kpt_cart(:,j),selector_shell);
    if (abs(norm(H(:,:,k)-H(:,:,k)'))>tiny) 
        fprintf('H is not hermitian\n')
        return
    end
    % bz.path(i).D(:,j) = sort(real(eig(H)));
end
end
% use eigenshuffle for nice band sepearation
[V,D]=eigenshuffle(H);
k=0;
for i = 1:bz.npaths
for j = 1:bz.ndivs
k=k+1;
bz.path(i).D(:,j) = real(D(:,k));
end
end


figure(1); set(gcf,'color','white'); clf; hold on;
X = [bz.path(:).x]';
Y_tb = [bz.path(:).D];
Y_dft= dr.E;
inds = round(linspace(1,length(bz.path)*bz.ndivs,30));
h_dft=plot(X,Y_dft([1:(end-band_skip)]+band_skip,:),'-');
% h_tb=plot(X,Y_tb,'-');
h_tb_marker=plot(X(inds),Y_tb(:,inds),'.');
colors = get(gca,'defaultAxesColorOrder');
markers= {'o','s','d','^','v','>','<','p','h'};
for i = 1:length(h_tb_marker)
%     set(h_tb(i),'color',colors(mod(i,size(colors,1))+1,:));
    set(h_dft(i),'color',[0.87 0.87 0.95],'linewidth',3);
    set(h_tb_marker(i),'color',colors(mod(i,size(colors,1))+1,:));
    set(h_tb_marker(i),'marker',markers{mod(i,length(markers))+1},'markersize',4);
    set(h_tb_marker(i),'MarkerFaceColor',colors(mod(i,size(colors,1))+1,:));
end
set(gca,'XMinorTick','on','YMinorTick','on')
hold off; axis tight; box on;
ylim([-10 6])

% get ticks
ticks(1) = 0;
for i = 1:bz.npaths
    ticks(i+1) = bz.path(i).x(end);
end
set(gca,'XTick',ticks);
set(gca,'Xticklabel',klabel);
ylabel('Energy [eV]')
xlabel('k')
%%
h=line([0,max([bz.path(:).x])],[0,0]);
set(h,'color',[0.1 0.1 0.1],'linewidth',2,'LineStyle',':');

set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'PaperSize',[10 10])
set(gcf,'PaperPosition',[0 0, 1.61803398875,1]*3)
print(gcf,'-depsc','band')
