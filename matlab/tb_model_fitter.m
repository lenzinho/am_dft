
clear;clc;
tiny=1.0E-6;

[dr,bz]=load_vasp_eigenval();
[pg]  = load_tb_point_group();
[pc]  = load_poscar('outfile.primitive');

bz.kpt = pc.recbas*bz.kpt;
%% nonlinear least squares
clc;
% X = zeros(18,1);
mask = false(bz.nkpts);
mask([1,20,50]) = true;
X = zeros(18,1);
fun = @(x)tb_model_fitter_compute_residual(bz,dr,pg,x,mask);
opts= optimoptions('lsqnonlin','Display','iter');
X = lsqnonlin(fun,X,[],[],opts);

%% genetic algorithm
nvars = 18;
fun  = @(x) sum(abs(tb_model_fitter_compute_residual(bz,dr,pg,x,mask)));
pfun = @(options,state,flag) tb_model_plot_dr(options,state,flag,bz,dr,pg);
opts= gaoptimset('display','iter','PopulationSize',100,'Generations',100,'PlotFcn',pfun);
ft_ = ga(fun,nvars,[],[],[],[],[],[],[],opts)

%% simulated annealing
nvars = 18;
fun  = @(x) sum(abs(tb_model_fitter_compute_residual(bz,dr,pg,x)));
opts= saoptimset('display','iter','InitialTemperature',100000000,'ReannealInterval',10);
ft_ = simulannealbnd(fun,zeros(nvars,1),[],[],opts)
%% particle swarm
nvars=18;
fun  = @(x) sum(abs(tb_model_fitter_compute_residual(bz,dr,pg,x)));
opts= optimoptions(@particleswarm,'display','iter');
ft_ = particleswarm(fun,nvars,-ones(nvars,1)*max([dr.E]),ones(nvars,1)*max([dr.E]),opts)


%%

for j = 1:bz.nkpts
    H = get_H_numeric_cart(pg, ft_, bz.kpt(:,j));
    tbdr.E(:,j) = sort(real(eig(H)));
end

plot(1:bz.nkpts,dr.E)
hold on;
plot(1:bz.nkpts,tbdr.E,':')
hold off;

%%
fun = @(x) norm(tb_model_fitter_compute_residual(bz,dr,pg,x));