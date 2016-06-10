
clear;clc;
tiny=1.0E-6;

[dr,bz]=load_vasp_eigenval();
[pg]  = load_tb_point_group();
[pc]  = load_poscar('POSCAR');

bz.kpt = pc.recbas*bz.kpt;
mask = false(bz.nkpts,1);
mask([1,30,74,120]) = true;


%% nonlinear least squares
clc;
% X = zeros(18,1);
mask = false(bz.nkpts,1);
mask([1,20,40,80,100]) = true;
x = rand(18,1);
x(1:5) = 1:5;
fun = @(x)tb_model_fitter_compute_residual(bz,dr,pg,x,mask);
opts= optimoptions('lsqnonlin','Display','iter','MaxIter',100);
x = lsqnonlin(fun,x,[],[],opts);
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
%% manual optimizaiton as described in Matthises
clear;clc;
tiny=1.0E-6;

[dr,bz]=load_vasp_eigenval();
[pg]  = load_tb_point_group();
[pc]  = load_poscar('POSCAR');

bz.kpt = pc.recbas*bz.kpt;
mask = false(bz.nkpts,1);
mask([1,80:5:120]) = true;

x =rand(18,1);
for i = 1:100
    [d,F1,A,B] = tb_model_fitter_compute_delta_vector(bz,dr,pg,x,mask);
    x = x - 10*d;
    fprintf('%5i %10.5f \n',i,norm(F1));
	% plot
    for j = 1:bz.nkpts
        H = get_H_numeric_cart(pg, x, bz.kpt(:,j));
        tbdr.E(:,j) = sort(real(eig(H)));
    end
    plot(1:bz.nkpts,dr.E,'-',1:bz.nkpts,tbdr.E,'--');
    drawnow;
end


%% best x values so far:
x =[ ...
   10.7054    9.1940
    8.8019    9.1628
   -6.0901   -5.8567
    5.5098    6.4349
   -0.5899    0.5587
    0.9720    1.1763
   -0.5930    0.8716
   -0.0210    0.4194
   -0.1461    0.2234
    0.5594    0.0936
    0.1699    0.3596
    0.3560    0.8212
   -0.1889   -0.3945
   -0.2038   -0.2159
   -0.0875    0.8687
   -0.0016    0.5045
   -0.7587   -0.0810
    0.2809   -0.2493];

%%

for j = 1:bz.nkpts
    H = get_H_numeric_cart(pg, x, bz.kpt(:,j));
    tbdr.E(:,j) = sort(real(eig(H)));
end

plot(1:bz.nkpts,dr.E)
hold on;
plot(1:bz.nkpts,tbdr.E,'--')
hold off;

%%

kstart=[
    0.500000   0.500000   0.500000
    0.000000   0.000000   0.000000
    1.000000   0.500000   0.500000
    ]';
kend=[
    0.000000   0.000000   0.000000
    0.000000   0.500000   0.500000
    0.000000   0.000000   0.000000
    ]';
klabel={'L','G','X','G'};

[pg] = load_tb_point_group();
[pc] = load_poscar('outfile.primitive');
[bz] = get_kpoint_path(pc.bas,kstart,kend,40);

for i = 1:bz.npaths
for j = 1:bz.ndivs
    H = get_H_model_frac(x, bz.path(i).kpt_frac(:,j)*2);
    if (abs(norm(H-H'))>tiny) 
        fprintf('H is not hermitian\n')
        return
    end
    bz.path(i).D(:,j) = sort(real(eig(H)));
end
end

for i = 1:bz.npaths
    plot(bz.path(i).x,bz.path(i).D,'-')
    hold on;
end
hold off; axis tight; box on;

% get ticks
ticks(1) = 0;
for i = 1:bz.npaths
    ticks(i+1) = bz.path(i).x(end);
end
set(gca,'XTick',ticks);
set(gca,'Xticklabel',klabel);