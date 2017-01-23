clear;clc;

% define figure properties
fig = @(h)   set(h,'color','white');
axs = @(h,qt,ql) set(h,'Box','on','XTick',qt,'Xticklabel',ql,'ylim',[-9 4.5]);

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

% set path divisions and start and end vectors for each path segment
N = 50; nqs = 3;
qs=[0.500000   0.500000   0.500000
    0.000000   0.000000   0.000000
    1.000000   0.500000   0.500000]';
qe=[0.000000   0.000000   0.000000
    0.000000   0.500000   0.500000
    0.000000   0.000000   0.000000]';
ql={'L','G','X','G'};

% load rotation matrices in tight binding representation
tb  = xml_read('../save.tb');
sym = reshape(str2num(tb.group.sym.value),tb.group.sym.shape);

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

% get eigenvalues
for ik = 1:nks
    H = get_H_numeric_frac(sym, v, k(1:3,ik), [1:100]);
    E(:,ik) = sort(real(eig(H)));
end

% plot band structure 
figure(1); fig(gcf);
plot(x,E,'-');
axs(gca,qt,ql); axis tight; 
ylabel('Energy E');
xlabel('Wavevector k');
