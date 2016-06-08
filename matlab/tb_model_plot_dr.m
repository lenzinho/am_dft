function state = tb_model_plot_dr(options,state,flag,bz,dr,pg)
[~,i] = min(state.Score);
x = state.Population(i,:);

for j = 1:bz.nkpts
    H = get_H_numeric_cart(pg, x, bz.kpt(:,j));
    tbdr.E(:,j) = sort(real(eig(H)));
end
figure(1)
plot(1:bz.nkpts,dr.E);       hold on;
plot(1:bz.nkpts,tbdr.E,':'); hold off;
drawnow;
end