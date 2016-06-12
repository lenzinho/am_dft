function spy_color(A)
    tiny = 1E-6;
    unique_values = unique(A);
    unique_values = unique_values(abs(unique_values)>=tiny);
    unique_values_len = length(unique_values);
    clist = parula(unique_values_len+1);
    figure(1); set(gcf,'color','w'); clf; hold on;
    for i = 1:unique_values_len
        spy(abs(A(:,:)-unique_values(i))<tiny);
    end
    hold off;
    h=get(gca,'children');
    for i = 1:unique_values_len
        set(h(i),'color',clist(i,:))
    end
    axis off;
end