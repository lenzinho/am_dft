function spy_color(A)
    tiny = 1E-6;
    unique_values = unique(A);
    unique_values = unique_values(abs(unique_values)>=tiny);
    unique_values_len = length(unique_values);
    clist = parula(unique_values_len+1);
    figure(1); 
    set(gcf,'color','w'); clf; 
    plot(0,0);
    hold on;
    for i = 1:unique_values_len
        spy(abs(A(:,:)-unique_values(i))<tiny);
    end
    hold off;
    h=get(gca,'children');
    for i = 1:unique_values_len
        set(h(i),'color',clist(i,:))
    end
    axis on;
    xlim([0.5, size(A,2)+0.5]);
    ylim([0.5, size(A,1)+0.5]);
    set(get(gca,'Xlabel'),'String','');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'XMinorGrid','on');
    set(gca,'YMinorGrid','on');
    box on;
end