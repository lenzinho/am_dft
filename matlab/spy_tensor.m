function spy_tensor(A)
    clist = parula(size(A,3));
    figure(1); set(gcf,'color','w'); clf; hold on;
    for i = 1:size(A,3)
        spy(A(:,:,i),'r')
    end
    hold off;
    h=get(gca,'children');
    for i = 1:size(A,3)
        set(h(i),'color',clist(i,:))
    end
    axis off;
end