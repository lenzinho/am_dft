
clear;clc;
tiny=1.0E-6;

% kpoint path [frac]
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


v(1) = -4.2 ; % E(s,a)
v(2) = 0    ; % 0
v(3) = 6.685; % E(s*,a)
v(4) = 1.715; % E(p,a)
v(5) =-8.300; % V(s,s)
v(6) = 0 	; % 0
v(7) = 0 	; % V(s*,s*) ??? did not find it in table
v(8) = 5.729; % V(sa,pc)
v(9) = 5.375; % V(s*a,pc)
v(10)= 4.575; % V(x,y)
v(11)= 1.715; % V(x,x)

[pg] = load_tb_point_group();
[pc] = load_poscar('outfile.primitive');
[bz] = get_kpoint_path(pc.bas,kstart,kend,40);

for i = 1:bz.npaths
for j = 1:bz.ndivs
    % WHY IS THIS FACTOR OF TWO NEEDED HRE?!
%     H = get_H_model_frac(v,bz.path(i).kpt_frac(:,j)*2);
    H = get_H_numeric_cart(pg,v,bz.path(i).kpt_cart(:,j));
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
