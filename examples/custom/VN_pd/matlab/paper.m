clear;clc;

tb=xml_read('../save.tb');
pg.sym=reshape(str2num(tb.group.sym.value),tb.group.sym.shape);
for c = [sqrt(3),2/sqrt(3),sqrt(3)/2]
pg.sym(abs(pg.sym-c)<1E-7) = c;
end

shell= [3];
v    = sym('v',[1,14],'real');
kpt  = sym('k',[1,3],'real');

H = get_H_symbolic_cart(pg,v,kpt,shell)
%%
rewrite(H,'sin')
%%
H=simplify(H,'steps',500)
% NOTE THAT IN THE PAPER d-1 and d0 states are switched to get E_g and T_2g together
% [                           0,                  0,                   0,                  0,                  0, 3^(1/2)*v5*sin(pi*a*k1)*2i,                   0, -3^(1/2)*v5*sin(pi*a*k3)*2i]
% [                           0,                  0,                   0,                  0,                  0,        -v4*sin(pi*a*k2)*2i, -v4*sin(pi*a*k1)*2i,                           0]
% [                           0,                  0,                   0,                  0,                  0,        -v5*sin(pi*a*k1)*2i,  v5*sin(pi*a*k2)*4i,         -v5*sin(pi*a*k3)*2i]
% [                           0,                  0,                   0,                  0,                  0,                          0, -v4*sin(pi*a*k3)*2i,         -v4*sin(pi*a*k2)*2i]
% [                           0,                  0,                   0,                  0,                  0,        -v4*sin(pi*a*k3)*2i,                   0,         -v4*sin(pi*a*k1)*2i]
% [ -3^(1/2)*v5*sin(pi*a*k1)*2i, v4*sin(pi*a*k2)*2i,  v5*sin(pi*a*k1)*2i,                  0, v4*sin(pi*a*k3)*2i,                          0,                   0,                           0]
% [                           0, v4*sin(pi*a*k1)*2i, -v5*sin(pi*a*k2)*4i, v4*sin(pi*a*k3)*2i,                  0,                          0,                   0,                           0]
% [  3^(1/2)*v5*sin(pi*a*k3)*2i,                  0,  v5*sin(pi*a*k3)*2i, v4*sin(pi*a*k2)*2i, v4*sin(pi*a*k1)*2i,                          0,                   0,                           0]
%  