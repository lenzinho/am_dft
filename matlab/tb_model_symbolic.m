clear;clc;
tiny=1.0E-6;

[pg] = load_tb_point_group();
[pc] = load_poscar('outfile.primitive');

k = sym('k',[1,3],'real');
v = sym('v',[1,10],'real');
H = getH_symb(pg, v, k);

syms a b c d
[H,a]=subexpr(H,a);
[H,b]=subexpr(H,b);
[H,c]=subexpr(H,c);
[H,d]=subexpr(H,d);
H=simplify(H,'steps',1000);
syms g1 g2 g3 g4 g5 g6 g7 g8
[H,g1]=subexpr(H,g1);
[H,g2]=subexpr(H,g2);
[H,g3]=subexpr(H,g3);
[H,g4]=subexpr(H,g4);
[H,g5]=subexpr(H,g5);
[H,g6]=subexpr(H,g6);
[H,g7]=subexpr(H,g7);
[H,g8]=subexpr(H,g8);
H

%
% Si Tight Binding model:
%
% t  = 1.35750 % [cart]
%
% a  = pi*(+ k2*t - k1*t + k3*t)*2i
% b  = pi*(+ k1*t - k2*t + k3*t)*2i
% c  = pi*(+ k1*t + k2*t - k3*t)*2i
% d  = pi*(- k1*t - k2*t - k3*t)*2i
%
% g1 = exp(-a) - exp(-b) + exp(-c) - exp(-d) % c.c. of g5
% g2 = exp(-a) + exp(-b) - exp(-c) - exp(-d) % c.c. of g6
% g3 = exp(-a) - exp(-b) - exp(-c) + exp(-d) % c.c. of g7
% g4 = exp(-a) + exp(-b) + exp(-c) + exp(-d) % c.c. of g8
%
% [    v1,      0,      0,      0,  g4*v3, g3*v4, -g1*v4, -g2*v4]
% [     0,     v2,      0,      0, -g3*v4, g4*v6,  g2*v5,  g1*v5]
% [     0,      0,     v2,      0,  g1*v4, g2*v5,  g4*v6, -g3*v5]
% [     0,      0,      0,     v2,  g2*v4, g1*v5, -g3*v5,  g4*v6]
% [ g8*v3, -g7*v4, -g5*v4, -g6*v4,     v1,     0,      0,      0]
% [ g7*v4,  g8*v6, -g6*v5, -g5*v5,      0,    v2,      0,      0]
% [ g5*v4, -g6*v5,  g8*v6, -g7*v5,      0,     0,     v2,      0]
% [ g6*v4, -g5*v5, -g7*v5,  g8*v6,      0,     0,      0,     v2]


