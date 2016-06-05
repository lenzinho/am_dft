clear;clc;
tiny=1.0E-6;

[pg] = load_tb_point_group();
[pc] = load_poscar('outfile.primitive');

k = sym('k',[1,3],'real');
v = sym('v',[1,11],'real');
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
% Si Tight Binding model (sps*):
%
%   t  = 1.35750 % [cart]
%
%   a  = pi*(+ k2*t - k1*t + k3*t)*2i
%   b  = pi*(+ k1*t - k2*t + k3*t)*2i
%   c  = pi*(+ k1*t + k2*t - k3*t)*2i
%   d  = pi*(- k1*t - k2*t - k3*t)*2i
%   
%   g1 = + exp(-a) + exp(-b) - exp(-c) - exp(-d)
%   g2 = + exp(-a) - exp(-b) + exp(-c) - exp(-d)
%   g3 = + exp(-a) - exp(-b) - exp(-c) + exp(-d)
%   g4 = + exp(-a) + exp(-b) + exp(-c) + exp(-d)
%   
%   g5 = + exp(+a) + exp(+b) - exp(+c) - exp(+d) % c.c. of g1
%   g6 = + exp(+a) - exp(+b) + exp(+c) - exp(+d) % c.c. of g2
%   g7 = + exp(+a) - exp(+b) - exp(+c) + exp(+d) % c.c. of g3
%   g8 = + exp(+a) + exp(+b) + exp(+c) + exp(+d) % c.c. of g4
%   
%   [    v1,    v2,       0,       0,       0,  g4*v5,  g4*v6,  g3*v8,  -g2*v8,  -g1*v8]
%   [    v2,    v3,       0,       0,       0,  g4*v6,  g4*v7,  g3*v9,  -g2*v9,  -g1*v9]
%   [     0,     0,      v4,       0,       0, -g3*v8, -g3*v9, g4*v11,  g1*v10,  g2*v10]
%   [     0,     0,       0,      v4,       0,  g2*v8,  g2*v9, g1*v10,  g4*v11, -g3*v10]
%   [     0,     0,       0,       0,      v4,  g1*v8,  g1*v9, g2*v10, -g3*v10,  g4*v11]
%   [ g8*v5, g8*v6,  -g7*v8,  -g6*v8,  -g5*v8,     v1,     v2,      0,       0,       0]
%   [ g8*v6, g8*v7,  -g7*v9,  -g6*v9,  -g5*v9,     v2,     v3,      0,       0,       0]
%   [ g7*v8, g7*v9,  g8*v11, -g5*v10, -g6*v10,      0,      0,     v4,       0,       0]
%   [ g6*v8, g6*v9, -g5*v10,  g8*v11, -g7*v10,      0,      0,      0,      v4,       0]
%   [ g5*v8, g5*v9, -g6*v10, -g7*v10,  g8*v11,      0,      0,      0,       0,      v4]
%
%   REORDERED MATRIX AS IN VOGL: 
%
%       s_a     s_c     px_a     py_a     pz_a    px_c     py_c     pz_c    s_a     s_c
%     -----   -----    -----    -----    -----   -----    -----    -----  -----   ----- 
%   [    v1,  g4*v5,       0,       0,       0,  g3*v8,  -g2*v8,  -g1*v8,    v2,  g4*v6]  s_a
%   [ g8*v5,     v1,  -g7*v8,  -g6*v8,  -g5*v8,      0,       0,       0, g8*v6,     v2]  s_c
%   [     0, -g3*v8,      v4,       0,       0, g4*v11,  g1*v10,  g2*v10,     0, -g3*v9]  px_a
%   [     0,  g2*v8,       0,      v4,       0, g1*v10,  g4*v11, -g3*v10,     0,  g2*v9]  py_a
%   [     0,  g1*v8,       0,       0,      v4, g2*v10, -g3*v10,  g4*v11,     0,  g1*v9]  pz_a
%   [ g7*v8,      0,  g8*v11, -g5*v10, -g6*v10,     v4,       0,       0, g7*v9,      0]  px_c
%   [ g6*v8,      0, -g5*v10,  g8*v11, -g7*v10,      0,      v4,       0, g6*v9,      0]  py_c
%   [ g5*v8,      0, -g6*v10, -g7*v10,  g8*v11,      0,       0,      v4, g5*v9,      0]  pz_c
%   [    v2,  g4*v6,       0,       0,       0,  g3*v9,  -g2*v9,  -g1*v9,    v3,  g4*v7]  s_a
%   [ g8*v6,     v2,  -g7*v9,  -g6*v9,  -g5*v9,      0,       0,       0, g8*v7,     v3]  s_c
%
%   OBTAINED USING THE TRANSFORMATION
%
%   T(1:10,1:10) = 0
%   T( 1, 1) = 1; T( 6, 8) = 1;
%   T( 2, 6) = 1; T( 7, 9) = 1;
%   T( 3, 3) = 1; T( 8,10) = 1;
%   T( 4, 4) = 1; T( 9, 2) = 1;
%   T( 5, 5) = 1; T(10, 7) = 1;
%
%   H_VOGL = T*H*T'
%



































