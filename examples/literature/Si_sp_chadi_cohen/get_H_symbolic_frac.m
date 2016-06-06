function [H] = get_H_symbolic_frac(pg,v,kpt)
i2pi = 2*sqrt(-1)*pi;
H(1:8,1:8) = 0;
H = sym(H);
% irreducible shell 1
H(1:4,1:4) = H(1:4,1:4) + pg.sym(1:4,1:4,1)' * getV1(v(1:2)) * pg.sym(1:4,1:4,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 1
H(5:8,5:8) = H(5:8,5:8) + pg.sym(5:8,5:8,1)' * getV1(v(1:2)) * pg.sym(5:8,5:8,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 2
H(1:4,5:8) = H(1:4,5:8) + pg.sym(1:4,1:4,1)' * getV2(v(3:6)) * pg.sym(5:8,5:8,1) * exp(i2pi*dot(kpt,[+0.12500, -0.37500, +0.12500]));
H(1:4,5:8) = H(1:4,5:8) + pg.sym(1:4,1:4,3)' * getV2(v(3:6)) * pg.sym(5:8,5:8,3) * exp(i2pi*dot(kpt,[+0.12500, +0.12500, -0.37500]));
H(1:4,5:8) = H(1:4,5:8) + pg.sym(1:4,1:4,4)' * getV2(v(3:6)) * pg.sym(5:8,5:8,4) * exp(i2pi*dot(kpt,[-0.37500, +0.12500, +0.12500]));
H(1:4,5:8) = H(1:4,5:8) + pg.sym(1:4,1:4,2)' * getV2(v(3:6)) * pg.sym(5:8,5:8,2) * exp(i2pi*dot(kpt,[+0.12500, +0.12500, +0.12500]));
% irreducible shell 2
H(5:8,1:4) = H(5:8,1:4) + pg.sym(5:8,5:8,1)' * getV2(v(3:6))' * pg.sym(1:4,1:4,1) * exp(i2pi*dot(kpt,[-0.12500, +0.37500, -0.12500]));
H(5:8,1:4) = H(5:8,1:4) + pg.sym(5:8,5:8,3)' * getV2(v(3:6))' * pg.sym(1:4,1:4,3) * exp(i2pi*dot(kpt,[-0.12500, -0.12500, +0.37500]));
H(5:8,1:4) = H(5:8,1:4) + pg.sym(5:8,5:8,4)' * getV2(v(3:6))' * pg.sym(1:4,1:4,4) * exp(i2pi*dot(kpt,[+0.37500, -0.12500, -0.12500]));
H(5:8,1:4) = H(5:8,1:4) + pg.sym(5:8,5:8,2)' * getV2(v(3:6))' * pg.sym(1:4,1:4,2) * exp(i2pi*dot(kpt,[-0.12500, -0.12500, -0.12500]));
end
function [a] = getV1(v)
a(1,1) = v(1);
a(4,4) = v(2);
a(2,2) = +1.00000*a(4,4);
a(3,3) = +1.00000*a(4,4);
a(2,1) = 0;
a(3,1) = 0;
a(4,1) = 0;
a(1,2) = 0;
a(3,2) = 0;
a(4,2) = 0;
a(1,3) = 0;
a(2,3) = 0;
a(4,3) = 0;
a(1,4) = 0;
a(2,4) = 0;
a(3,4) = 0;
end
function [a] = getV2(v)
a(1,1) = v(1);
a(1,4) = v(2);
a(3,4) = v(3);
a(4,4) = v(4);
a(2,1) = +1.00000*a(1,4);
a(3,1) = -1.00000*a(1,4);
a(4,1) = -1.00000*a(1,4);
a(1,2) = -1.00000*a(1,4);
a(2,2) = +1.00000*a(4,4);
a(3,2) = -1.00000*a(3,4);
a(4,2) = -1.00000*a(3,4);
a(1,3) = +1.00000*a(1,4);
a(2,3) = -1.00000*a(3,4);
a(3,3) = +1.00000*a(4,4);
a(4,3) = +1.00000*a(3,4);
a(2,4) = -1.00000*a(3,4);
end
