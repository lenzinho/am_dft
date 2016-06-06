function [H] = get_H_symbolic_frac(pg,v,kpt)
i2pi = 2*sqrt(-1)*pi;
H(1:10,1:10) = 0;
H = sym(H);
% irreducible shell 1
H(1:5,1:5) = H(1:5,1:5) + pg.sym(1:5,1:5,1)' * getV1(v(1:4)) * pg.sym(1:5,1:5,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 1
H(6:10,6:10) = H(6:10,6:10) + pg.sym(6:10,6:10,1)' * getV1(v(1:4)) * pg.sym(6:10,6:10,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 2
H(1:5,6:10) = H(1:5,6:10) + pg.sym(1:5,1:5,1)' * getV2(v(5:11)) * pg.sym(6:10,6:10,1) * exp(i2pi*dot(kpt,[+0.12500, -0.37500, +0.12500]));
H(1:5,6:10) = H(1:5,6:10) + pg.sym(1:5,1:5,3)' * getV2(v(5:11)) * pg.sym(6:10,6:10,3) * exp(i2pi*dot(kpt,[+0.12500, +0.12500, -0.37500]));
H(1:5,6:10) = H(1:5,6:10) + pg.sym(1:5,1:5,4)' * getV2(v(5:11)) * pg.sym(6:10,6:10,4) * exp(i2pi*dot(kpt,[-0.37500, +0.12500, +0.12500]));
H(1:5,6:10) = H(1:5,6:10) + pg.sym(1:5,1:5,2)' * getV2(v(5:11)) * pg.sym(6:10,6:10,2) * exp(i2pi*dot(kpt,[+0.12500, +0.12500, +0.12500]));
% irreducible shell 2
H(6:10,1:5) = H(6:10,1:5) + pg.sym(6:10,6:10,1)' * getV2(v(5:11))' * pg.sym(1:5,1:5,1) * exp(i2pi*dot(kpt,[-0.12500, +0.37500, -0.12500]));
H(6:10,1:5) = H(6:10,1:5) + pg.sym(6:10,6:10,3)' * getV2(v(5:11))' * pg.sym(1:5,1:5,3) * exp(i2pi*dot(kpt,[-0.12500, -0.12500, +0.37500]));
H(6:10,1:5) = H(6:10,1:5) + pg.sym(6:10,6:10,4)' * getV2(v(5:11))' * pg.sym(1:5,1:5,4) * exp(i2pi*dot(kpt,[+0.37500, -0.12500, -0.12500]));
H(6:10,1:5) = H(6:10,1:5) + pg.sym(6:10,6:10,2)' * getV2(v(5:11))' * pg.sym(1:5,1:5,2) * exp(i2pi*dot(kpt,[-0.12500, -0.12500, -0.12500]));
end
function [a] = getV1(v)
a(1,1) = v(1);
a(1,2) = v(2);
a(2,2) = v(3);
a(5,5) = v(4);
a(2,1) = +1.00000*a(1,2);
a(3,3) = +1.00000*a(5,5);
a(4,4) = +1.00000*a(5,5);
a(3,1) = 0;
a(4,1) = 0;
a(5,1) = 0;
a(3,2) = 0;
a(4,2) = 0;
a(5,2) = 0;
a(1,3) = 0;
a(2,3) = 0;
a(4,3) = 0;
a(5,3) = 0;
a(1,4) = 0;
a(2,4) = 0;
a(3,4) = 0;
a(5,4) = 0;
a(1,5) = 0;
a(2,5) = 0;
a(3,5) = 0;
a(4,5) = 0;
end
function [a] = getV2(v)
a(1,1) = v(1);
a(1,2) = v(2);
a(2,2) = v(3);
a(1,5) = v(4);
a(2,5) = v(5);
a(4,5) = v(6);
a(5,5) = v(7);
a(2,1) = +1.00000*a(1,2);
a(3,1) = +1.00000*a(1,5);
a(4,1) = -1.00000*a(1,5);
a(5,1) = -1.00000*a(1,5);
a(3,2) = +1.00000*a(2,5);
a(4,2) = -1.00000*a(2,5);
a(5,2) = -1.00000*a(2,5);
a(1,3) = -1.00000*a(1,5);
a(2,3) = -1.00000*a(2,5);
a(3,3) = +1.00000*a(5,5);
a(4,3) = -1.00000*a(4,5);
a(5,3) = -1.00000*a(4,5);
a(1,4) = +1.00000*a(1,5);
a(2,4) = +1.00000*a(2,5);
a(3,4) = -1.00000*a(4,5);
a(4,4) = +1.00000*a(5,5);
a(5,4) = +1.00000*a(4,5);
a(3,5) = -1.00000*a(4,5);
end
