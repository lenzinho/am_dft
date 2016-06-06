function [H] = get_H_symbolic_cart(pg,v,kpt)
i2pi = 2*sqrt(-1)*pi;
H(1:9,1:9) = 0;
H = sym(H);
% irreducible shell 1
H(1:5,1:5) = H(1:5,1:5) + pg.sym(1:5,1:5,1)' * getV1(v(1:2)) * pg.sym(1:5,1:5,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 2
H(1:5,6:9) = H(1:5,6:9) + pg.sym(1:5,1:5,1)' * getV2(v(3:5)) * pg.sym(6:9,6:9,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, -2.03524]));
H(1:5,6:9) = H(1:5,6:9) + pg.sym(1:5,1:5,5)' * getV2(v(3:5)) * pg.sym(6:9,6:9,5) * exp(i2pi*dot(kpt,[+0.00000, -2.03524, +0.00000]));
H(1:5,6:9) = H(1:5,6:9) + pg.sym(1:5,1:5,7)' * getV2(v(3:5)) * pg.sym(6:9,6:9,7) * exp(i2pi*dot(kpt,[-2.03524, +0.00000, +0.00000]));
H(1:5,6:9) = H(1:5,6:9) + pg.sym(1:5,1:5,9)' * getV2(v(3:5)) * pg.sym(6:9,6:9,9) * exp(i2pi*dot(kpt,[+2.03524, +0.00000, +0.00000]));
H(1:5,6:9) = H(1:5,6:9) + pg.sym(1:5,1:5,8)' * getV2(v(3:5)) * pg.sym(6:9,6:9,8) * exp(i2pi*dot(kpt,[+0.00000, +2.03524, +0.00000]));
H(1:5,6:9) = H(1:5,6:9) + pg.sym(1:5,1:5,2)' * getV2(v(3:5)) * pg.sym(6:9,6:9,2) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +2.03524]));
% irreducible shell 2
H(6:9,1:5) = H(6:9,1:5) + pg.sym(6:9,6:9,1)' * getV2(v(3:5))' * pg.sym(1:5,1:5,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +2.03524]));
H(6:9,1:5) = H(6:9,1:5) + pg.sym(6:9,6:9,5)' * getV2(v(3:5))' * pg.sym(1:5,1:5,5) * exp(i2pi*dot(kpt,[+0.00000, +2.03524, +0.00000]));
H(6:9,1:5) = H(6:9,1:5) + pg.sym(6:9,6:9,7)' * getV2(v(3:5))' * pg.sym(1:5,1:5,7) * exp(i2pi*dot(kpt,[+2.03524, +0.00000, +0.00000]));
H(6:9,1:5) = H(6:9,1:5) + pg.sym(6:9,6:9,9)' * getV2(v(3:5))' * pg.sym(1:5,1:5,9) * exp(i2pi*dot(kpt,[-2.03524, +0.00000, +0.00000]));
H(6:9,1:5) = H(6:9,1:5) + pg.sym(6:9,6:9,8)' * getV2(v(3:5))' * pg.sym(1:5,1:5,8) * exp(i2pi*dot(kpt,[+0.00000, -2.03524, +0.00000]));
H(6:9,1:5) = H(6:9,1:5) + pg.sym(6:9,6:9,2)' * getV2(v(3:5))' * pg.sym(1:5,1:5,2) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, -2.03524]));
% irreducible shell 3
H(6:9,6:9) = H(6:9,6:9) + pg.sym(6:9,6:9,1)' * getV3(v(6:7)) * pg.sym(6:9,6:9,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
end
function [a] = getV1(v)
a(3,3) = v(1);
a(5,5) = v(2);
a(1,1) = +1.00000*a(3,3);
a(2,2) = +1.00000*a(5,5);
a(4,4) = +1.00000*a(5,5);
a(2,1) = 0;
a(3,1) = 0;
a(4,1) = 0;
a(5,1) = 0;
a(1,2) = 0;
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
a(3,1) = v(1);
a(4,3) = v(2);
a(3,4) = v(3);
a(1,1) = +1.73205*a(3,1);
a(5,2) = +1.00000*a(4,3);
a(1,4) = +1.73205*a(3,4);
a(2,1) = 0;
a(4,1) = 0;
a(5,1) = 0;
a(1,2) = 0;
a(2,2) = 0;
a(3,2) = 0;
a(4,2) = 0;
a(1,3) = 0;
a(2,3) = 0;
a(3,3) = 0;
a(5,3) = 0;
a(2,4) = 0;
a(4,4) = 0;
a(5,4) = 0;
end
function [a] = getV3(v)
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
