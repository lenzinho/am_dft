function [H] = get_H_symbolic_frac(pg,v,kpt,shell)
i2pi = 2*sqrt(-1)*pi;
H(1:8,1:8) = 0;
H = sym(H);
if any(shell==1)
H(1:5,1:5) = H(1:5,1:5) + pg.sym(1:5,1:5,1)' * getV1(v(1:2)) * pg.sym(1:5,1:5,1) * exp(i2pi*dot(kpt,[+0.00000 +0.00000 +0.00000]));
end
if any(shell==2)
H(6:8,6:8) = H(6:8,6:8) + pg.sym(6:8,6:8,1)' * getV2(v(3:3)) * pg.sym(6:8,6:8,1) * exp(i2pi*dot(kpt,[+0.00000 +0.00000 +0.00000]));
end
if any(shell==3)
H(1:5,6:8) = H(1:5,6:8) + pg.sym(1:5,1:5,1)' * getV3(v(4:5)) * pg.sym(6:8,6:8,1) * exp(i2pi*dot(kpt,[+0.50000 -0.50000 -0.50000]));
H(1:5,6:8) = H(1:5,6:8) + pg.sym(1:5,1:5,5)' * getV3(v(4:5)) * pg.sym(6:8,6:8,5) * exp(i2pi*dot(kpt,[-0.50000 +0.50000 -0.50000]));
H(1:5,6:8) = H(1:5,6:8) + pg.sym(1:5,1:5,7)' * getV3(v(4:5)) * pg.sym(6:8,6:8,7) * exp(i2pi*dot(kpt,[-0.50000 -0.50000 +0.50000]));
H(1:5,6:8) = H(1:5,6:8) + pg.sym(1:5,1:5,9)' * getV3(v(4:5)) * pg.sym(6:8,6:8,9) * exp(i2pi*dot(kpt,[+0.50000 +0.50000 -0.50000]));
H(1:5,6:8) = H(1:5,6:8) + pg.sym(1:5,1:5,8)' * getV3(v(4:5)) * pg.sym(6:8,6:8,8) * exp(i2pi*dot(kpt,[+0.50000 -0.50000 +0.50000]));
H(1:5,6:8) = H(1:5,6:8) + pg.sym(1:5,1:5,2)' * getV3(v(4:5)) * pg.sym(6:8,6:8,2) * exp(i2pi*dot(kpt,[-0.50000 +0.50000 +0.50000]));
H(6:8,1:5) = H(6:8,1:5) + pg.sym(6:8,6:8,1)' * getV3(v(4:5))' * pg.sym(1:5,1:5,1) * exp(i2pi*dot(kpt,[-0.50000 +0.50000 +0.50000]));
H(6:8,1:5) = H(6:8,1:5) + pg.sym(6:8,6:8,5)' * getV3(v(4:5))' * pg.sym(1:5,1:5,5) * exp(i2pi*dot(kpt,[+0.50000 -0.50000 +0.50000]));
H(6:8,1:5) = H(6:8,1:5) + pg.sym(6:8,6:8,7)' * getV3(v(4:5))' * pg.sym(1:5,1:5,7) * exp(i2pi*dot(kpt,[+0.50000 +0.50000 -0.50000]));
H(6:8,1:5) = H(6:8,1:5) + pg.sym(6:8,6:8,9)' * getV3(v(4:5))' * pg.sym(1:5,1:5,9) * exp(i2pi*dot(kpt,[-0.50000 -0.50000 +0.50000]));
H(6:8,1:5) = H(6:8,1:5) + pg.sym(6:8,6:8,8)' * getV3(v(4:5))' * pg.sym(1:5,1:5,8) * exp(i2pi*dot(kpt,[-0.50000 +0.50000 -0.50000]));
H(6:8,1:5) = H(6:8,1:5) + pg.sym(6:8,6:8,2)' * getV3(v(4:5))' * pg.sym(1:5,1:5,2) * exp(i2pi*dot(kpt,[+0.50000 -0.50000 -0.50000]));
end
if any(shell==4)
H(1:5,1:5) = H(1:5,1:5) + pg.sym(1:5,1:5,1)' * getV4(v(6:11)) * pg.sym(1:5,1:5,1) * exp(i2pi*dot(kpt,[+0.00000 +0.00000 +1.00000]));
H(1:5,1:5) = H(1:5,1:5) + pg.sym(1:5,1:5,14)' * getV4(v(6:11)) * pg.sym(1:5,1:5,14) * exp(i2pi*dot(kpt,[+0.00000 +1.00000 +0.00000]));
H(1:5,1:5) = H(1:5,1:5) + pg.sym(1:5,1:5,15)' * getV4(v(6:11)) * pg.sym(1:5,1:5,15) * exp(i2pi*dot(kpt,[-1.00000 +0.00000 +1.00000]));
H(1:5,1:5) = H(1:5,1:5) + pg.sym(1:5,1:5,4)' * getV4(v(6:11)) * pg.sym(1:5,1:5,4) * exp(i2pi*dot(kpt,[-1.00000 +1.00000 +0.00000]));
H(1:5,1:5) = H(1:5,1:5) + pg.sym(1:5,1:5,16)' * getV4(v(6:11)) * pg.sym(1:5,1:5,16) * exp(i2pi*dot(kpt,[+1.00000 +0.00000 +0.00000]));
H(1:5,1:5) = H(1:5,1:5) + pg.sym(1:5,1:5,12)' * getV4(v(6:11)) * pg.sym(1:5,1:5,12) * exp(i2pi*dot(kpt,[+0.00000 -1.00000 +1.00000]));
H(1:5,1:5) = H(1:5,1:5) + pg.sym(1:5,1:5,7)' * getV4(v(6:11)) * pg.sym(1:5,1:5,7) * exp(i2pi*dot(kpt,[+0.00000 +1.00000 -1.00000]));
H(1:5,1:5) = H(1:5,1:5) + pg.sym(1:5,1:5,9)' * getV4(v(6:11)) * pg.sym(1:5,1:5,9) * exp(i2pi*dot(kpt,[-1.00000 +0.00000 +0.00000]));
H(1:5,1:5) = H(1:5,1:5) + pg.sym(1:5,1:5,2)' * getV4(v(6:11)) * pg.sym(1:5,1:5,2) * exp(i2pi*dot(kpt,[+1.00000 -1.00000 +0.00000]));
H(1:5,1:5) = H(1:5,1:5) + pg.sym(1:5,1:5,10)' * getV4(v(6:11)) * pg.sym(1:5,1:5,10) * exp(i2pi*dot(kpt,[+1.00000 +0.00000 -1.00000]));
H(1:5,1:5) = H(1:5,1:5) + pg.sym(1:5,1:5,6)' * getV4(v(6:11)) * pg.sym(1:5,1:5,6) * exp(i2pi*dot(kpt,[+0.00000 -1.00000 +0.00000]));
H(1:5,1:5) = H(1:5,1:5) + pg.sym(1:5,1:5,3)' * getV4(v(6:11)) * pg.sym(1:5,1:5,3) * exp(i2pi*dot(kpt,[+0.00000 +0.00000 -1.00000]));
end
if any(shell==5)
H(6:8,6:8) = H(6:8,6:8) + pg.sym(6:8,6:8,1)' * getV5(v(12:14)) * pg.sym(6:8,6:8,1) * exp(i2pi*dot(kpt,[+0.00000 +0.00000 +1.00000]));
H(6:8,6:8) = H(6:8,6:8) + pg.sym(6:8,6:8,14)' * getV5(v(12:14)) * pg.sym(6:8,6:8,14) * exp(i2pi*dot(kpt,[+0.00000 +1.00000 +0.00000]));
H(6:8,6:8) = H(6:8,6:8) + pg.sym(6:8,6:8,15)' * getV5(v(12:14)) * pg.sym(6:8,6:8,15) * exp(i2pi*dot(kpt,[-1.00000 +0.00000 +1.00000]));
H(6:8,6:8) = H(6:8,6:8) + pg.sym(6:8,6:8,4)' * getV5(v(12:14)) * pg.sym(6:8,6:8,4) * exp(i2pi*dot(kpt,[-1.00000 +1.00000 +0.00000]));
H(6:8,6:8) = H(6:8,6:8) + pg.sym(6:8,6:8,16)' * getV5(v(12:14)) * pg.sym(6:8,6:8,16) * exp(i2pi*dot(kpt,[+1.00000 +0.00000 +0.00000]));
H(6:8,6:8) = H(6:8,6:8) + pg.sym(6:8,6:8,12)' * getV5(v(12:14)) * pg.sym(6:8,6:8,12) * exp(i2pi*dot(kpt,[+0.00000 -1.00000 +1.00000]));
H(6:8,6:8) = H(6:8,6:8) + pg.sym(6:8,6:8,7)' * getV5(v(12:14)) * pg.sym(6:8,6:8,7) * exp(i2pi*dot(kpt,[+0.00000 +1.00000 -1.00000]));
H(6:8,6:8) = H(6:8,6:8) + pg.sym(6:8,6:8,9)' * getV5(v(12:14)) * pg.sym(6:8,6:8,9) * exp(i2pi*dot(kpt,[-1.00000 +0.00000 +0.00000]));
H(6:8,6:8) = H(6:8,6:8) + pg.sym(6:8,6:8,2)' * getV5(v(12:14)) * pg.sym(6:8,6:8,2) * exp(i2pi*dot(kpt,[+1.00000 -1.00000 +0.00000]));
H(6:8,6:8) = H(6:8,6:8) + pg.sym(6:8,6:8,10)' * getV5(v(12:14)) * pg.sym(6:8,6:8,10) * exp(i2pi*dot(kpt,[+1.00000 +0.00000 -1.00000]));
H(6:8,6:8) = H(6:8,6:8) + pg.sym(6:8,6:8,6)' * getV5(v(12:14)) * pg.sym(6:8,6:8,6) * exp(i2pi*dot(kpt,[+0.00000 -1.00000 +0.00000]));
H(6:8,6:8) = H(6:8,6:8) + pg.sym(6:8,6:8,3)' * getV5(v(12:14)) * pg.sym(6:8,6:8,3) * exp(i2pi*dot(kpt,[+0.00000 +0.00000 -1.00000]));
end
end
function [a] = getV1(v)
a(3,3) = v(1);
a(5,5) = v(2);
a(1,1) = +1.0000000000000000*a(3,3);
a(2,2) = +1.0000000000000000*a(5,5);
a(4,4) = +1.0000000000000000*a(5,5);
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
a(3,3) = v(1);
a(1,1) = +1.0000000000000000*a(3,3);
a(2,2) = +1.0000000000000000*a(3,3);
a(2,1) = 0;
a(3,1) = 0;
a(1,2) = 0;
a(3,2) = 0;
a(1,3) = 0;
a(2,3) = 0;
end
function [a] = getV3(v)
a(4,2) = v(1);
a(3,3) = v(2);
a(5,1) = +1.0000000000000000*a(4,2);
a(1,3) = +1.7320508075688770*a(3,3);
a(1,1) = 0;
a(2,1) = 0;
a(3,1) = 0;
a(4,1) = 0;
a(1,2) = 0;
a(2,2) = 0;
a(3,2) = 0;
a(5,2) = 0;
a(2,3) = 0;
a(4,3) = 0;
a(5,3) = 0;
end
function [a] = getV4(v)
a(1,3) = v(1);
a(3,3) = v(2);
a(3,4) = v(3);
a(4,4) = v(4);
a(2,5) = v(5);
a(5,5) = v(6);
a(1,1) = -1.1547005383792519*a(1,3)+1.0000000000000000*a(3,3);
a(3,1) = +1.0000000000000000*a(1,3);
a(4,1) = -1.7320508075688770*a(3,4);
a(2,2) = +1.0000000000000000*a(5,5);
a(5,2) = +1.0000000000000000*a(2,5);
a(4,3) = +1.0000000000000000*a(3,4);
a(1,4) = -1.7320508075688770*a(3,4);
a(2,1) = 0;
a(5,1) = 0;
a(1,2) = 0;
a(3,2) = 0;
a(4,2) = 0;
a(2,3) = 0;
a(5,3) = 0;
a(2,4) = 0;
a(5,4) = 0;
a(1,5) = 0;
a(3,5) = 0;
a(4,5) = 0;
end
function [a] = getV5(v)
a(1,1) = v(1);
a(2,3) = v(2);
a(3,3) = v(3);
a(2,2) = +1.0000000000000000*a(3,3);
a(3,2) = +1.0000000000000000*a(2,3);
a(2,1) = 0;
a(3,1) = 0;
a(1,2) = 0;
a(1,3) = 0;
end
