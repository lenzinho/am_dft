function [a] = getV2(v)
a(2,3) = v(1);
a(3,3) = v(2);
a(3,6) = v(3);
a(6,6) = v(4);
a(2,8) = v(5);
a(3,8) = v(6);
a(6,8) = v(7);
a(7,8) = v(8);
a(8,8) = v(9);
a(1,1) = +1.00000*a(3,3);
a(2,1) = -1.00000*a(2,3);
a(3,1) = +1.00000*a(2,3);
a(4,1) = -1.73205*a(3,6);
a(5,1) = -1.00000*a(3,8);
a(6,1) = +1.00000*a(3,6);
a(7,1) = -1.00000*a(2,8);
a(8,1) = +1.00000*a(3,8);
a(1,2) = -1.00000*a(2,3);
a(2,2) = +1.00000*a(3,3);
a(3,2) = +1.00000*a(2,3);
a(5,2) = -1.00000*a(3,8);
a(6,2) = -2.00000*a(3,6);
a(7,2) = +1.00000*a(3,8);
a(8,2) = -1.00000*a(2,8);
a(1,3) = +1.00000*a(2,3);
a(4,3) = -1.73205*a(3,6);
a(5,3) = -1.00000*a(2,8);
a(6,3) = -1.00000*a(3,6);
a(7,3) = -1.00000*a(3,8);
a(8,3) = -1.00000*a(3,8);
a(1,4) = +1.73205*a(3,6);
a(3,4) = +1.73205*a(3,6);
a(4,4) = +1.00000*a(6,6);
a(5,4) = +0.86603*a(6,8);
a(7,4) = +0.86603*a(6,8);
a(1,5) = +1.00000*a(3,8);
a(2,5) = +1.00000*a(3,8);
a(3,5) = +1.00000*a(2,8);
a(4,5) = +0.86603*a(6,8);
a(5,5) = +1.00000*a(8,8);
a(6,5) = +0.50000*a(6,8);
a(7,5) = -1.00000*a(7,8);
a(8,5) = -1.00000*a(7,8);
a(1,6) = -1.00000*a(3,6);
a(2,6) = +2.00000*a(3,6);
a(5,6) = +0.50000*a(6,8);
a(7,6) = -0.50000*a(6,8);
a(8,6) = +1.00000*a(6,8);
a(1,7) = +1.00000*a(2,8);
a(2,7) = -1.00000*a(3,8);
a(3,7) = +1.00000*a(3,8);
a(4,7) = +0.86603*a(6,8);
a(5,7) = -1.00000*a(7,8);
a(6,7) = -0.50000*a(6,8);
a(7,7) = +1.00000*a(8,8);
a(8,7) = +1.00000*a(7,8);
a(1,8) = -1.00000*a(3,8);
a(5,8) = -1.00000*a(7,8);
a(4,2) = 0;
a(2,4) = 0;
a(6,4) = 0;
a(8,4) = 0;
a(4,6) = 0;
a(4,8) = 0;
end
