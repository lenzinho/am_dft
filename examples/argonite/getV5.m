function [a] = getV5(v)
a(1,1) = v(1);
a(2,2) = v(2);
a(3,3) = v(3);
a(4,4) = v(4);
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
