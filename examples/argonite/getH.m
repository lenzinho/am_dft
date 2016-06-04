function [H] = getH(pg,v,kpt)
i2pi = 2*sqrt(-1)*pi;
H(1:80,1:80) = 0;
% irreducible shell 1
H(1:4,1:4) = H(1:4,1:4) + pg.sym(1:4,1:4,1)' * getV1(v(1:4)) * pg.sym(1:4,1:4,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 1
H(5:8,5:8) = H(5:8,5:8) + pg.sym(5:8,5:8,1)' * getV1(v(1:4)) * pg.sym(5:8,5:8,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 1
H(9:12,9:12) = H(9:12,9:12) + pg.sym(9:12,9:12,1)' * getV1(v(1:4)) * pg.sym(9:12,9:12,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 1
H(13:16,13:16) = H(13:16,13:16) + pg.sym(13:16,13:16,1)' * getV1(v(1:4)) * pg.sym(13:16,13:16,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 2
H(17:20,17:20) = H(17:20,17:20) + pg.sym(17:20,17:20,1)' * getV2(v(4:7)) * pg.sym(17:20,17:20,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 2
H(21:24,21:24) = H(21:24,21:24) + pg.sym(21:24,21:24,1)' * getV2(v(4:7)) * pg.sym(21:24,21:24,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 2
H(25:28,25:28) = H(25:28,25:28) + pg.sym(25:28,25:28,1)' * getV2(v(4:7)) * pg.sym(25:28,25:28,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 2
H(29:32,29:32) = H(29:32,29:32) + pg.sym(29:32,29:32,1)' * getV2(v(4:7)) * pg.sym(29:32,29:32,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 3
H(17:20,33:36) = H(17:20,33:36) + pg.sym(17:20,17:20,1)' * getV3(v(7:16)) * pg.sym(33:36,33:36,1) * exp(i2pi*dot(kpt,[+0.00000, +1.28188, -0.04371]));
% irreducible shell 3
H(21:24,37:40) = H(21:24,37:40) + pg.sym(21:24,21:24,1)' * getV3(v(7:16)) * pg.sym(37:40,37:40,1) * exp(i2pi*dot(kpt,[+0.00000, -1.28188, +0.04371]));
% irreducible shell 3
H(25:28,41:44) = H(25:28,41:44) + pg.sym(25:28,25:28,1)' * getV3(v(7:16)) * pg.sym(41:44,41:44,1) * exp(i2pi*dot(kpt,[+0.00000, +1.28188, +0.04371]));
% irreducible shell 3
H(29:32,45:48) = H(29:32,45:48) + pg.sym(29:32,29:32,1)' * getV3(v(7:16)) * pg.sym(45:48,45:48,1) * exp(i2pi*dot(kpt,[+0.00000, -1.28188, -0.04371]));
% irreducible shell 3
H(33:36,17:20) = H(33:36,17:20) + pg.sym(33:36,33:36,1)' * getV3(v(7:16))' * pg.sym(17:20,17:20,1) * exp(i2pi*dot(kpt,[+0.00000, -1.28188, +0.04371]));
% irreducible shell 3
H(37:40,21:24) = H(37:40,21:24) + pg.sym(37:40,37:40,1)' * getV3(v(7:16))' * pg.sym(21:24,21:24,1) * exp(i2pi*dot(kpt,[+0.00000, +1.28188, -0.04371]));
% irreducible shell 3
H(41:44,25:28) = H(41:44,25:28) + pg.sym(41:44,41:44,1)' * getV3(v(7:16))' * pg.sym(25:28,25:28,1) * exp(i2pi*dot(kpt,[+0.00000, -1.28188, -0.04371]));
% irreducible shell 3
H(45:48,29:32) = H(45:48,29:32) + pg.sym(45:48,45:48,1)' * getV3(v(7:16))' * pg.sym(29:32,29:32,1) * exp(i2pi*dot(kpt,[+0.00000, +1.28188, +0.04371]));
% irreducible shell 4
H(17:20,49:52) = H(17:20,49:52) + pg.sym(17:20,17:20,1)' * getV4(v(16:31)) * pg.sym(49:52,49:52,1) * exp(i2pi*dot(kpt,[+1.11299, -0.65028, +0.01286]));
H(17:20,49:52) = H(17:20,49:52) + pg.sym(17:20,17:20,6)' * getV4(v(16:31)) * pg.sym(49:52,49:52,6) * exp(i2pi*dot(kpt,[-1.11299, -0.65028, +0.01286]));
% irreducible shell 4
H(21:24,65:68) = H(21:24,65:68) + pg.sym(21:24,21:24,1)' * getV4(v(16:31)) * pg.sym(65:68,65:68,1) * exp(i2pi*dot(kpt,[+1.11299, +0.65028, -0.01286]));
H(21:24,65:68) = H(21:24,65:68) + pg.sym(21:24,21:24,6)' * getV4(v(16:31)) * pg.sym(65:68,65:68,6) * exp(i2pi*dot(kpt,[-1.11299, +0.65028, -0.01286]));
% irreducible shell 4
H(25:28,77:80) = H(25:28,77:80) + pg.sym(25:28,25:28,1)' * getV4(v(16:31)) * pg.sym(77:80,77:80,1) * exp(i2pi*dot(kpt,[+1.11299, -0.65028, -0.01286]));
H(25:28,77:80) = H(25:28,77:80) + pg.sym(25:28,25:28,6)' * getV4(v(16:31)) * pg.sym(77:80,77:80,6) * exp(i2pi*dot(kpt,[-1.11299, -0.65028, -0.01286]));
% irreducible shell 4
H(29:32,61:64) = H(29:32,61:64) + pg.sym(29:32,29:32,1)' * getV4(v(16:31)) * pg.sym(61:64,61:64,1) * exp(i2pi*dot(kpt,[+1.11299, +0.65028, +0.01286]));
H(29:32,61:64) = H(29:32,61:64) + pg.sym(29:32,29:32,6)' * getV4(v(16:31)) * pg.sym(61:64,61:64,6) * exp(i2pi*dot(kpt,[-1.11299, +0.65028, +0.01286]));
% irreducible shell 5
H(33:36,33:36) = H(33:36,33:36) + pg.sym(33:36,33:36,1)' * getV5(v(31:34)) * pg.sym(33:36,33:36,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 5
H(37:40,37:40) = H(37:40,37:40) + pg.sym(37:40,37:40,1)' * getV5(v(31:34)) * pg.sym(37:40,37:40,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 5
H(41:44,41:44) = H(41:44,41:44) + pg.sym(41:44,41:44,1)' * getV5(v(31:34)) * pg.sym(41:44,41:44,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 5
H(45:48,45:48) = H(45:48,45:48) + pg.sym(45:48,45:48,1)' * getV5(v(31:34)) * pg.sym(45:48,45:48,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 5
H(49:52,49:52) = H(49:52,49:52) + pg.sym(49:52,49:52,1)' * getV5(v(31:34)) * pg.sym(49:52,49:52,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 5
H(53:56,53:56) = H(53:56,53:56) + pg.sym(53:56,53:56,1)' * getV5(v(31:34)) * pg.sym(53:56,53:56,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 5
H(57:60,57:60) = H(57:60,57:60) + pg.sym(57:60,57:60,1)' * getV5(v(31:34)) * pg.sym(57:60,57:60,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 5
H(61:64,61:64) = H(61:64,61:64) + pg.sym(61:64,61:64,1)' * getV5(v(31:34)) * pg.sym(61:64,61:64,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 5
H(65:68,65:68) = H(65:68,65:68) + pg.sym(65:68,65:68,1)' * getV5(v(31:34)) * pg.sym(65:68,65:68,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 5
H(69:72,69:72) = H(69:72,69:72) + pg.sym(69:72,69:72,1)' * getV5(v(31:34)) * pg.sym(69:72,69:72,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 5
H(73:76,73:76) = H(73:76,73:76) + pg.sym(73:76,73:76,1)' * getV5(v(31:34)) * pg.sym(73:76,73:76,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 5
H(77:80,77:80) = H(77:80,77:80) + pg.sym(77:80,77:80,1)' * getV5(v(31:34)) * pg.sym(77:80,77:80,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
% irreducible shell 6
H(33:36,49:52) = H(33:36,49:52) + pg.sym(33:36,33:36,1)' * getV6(v(34:49)) * pg.sym(49:52,49:52,1) * exp(i2pi*dot(kpt,[+1.11299, -1.93217, +0.05657]));
H(33:36,49:52) = H(33:36,49:52) + pg.sym(33:36,33:36,6)' * getV6(v(34:49)) * pg.sym(49:52,49:52,6) * exp(i2pi*dot(kpt,[-1.11299, -1.93217, +0.05657]));
% irreducible shell 6
H(37:40,65:68) = H(37:40,65:68) + pg.sym(37:40,37:40,1)' * getV6(v(34:49)) * pg.sym(65:68,65:68,1) * exp(i2pi*dot(kpt,[+1.11299, +1.93217, -0.05657]));
H(37:40,65:68) = H(37:40,65:68) + pg.sym(37:40,37:40,6)' * getV6(v(34:49)) * pg.sym(65:68,65:68,6) * exp(i2pi*dot(kpt,[-1.11299, +1.93217, -0.05657]));
% irreducible shell 6
H(41:44,77:80) = H(41:44,77:80) + pg.sym(41:44,41:44,1)' * getV6(v(34:49)) * pg.sym(77:80,77:80,1) * exp(i2pi*dot(kpt,[+1.11299, -1.93217, -0.05657]));
H(41:44,77:80) = H(41:44,77:80) + pg.sym(41:44,41:44,6)' * getV6(v(34:49)) * pg.sym(77:80,77:80,6) * exp(i2pi*dot(kpt,[-1.11299, -1.93217, -0.05657]));
% irreducible shell 6
H(45:48,61:64) = H(45:48,61:64) + pg.sym(45:48,45:48,1)' * getV6(v(34:49)) * pg.sym(61:64,61:64,1) * exp(i2pi*dot(kpt,[+1.11299, +1.93217, +0.05657]));
H(45:48,61:64) = H(45:48,61:64) + pg.sym(45:48,45:48,6)' * getV6(v(34:49)) * pg.sym(61:64,61:64,6) * exp(i2pi*dot(kpt,[-1.11299, +1.93217, +0.05657]));
% irreducible shell 7
H(49:52,17:20) = H(49:52,17:20) + pg.sym(49:52,49:52,1)' * getV7(v(49:64))' * pg.sym(17:20,17:20,1) * exp(i2pi*dot(kpt,[-1.11299, +0.65028, -0.01286]));
% irreducible shell 7
H(53:56,21:24) = H(53:56,21:24) + pg.sym(53:56,53:56,1)' * getV7(v(49:64))' * pg.sym(21:24,21:24,1) * exp(i2pi*dot(kpt,[+1.11299, -0.65028, +0.01286]));
% irreducible shell 7
H(57:60,25:28) = H(57:60,25:28) + pg.sym(57:60,57:60,1)' * getV7(v(49:64))' * pg.sym(25:28,25:28,1) * exp(i2pi*dot(kpt,[+1.11299, +0.65028, +0.01286]));
% irreducible shell 7
H(61:64,29:32) = H(61:64,29:32) + pg.sym(61:64,61:64,1)' * getV7(v(49:64))' * pg.sym(29:32,29:32,1) * exp(i2pi*dot(kpt,[-1.11299, -0.65028, -0.01286]));
% irreducible shell 7
H(65:68,21:24) = H(65:68,21:24) + pg.sym(65:68,65:68,1)' * getV7(v(49:64))' * pg.sym(21:24,21:24,1) * exp(i2pi*dot(kpt,[-1.11299, -0.65028, +0.01286]));
% irreducible shell 7
H(69:72,17:20) = H(69:72,17:20) + pg.sym(69:72,69:72,1)' * getV7(v(49:64))' * pg.sym(17:20,17:20,1) * exp(i2pi*dot(kpt,[+1.11299, +0.65028, -0.01286]));
% irreducible shell 7
H(73:76,29:32) = H(73:76,29:32) + pg.sym(73:76,73:76,1)' * getV7(v(49:64))' * pg.sym(29:32,29:32,1) * exp(i2pi*dot(kpt,[+1.11299, -0.65028, -0.01286]));
% irreducible shell 7
H(77:80,25:28) = H(77:80,25:28) + pg.sym(77:80,77:80,1)' * getV7(v(49:64))' * pg.sym(25:28,25:28,1) * exp(i2pi*dot(kpt,[-1.11299, +0.65028, +0.01286]));
% irreducible shell 8
H(49:52,69:72) = H(49:52,69:72) + pg.sym(49:52,49:52,1)' * getV8(v(64:68)) * pg.sym(69:72,69:72,1) * exp(i2pi*dot(kpt,[-2.22598, +0.00000, +0.00000]));
% irreducible shell 8
H(53:56,65:68) = H(53:56,65:68) + pg.sym(53:56,53:56,1)' * getV8(v(64:68)) * pg.sym(65:68,65:68,1) * exp(i2pi*dot(kpt,[+2.22598, +0.00000, +0.00000]));
% irreducible shell 8
H(57:60,77:80) = H(57:60,77:80) + pg.sym(57:60,57:60,1)' * getV8(v(64:68)) * pg.sym(77:80,77:80,1) * exp(i2pi*dot(kpt,[+2.22598, +0.00000, +0.00000]));
% irreducible shell 8
H(61:64,73:76) = H(61:64,73:76) + pg.sym(61:64,61:64,1)' * getV8(v(64:68)) * pg.sym(73:76,73:76,1) * exp(i2pi*dot(kpt,[-2.22598, +0.00000, +0.00000]));
% irreducible shell 8
H(65:68,53:56) = H(65:68,53:56) + pg.sym(65:68,65:68,1)' * getV8(v(64:68))' * pg.sym(53:56,53:56,1) * exp(i2pi*dot(kpt,[-2.22598, +0.00000, +0.00000]));
% irreducible shell 8
H(69:72,49:52) = H(69:72,49:52) + pg.sym(69:72,69:72,1)' * getV8(v(64:68))' * pg.sym(49:52,49:52,1) * exp(i2pi*dot(kpt,[+2.22598, +0.00000, +0.00000]));
% irreducible shell 8
H(73:76,61:64) = H(73:76,61:64) + pg.sym(73:76,73:76,1)' * getV8(v(64:68))' * pg.sym(61:64,61:64,1) * exp(i2pi*dot(kpt,[+2.22598, +0.00000, +0.00000]));
% irreducible shell 8
H(77:80,57:60) = H(77:80,57:60) + pg.sym(77:80,77:80,1)' * getV8(v(64:68))' * pg.sym(57:60,57:60,1) * exp(i2pi*dot(kpt,[-2.22598, +0.00000, +0.00000]));
% irreducible shell 9
H(49:52,33:36) = H(49:52,33:36) + pg.sym(49:52,49:52,1)' * getV9(v(68:83))' * pg.sym(33:36,33:36,1) * exp(i2pi*dot(kpt,[-1.11299, +1.93217, -0.05657]));
% irreducible shell 9
H(53:56,37:40) = H(53:56,37:40) + pg.sym(53:56,53:56,1)' * getV9(v(68:83))' * pg.sym(37:40,37:40,1) * exp(i2pi*dot(kpt,[+1.11299, -1.93217, +0.05657]));
% irreducible shell 9
H(57:60,41:44) = H(57:60,41:44) + pg.sym(57:60,57:60,1)' * getV9(v(68:83))' * pg.sym(41:44,41:44,1) * exp(i2pi*dot(kpt,[+1.11299, +1.93217, +0.05657]));
% irreducible shell 9
H(61:64,45:48) = H(61:64,45:48) + pg.sym(61:64,61:64,1)' * getV9(v(68:83))' * pg.sym(45:48,45:48,1) * exp(i2pi*dot(kpt,[-1.11299, -1.93217, -0.05657]));
% irreducible shell 9
H(65:68,37:40) = H(65:68,37:40) + pg.sym(65:68,65:68,1)' * getV9(v(68:83))' * pg.sym(37:40,37:40,1) * exp(i2pi*dot(kpt,[-1.11299, -1.93217, +0.05657]));
% irreducible shell 9
H(69:72,33:36) = H(69:72,33:36) + pg.sym(69:72,69:72,1)' * getV9(v(68:83))' * pg.sym(33:36,33:36,1) * exp(i2pi*dot(kpt,[+1.11299, +1.93217, -0.05657]));
% irreducible shell 9
H(73:76,45:48) = H(73:76,45:48) + pg.sym(73:76,73:76,1)' * getV9(v(68:83))' * pg.sym(45:48,45:48,1) * exp(i2pi*dot(kpt,[+1.11299, -1.93217, -0.05657]));
% irreducible shell 9
H(77:80,41:44) = H(77:80,41:44) + pg.sym(77:80,77:80,1)' * getV9(v(68:83))' * pg.sym(41:44,41:44,1) * exp(i2pi*dot(kpt,[-1.11299, +1.93217, +0.05657]));
end
