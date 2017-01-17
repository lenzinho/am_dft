clear;clc;addpath('../essentials');

% load data
read_ = @(fname) fscanf(fopen(fname),'%f');
th_1d = read_('th.itx'); 
% phi_1d = read_('phi.itx');
N = 450;
th_1d  = [0:N]/N*2.26;
phi_1d = [0:N]/N*360;

% load rotation matrices in tight binding representation
tb  = xml_read('../save.tb');
sym = reshape(str2num(tb.group.sym.value),tb.group.sym.shape);

% set tight binding matrix elements
v = [1.158,-0.476,-4.085,1.045,1.038,-0.130,-0.220,-0.074,-0.533,0.055,0.112,-0.091,0.432,0.317];

% set units
hbar = 1; m = 1; eV = 1;
units = 3.622627*sqrt(eV*m)/hbar; % [1/nm] = sqrt(E [eV] * m_e )/hbar

% set parameters and constants
E_photon = 137.85;
E_binding = 142;
for E_inner = 30:0.5:50

    % resample
    rs_ = @(x,n) x(1:n:end);

    % compute k points
    [phi,th,E] = meshgrid(rs_(th_1d,15),rs_(phi_1d,3),E_photon);
    x = units*sqrt(2*m*E).*sind(phi).*cosd(th);
    y = units*sqrt(2*m*E).*sind(phi).*sind(th);
    z = units*sqrt(2*m*(E.*cosd(phi).^2+E_inner));

    % % plot kpoints
    % surf(x,y,z,'edgecolor','none'); view([0 0 1]);

    % change orientation to 011
    R_ = @(th) [1 0 0; 0 cosd(th) -sind(th); 0 sind(th) cosd(th)];
    
    % centering matrix 
    bas = 0.4134 * R_(45) * (ones(3)-eye(3))/2;

    % get kpoints
    nks = numel(x); k(1:nks,1:3) = [x(:),y(:),z(:)]*bas;

    % define gaussian function
    degauss = 1; gauss_ = @(x) exp(-abs(x).^2)./sqrt(pi);

    % get eigenvalues and spectral function
    for ik = 1:nks
        H = get_H_numeric_frac(sym, v, k(ik,1:3), [1:100]);
        A(:,ik) = sum( gauss_(eig(H)/degauss)/degauss, 1);
    end

    h=figure(1); set(h,'color','w'); clf;
    surf(x,y,reshape(A,size(x,1),size(x,2)),'edgecolor','none')
    view([0 0 1]); daspect([1 1 1]); grid off; box on; drawnow;
    saveas(h,sprintf('%2.1f.jpeg',E_inner));
end

%%
function [H] = get_H_numeric_frac(sym,v,kpt,shell)
    i2pi = 2*sqrt(-1)*pi;
    H(1:8,1:8) = 0;
    if any(shell==1)
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,1)' * getV1(v(1:2)) * sym(1:5,1:5,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
    end
    if any(shell==2)
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,1)' * getV2(v(3:3)) * sym(6:8,6:8,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
    end
    if any(shell==3)
    H(1:5,6:8) = H(1:5,6:8) + sym(1:5,1:5,1)' * getV3(v(4:5)) * sym(6:8,6:8,1) * exp(i2pi*dot(kpt,[-0.50000, -0.50000, +0.50000]));
    H(1:5,6:8) = H(1:5,6:8) + sym(1:5,1:5,10)' * getV3(v(4:5)) * sym(6:8,6:8,10) * exp(i2pi*dot(kpt,[-0.50000, +0.50000, -0.50000]));
    H(1:5,6:8) = H(1:5,6:8) + sym(1:5,1:5,7)' * getV3(v(4:5)) * sym(6:8,6:8,7) * exp(i2pi*dot(kpt,[+0.50000, -0.50000, -0.50000]));
    H(1:5,6:8) = H(1:5,6:8) + sym(1:5,1:5,9)' * getV3(v(4:5)) * sym(6:8,6:8,9) * exp(i2pi*dot(kpt,[-0.50000, +0.50000, +0.50000]));
    H(1:5,6:8) = H(1:5,6:8) + sym(1:5,1:5,6)' * getV3(v(4:5)) * sym(6:8,6:8,6) * exp(i2pi*dot(kpt,[+0.50000, -0.50000, +0.50000]));
    H(1:5,6:8) = H(1:5,6:8) + sym(1:5,1:5,2)' * getV3(v(4:5)) * sym(6:8,6:8,2) * exp(i2pi*dot(kpt,[+0.50000, +0.50000, -0.50000]));
    H(6:8,1:5) = H(6:8,1:5) + sym(6:8,6:8,1)' * getV3(v(4:5))' * sym(1:5,1:5,1) * exp(i2pi*dot(kpt,[+0.50000, +0.50000, -0.50000]));
    H(6:8,1:5) = H(6:8,1:5) + sym(6:8,6:8,10)' * getV3(v(4:5))' * sym(1:5,1:5,10) * exp(i2pi*dot(kpt,[+0.50000, -0.50000, +0.50000]));
    H(6:8,1:5) = H(6:8,1:5) + sym(6:8,6:8,7)' * getV3(v(4:5))' * sym(1:5,1:5,7) * exp(i2pi*dot(kpt,[-0.50000, +0.50000, +0.50000]));
    H(6:8,1:5) = H(6:8,1:5) + sym(6:8,6:8,9)' * getV3(v(4:5))' * sym(1:5,1:5,9) * exp(i2pi*dot(kpt,[+0.50000, -0.50000, -0.50000]));
    H(6:8,1:5) = H(6:8,1:5) + sym(6:8,6:8,6)' * getV3(v(4:5))' * sym(1:5,1:5,6) * exp(i2pi*dot(kpt,[-0.50000, +0.50000, -0.50000]));
    H(6:8,1:5) = H(6:8,1:5) + sym(6:8,6:8,2)' * getV3(v(4:5))' * sym(1:5,1:5,2) * exp(i2pi*dot(kpt,[-0.50000, -0.50000, +0.50000]));
    end
    if any(shell==4)
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,1)' * getV4(v(6:11)) * sym(1:5,1:5,1) * exp(i2pi*dot(kpt,[+1.00000, +0.00000, +0.00000]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,16)' * getV4(v(6:11)) * sym(1:5,1:5,16) * exp(i2pi*dot(kpt,[+0.00000, +1.00000, +0.00000]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,17)' * getV4(v(6:11)) * sym(1:5,1:5,17) * exp(i2pi*dot(kpt,[+1.00000, +0.00000, -1.00000]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,3)' * getV4(v(6:11)) * sym(1:5,1:5,3) * exp(i2pi*dot(kpt,[+0.00000, +1.00000, -1.00000]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,14)' * getV4(v(6:11)) * sym(1:5,1:5,14) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +1.00000]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,18)' * getV4(v(6:11)) * sym(1:5,1:5,18) * exp(i2pi*dot(kpt,[+1.00000, -1.00000, +0.00000]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,7)' * getV4(v(6:11)) * sym(1:5,1:5,7) * exp(i2pi*dot(kpt,[-1.00000, +1.00000, +0.00000]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,9)' * getV4(v(6:11)) * sym(1:5,1:5,9) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, -1.00000]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,2)' * getV4(v(6:11)) * sym(1:5,1:5,2) * exp(i2pi*dot(kpt,[+0.00000, -1.00000, +1.00000]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,5)' * getV4(v(6:11)) * sym(1:5,1:5,5) * exp(i2pi*dot(kpt,[-1.00000, +0.00000, +1.00000]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,8)' * getV4(v(6:11)) * sym(1:5,1:5,8) * exp(i2pi*dot(kpt,[+0.00000, -1.00000, +0.00000]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,4)' * getV4(v(6:11)) * sym(1:5,1:5,4) * exp(i2pi*dot(kpt,[-1.00000, +0.00000, +0.00000]));
    end
    if any(shell==5)
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,1)' * getV5(v(12:14)) * sym(6:8,6:8,1) * exp(i2pi*dot(kpt,[+1.00000, +0.00000, +0.00000]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,16)' * getV5(v(12:14)) * sym(6:8,6:8,16) * exp(i2pi*dot(kpt,[+0.00000, +1.00000, +0.00000]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,17)' * getV5(v(12:14)) * sym(6:8,6:8,17) * exp(i2pi*dot(kpt,[+1.00000, +0.00000, -1.00000]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,3)' * getV5(v(12:14)) * sym(6:8,6:8,3) * exp(i2pi*dot(kpt,[+0.00000, +1.00000, -1.00000]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,14)' * getV5(v(12:14)) * sym(6:8,6:8,14) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +1.00000]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,18)' * getV5(v(12:14)) * sym(6:8,6:8,18) * exp(i2pi*dot(kpt,[+1.00000, -1.00000, +0.00000]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,7)' * getV5(v(12:14)) * sym(6:8,6:8,7) * exp(i2pi*dot(kpt,[-1.00000, +1.00000, +0.00000]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,9)' * getV5(v(12:14)) * sym(6:8,6:8,9) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, -1.00000]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,2)' * getV5(v(12:14)) * sym(6:8,6:8,2) * exp(i2pi*dot(kpt,[+0.00000, -1.00000, +1.00000]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,5)' * getV5(v(12:14)) * sym(6:8,6:8,5) * exp(i2pi*dot(kpt,[-1.00000, +0.00000, +1.00000]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,8)' * getV5(v(12:14)) * sym(6:8,6:8,8) * exp(i2pi*dot(kpt,[+0.00000, -1.00000, +0.00000]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,4)' * getV5(v(12:14)) * sym(6:8,6:8,4) * exp(i2pi*dot(kpt,[-1.00000, +0.00000, +0.00000]));
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
end
function [H] = get_H_numeric_cart(sym,v,kpt,shell)
    i2pi = 2*sqrt(-1)*pi;
    H(1:8,1:8) = 0;
    if any(shell==1)
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,1)' * getV1(v(1:2)) * sym(1:5,1:5,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
    end
    if any(shell==2)
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,1)' * getV2(v(3:3)) * sym(6:8,6:8,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +0.00000]));
    end
    if any(shell==3)
    H(1:5,6:8) = H(1:5,6:8) + sym(1:5,1:5,1)' * getV3(v(4:5)) * sym(6:8,6:8,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, -2.03558]));
    H(1:5,6:8) = H(1:5,6:8) + sym(1:5,1:5,10)' * getV3(v(4:5)) * sym(6:8,6:8,10) * exp(i2pi*dot(kpt,[+0.00000, -2.03558, +0.00000]));
    H(1:5,6:8) = H(1:5,6:8) + sym(1:5,1:5,7)' * getV3(v(4:5)) * sym(6:8,6:8,7) * exp(i2pi*dot(kpt,[-2.03558, +0.00000, +0.00000]));
    H(1:5,6:8) = H(1:5,6:8) + sym(1:5,1:5,9)' * getV3(v(4:5)) * sym(6:8,6:8,9) * exp(i2pi*dot(kpt,[+2.03558, +0.00000, +0.00000]));
    H(1:5,6:8) = H(1:5,6:8) + sym(1:5,1:5,6)' * getV3(v(4:5)) * sym(6:8,6:8,6) * exp(i2pi*dot(kpt,[+0.00000, +2.03558, +0.00000]));
    H(1:5,6:8) = H(1:5,6:8) + sym(1:5,1:5,2)' * getV3(v(4:5)) * sym(6:8,6:8,2) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +2.03558]));
    H(6:8,1:5) = H(6:8,1:5) + sym(6:8,6:8,1)' * getV3(v(4:5))' * sym(1:5,1:5,1) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, +2.03558]));
    H(6:8,1:5) = H(6:8,1:5) + sym(6:8,6:8,10)' * getV3(v(4:5))' * sym(1:5,1:5,10) * exp(i2pi*dot(kpt,[+0.00000, +2.03558, +0.00000]));
    H(6:8,1:5) = H(6:8,1:5) + sym(6:8,6:8,7)' * getV3(v(4:5))' * sym(1:5,1:5,7) * exp(i2pi*dot(kpt,[+2.03558, +0.00000, +0.00000]));
    H(6:8,1:5) = H(6:8,1:5) + sym(6:8,6:8,9)' * getV3(v(4:5))' * sym(1:5,1:5,9) * exp(i2pi*dot(kpt,[-2.03558, +0.00000, +0.00000]));
    H(6:8,1:5) = H(6:8,1:5) + sym(6:8,6:8,6)' * getV3(v(4:5))' * sym(1:5,1:5,6) * exp(i2pi*dot(kpt,[+0.00000, -2.03558, +0.00000]));
    H(6:8,1:5) = H(6:8,1:5) + sym(6:8,6:8,2)' * getV3(v(4:5))' * sym(1:5,1:5,2) * exp(i2pi*dot(kpt,[+0.00000, +0.00000, -2.03558]));
    end
    if any(shell==4)
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,1)' * getV4(v(6:11)) * sym(1:5,1:5,1) * exp(i2pi*dot(kpt,[+0.00000, +2.03558, +2.03558]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,16)' * getV4(v(6:11)) * sym(1:5,1:5,16) * exp(i2pi*dot(kpt,[+2.03558, +0.00000, +2.03558]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,17)' * getV4(v(6:11)) * sym(1:5,1:5,17) * exp(i2pi*dot(kpt,[-2.03558, +0.00000, +2.03558]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,3)' * getV4(v(6:11)) * sym(1:5,1:5,3) * exp(i2pi*dot(kpt,[+0.00000, -2.03558, +2.03558]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,14)' * getV4(v(6:11)) * sym(1:5,1:5,14) * exp(i2pi*dot(kpt,[+2.03558, +2.03558, +0.00000]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,18)' * getV4(v(6:11)) * sym(1:5,1:5,18) * exp(i2pi*dot(kpt,[-2.03558, +2.03558, +0.00000]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,7)' * getV4(v(6:11)) * sym(1:5,1:5,7) * exp(i2pi*dot(kpt,[+2.03558, -2.03558, +0.00000]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,9)' * getV4(v(6:11)) * sym(1:5,1:5,9) * exp(i2pi*dot(kpt,[-2.03558, -2.03558, +0.00000]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,2)' * getV4(v(6:11)) * sym(1:5,1:5,2) * exp(i2pi*dot(kpt,[+0.00000, +2.03558, -2.03558]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,5)' * getV4(v(6:11)) * sym(1:5,1:5,5) * exp(i2pi*dot(kpt,[+2.03558, +0.00000, -2.03558]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,8)' * getV4(v(6:11)) * sym(1:5,1:5,8) * exp(i2pi*dot(kpt,[-2.03558, +0.00000, -2.03558]));
    H(1:5,1:5) = H(1:5,1:5) + sym(1:5,1:5,4)' * getV4(v(6:11)) * sym(1:5,1:5,4) * exp(i2pi*dot(kpt,[+0.00000, -2.03558, -2.03558]));
    end
    if any(shell==5)
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,1)' * getV5(v(12:14)) * sym(6:8,6:8,1) * exp(i2pi*dot(kpt,[+0.00000, +2.03558, +2.03558]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,16)' * getV5(v(12:14)) * sym(6:8,6:8,16) * exp(i2pi*dot(kpt,[+2.03558, +0.00000, +2.03558]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,17)' * getV5(v(12:14)) * sym(6:8,6:8,17) * exp(i2pi*dot(kpt,[-2.03558, +0.00000, +2.03558]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,3)' * getV5(v(12:14)) * sym(6:8,6:8,3) * exp(i2pi*dot(kpt,[+0.00000, -2.03558, +2.03558]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,14)' * getV5(v(12:14)) * sym(6:8,6:8,14) * exp(i2pi*dot(kpt,[+2.03558, +2.03558, +0.00000]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,18)' * getV5(v(12:14)) * sym(6:8,6:8,18) * exp(i2pi*dot(kpt,[-2.03558, +2.03558, +0.00000]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,7)' * getV5(v(12:14)) * sym(6:8,6:8,7) * exp(i2pi*dot(kpt,[+2.03558, -2.03558, +0.00000]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,9)' * getV5(v(12:14)) * sym(6:8,6:8,9) * exp(i2pi*dot(kpt,[-2.03558, -2.03558, +0.00000]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,2)' * getV5(v(12:14)) * sym(6:8,6:8,2) * exp(i2pi*dot(kpt,[+0.00000, +2.03558, -2.03558]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,5)' * getV5(v(12:14)) * sym(6:8,6:8,5) * exp(i2pi*dot(kpt,[+2.03558, +0.00000, -2.03558]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,8)' * getV5(v(12:14)) * sym(6:8,6:8,8) * exp(i2pi*dot(kpt,[-2.03558, +0.00000, -2.03558]));
    H(6:8,6:8) = H(6:8,6:8) + sym(6:8,6:8,4)' * getV5(v(12:14)) * sym(6:8,6:8,4) * exp(i2pi*dot(kpt,[+0.00000, -2.03558, -2.03558]));
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
end

