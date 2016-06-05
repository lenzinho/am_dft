clear;clc;
tiny=1.0E-6;

%% ------------------------------------------------------------------------
% set options (works best with frac)
maxvars = 1000;
flags = 'frac';












%% ------------------------------------------------------------------------
fprintf('Building get_H_model:\n');
%% ------------------------------------------------------------------------
fprintf(' ... loading symmetries \n');
[pg] = load_tb_point_group();
%% ------------------------------------------------------------------------
fprintf(' ... loading poscar \n');
[pc] = load_poscar('outfile.primitive');
%% ------------------------------------------------------------------------
fprintf(' ... substituting kpoints \n');
k = sym('k%d',[1,3],'real');
%% ------------------------------------------------------------------------
fprintf(' ... getting symbolic hamiltonian\n');
for nVs = 1:maxvars
    v = sym('v%d',[1,nVs],'real');
    try 
        if     ~isempty(strfind(flags,'frac'))
            H = getH_symb_frac(pg, v, k);
        elseif ~isempty(strfind(flags,'cart'))
            H = getH_symb(pg, v, k);
        else
            fprintf('ERROR: frac/cart?\n')
            return
        end
        break
    end
end
%% ------------------------------------------------------------------------
fprintf(' ... substituting exponents in hamiltonian \n');
syms z real
a = sym('a%d',[1,maxvars],'real');
amax = 0;
for i = 1:maxvars
    [~,z]=subexpr(H,z);
    if isempty(findstr(char(z),'exp'))
        amax = amax+1;
        [H,a(amax)]=subexpr(H,sprintf('a(%i)',amax));
        a(amax) = simplify(a(amax),'steps',50);
    else
        break
    end
end
%% ------------------------------------------------------------------------
fprintf(' ... factoring v in hamiltonian \n');
H=simplify(H,'steps',1000);
%% ------------------------------------------------------------------------
fprintf(' ... substitute exponentials in hamiltonian \n');
g = sym('g%d',[1,maxvars],'real');
gmax = 0;
for i = 1:maxvars
    [~,z]=subexpr(H,z);
    if ~isempty(findstr(char(z),'exp'))
        gmax = gmax+1;
        [H,g(gmax)]=subexpr(H,sprintf('g(%i)',gmax));
    else
        break
    end
end
%% ------------------------------------------------------------------------
fprintf(' ... converting symbolic array k to accept numerical vectors\n');
for i = 1:amax
for j = 1:3
a(i) = subs(a(i),k(j),sprintf('k(%i)',j));
end
end
%% ------------------------------------------------------------------------
fprintf(' ... converting symbolic array v to accept numerical vectors\n');
for i = 1:amax
for j = 1:nVs
a(i) = subs(a(i),v(j),sprintf('v(%i)',j));
end
end
%% ------------------------------------------------------------------------
fprintf(' ... converting symbolic hamiltonian v to accept numerical vectors\n');
for i = 1:size(H,1)
for j = 1:size(H,2)
for m = 1:nVs
H(i,j) = subs(H(i,j),v(m),sprintf('v(%i)',m));
end
end
end
%% ------------------------------------------------------------------------
if     ~isempty(strfind(flags,'frac'))
    fprintf(' ... writting getH_model_frac.m\n');
    fid = fopen('getH_model_frac.m','w');
    fprintf(fid,'function [H] = getH_model_frac(v,k)\n');
elseif ~isempty(strfind(flags,'cart'))
    fprintf(' ... writting getH_model_cart.m\n');
    fid = fopen('getH_model_cart.m','w');
    fprintf(fid,'function [H] = getH_model_cart(v,k)\n');
else
    fprintf('ERROR: frac/cart?\n')
    return
end
fprintf(fid,'\n');
% print a
for i = 1:amax
    fprintf(fid,'a(%i) = ',i);
    fprintf(fid,'%s;',char(a(i)));
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
% print g
for i = 1:gmax
    fprintf(fid,'g(%i) = ',i);
    fprintf(fid,'%s;',char(g(i)));
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
% get widest field in H
fmax = 0;
for i = 1:size(H,1)
for j = 1:size(H,2)
if (length(char(H(i,j)))>fmax) 
    fmax = length(char(H(i,j)));
end
end
end
fmt = sprintf('%%%is',fmax);
% print H
fprintf(fid,'H = [');
for i = 1:size(H,1)
    if i~=1
        fprintf(fid,'     [');
    else
        fprintf(fid,'[');
    end
    for j = 1:size(H,2)
        fprintf(fid,fmt,char(H(i,j)));
        if (j==size(H,2)) 
            fprintf(fid,']');
        else
            fprintf(fid,',');
        end
    end
    if (i==size(H,1))
        fprintf(fid,'];\n');
    else
        fprintf(fid,';\n');
    end
end
fprintf(fid,'\n');
fprintf(fid,'end');

