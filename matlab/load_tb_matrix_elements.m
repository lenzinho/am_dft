function [pp] = load_tb_matrix_elements()
    fid=fopen('outfile.tb_matrix_elements');
    % number of atoms
    pp.nshells = sscanf(fgetl(fid),'%i');
    fprintf('shells %i\n',pp.nshells)
    for i = 1:pp.nshells
        % divider
        fgetl(fid);
        k = sscanf(fgetl(fid),'%i');
        pp.shell(k).m = sscanf(fgetl(fid),'%i');
        pp.shell(k).n = sscanf(fgetl(fid),'%i');
        SE(1:2) = sscanf(fgetl(fid),'%i'); % [m]
        pp.S(pp.shell(k).m) = SE(1); 
        pp.E(pp.shell(k).m) = SE(2); 
        SE(1:2) = sscanf(fgetl(fid),'%i'); % [n]
        pp.S(pp.shell(k).n) = SE(1); 
        pp.E(pp.shell(k).n) = SE(2); 
        pp.shell(k).npairs = sscanf(fgetl(fid),'%i');
        for p = 1:pp.shell(k).npairs
        pp.shell(k).pair(p).tau(1:3) = sscanf(fgetl(fid),'%f'); % [cart]
        fprintf(' ... shell %4i pair %3i: %10.5f %10.5f %10.5f\n',k,p,pp.shell(k).pair(p).tau)
        for j = 1:(pp.E(pp.shell(k).m)-pp.S(pp.shell(k).m)+1)
        pp.shell(k).pair(p).V(j,:) = sscanf(fgetl(fid),'%f');
        end
        end
    end
end



