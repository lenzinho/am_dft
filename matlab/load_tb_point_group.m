function [pg] = load_tb_point_group()
    fid=fopen('outfile.tb_pointgroup');
    % number of atoms
    pg.nsyms = sscanf(fgetl(fid),'%i');
    fprintf('symmeties %i\n',pg.nsyms);
    pg.nbases = sscanf(fgetl(fid),'%i');
    pg.sym(1:pg.nbases,1:pg.nbases,1:pg.nsyms) = 0;
    for j = 1:pg.nsyms
        i = sscanf(fgetl(fid),'%i');
        fprintf(' ... symmetry %3i\n',i);
        fgetl(fid); % divider
        for j = 1:pg.nbases
        pg.sym(j,:,i) = sscanf(fgetl(fid),'%f');
        end
    end
end



