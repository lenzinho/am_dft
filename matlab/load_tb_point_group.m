function [pg] = load_tb_point_group()
    fid=fopen('outfile.tb_pointgroup');
    % number of atoms
    pg.nsyms = sscanf(fgetl(fid),'%i');
    fprintf('symmetries: \n',pg.nsyms);
    pg.nbases = sscanf(fgetl(fid),'%i');
    pg.sym(1:pg.nbases,1:pg.nbases,1:pg.nsyms) = 0;
    for j = 1:pg.nsyms
        i = sscanf(fgetl(fid),'%i');
        fprintf(' %3i',i);
        fgetl(fid); % divider
        for k = 1:pg.nbases
        pg.sym(k,:,i) = sscanf(fgetl(fid),'%f');
        end
        if (mod(j,8)==0) 
            fprintf('\n');
        elseif j == pg.nsyms
            fprintf('\n');
        end
    end
    fclose(fid);
end



