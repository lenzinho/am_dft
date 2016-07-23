function [V] = load_tb_matrix_elements_irre(fname)
    % fname = 'outfile.tb_matrix_elements_irreducible'
    fid=fopen(fname);
        fgetl(fid); % header
        fprintf('loading irreducible matrix elements:\n');
        for i = 1:1000
            buffer = fgetl(fid);
            word   = strsplit(strtrim(buffer));
            j      = sscanf(word{1},'%i');
            V(j)   = sscanf(word{2},'%f');
            if (feof(fid))
                break
            end
        end
        fprintf(' ... matrix elements = %i\n',j);
        fprintf(' ... done\n');
    fclose(fid);
end



