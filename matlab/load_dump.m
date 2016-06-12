function [A] = load_dump(fname)
    fid=fopen(fname);
    fprintf('Loading from %s\n',fname);
        % number of atoms
        shape = sscanf(fgetl(fid),'%i');
        if  (length(shape)==1)
            for i = 1:shape(1)
                A(i) = str2num(fgets(fid));
            end
        elseif (length(shape)==2)
            for i = 1:shape(1)
                A(i,:) = str2num(fgets(fid));
            end
        elseif (length(shape)==3)
            for j = 1:shape(3)
            for i = 1:shape(1)
                A(i,:,j) = str2num(fgets(fid));
            end
            end
            
        end
    fclose(fid);
    fprintf(' ... done\n');
end