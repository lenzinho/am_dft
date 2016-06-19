function [A] = load_dump(fname)
    fprintf('Loading: %-30s ',fname);
    fid=fopen(fname);
        % number of atoms
        shape = sscanf(fgetl(fid),'%i');
        if  isempty(shape)
                A = str2num(fgets(fid));
        elseif (length(shape)==1)
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
        elseif (length(shape)==4)
            for k = 1:shape(4)
            for j = 1:shape(3)
            for i = 1:shape(1)
                A(i,:,j,k) = str2num(fgets(fid));
            end
            end
            end
        end
    fclose(fid);
    fprintf(' ... done\n');
end