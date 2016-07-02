function [dr,bz] = load_vasp_eigenval(fname)
    fprintf(' ... loading dispersion from: %s \n',fname);
    fid=fopen(fname);
    % skip first five lines
    for i = 1:5
    fgetl(fid);
    end
    buffer = strsplit(strtrim(fgetl(fid)));
    dr.nelecs = sscanf(buffer{1},'%i');
    bz.nkpts  = sscanf(buffer{2},'%i');
    dr.nbands = sscanf(buffer{3},'%i');
    fprintf(' ... electrons = %i \n',dr.nelecs);
    fprintf(' ... kpoints = %i \n',bz.nkpts);
    fprintf(' ... bands = %i \n',dr.nbands);
    for i = 1:bz.nkpts
        % skip line
        fgetl(fid);
        % get kpnts
        buffer = strsplit(strtrim(fgetl(fid)));
        bz.kpt(1,i) = sscanf(buffer{1},'%f');
        bz.kpt(2,i) = sscanf(buffer{2},'%f');
        bz.kpt(3,i) = sscanf(buffer{3},'%f');
        % loop over bands
        for j = 1:dr.nbands
            buffer = strsplit(strtrim(fgetl(fid)));
            dr.E(j,i)  = sscanf(buffer{2},'%f');
        end
        dr.E(:,i) = sort(dr.E(:,i));
    end
    fprintf(' ... done\n');
    fclose(fid);
end



