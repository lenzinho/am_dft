clear;clc
tiny = 1E-3;
% load all files into workspace
f = dir('outfile.*');
for i = 1:length(f)
    str = sprintf('%s = load_dump(''%s'');', strrep(f(i).name,'outfile.',''),f(i).name);
    eval(str)
end
clear f i str;