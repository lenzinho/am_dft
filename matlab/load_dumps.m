% load all files into workspace
tiny = 1E-6;
if ~exist('dumpdir','var')
    dumpdir='.';
end
f = dir([dumpdir,'/outfile.*']);
for i = 1:length(f)
    str = sprintf('%s = load_dump(''%s'');', strrep(f(i).name,'outfile.',''),[dumpdir,'/',f(i).name]);
    eval(str)
end
clear f i str;
clear dumpdir;