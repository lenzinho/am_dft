#!/usr/local/bin/julia

# setup original directory (from where the script is run)
odir = pwd()
# setup base directory: vasp run folder
bdir = ARGS[1]
@printf("dir: %s\n",bdir)
# setup directories
src = odir
dst = bdir
rel = relpath( odir , bdir )
# prepare input
for file in ["/save.pp","/save.tb","/save.ibz"]
    if ~islink(dst*file) 
        symlink(rel*file,dst*file)
    end
end
cp(dst*"/outfile.tb_vsk",dst*"/infile.tb_vsk",remove_destination=true)
# run
cd(bdir)
cmdi=`tbdos -nelecs 8`
cmdo=`tee stdout.tbdos`
run(pipeline(cmdi,cmdo))