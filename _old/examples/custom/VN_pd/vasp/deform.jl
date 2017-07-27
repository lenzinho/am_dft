#!/usr/local/bin/julia

def=7
maxdef=0.05
ndefs=10

# run shell command
cmd = `uc -strain $def $maxdef $ndefs`
run(cmd)
# setup original directory (from where the script is run)
odir = pwd()
# setup base directory containing output of shell command
bdir = odir * "/out/uc"
# list of vasp fils
vlist = ["INCAR","KPOINTS","POTCAR","QSUB"]
# list of deformation files
flist = readdir(bdir)
# loop over deformation files
for f in flist
	@printf "POSCAR: %s\n" f
	# get deformation numbr
	i = replace(f, "outfile.def_", "")
	@printf " ... i: %s\n" i
	# get working directory name
	wdir = odir * "/def/" * i
	@printf " ... wdir: %s\n" wdir
	# create working directory
	mkpath(wdir)
	# copy poscar
	src = bdir * "/" * f
	dst = wdir * "/" * "POSCAR"
	@printf " ... copying: %s\n" src
	@printf " ...          %s\n" dst
	cp( src, dst, remove_destination=true, follow_symlinks=true )
	# create symbolic links
	for v in vlist
		src = odir * "/"
		dst = wdir * "/"
		rel = relpath( src, dst )
		@printf " ... linking: %s\n" src
		@printf "              %s\n" dst
		@printf "              %s\n" rel
		symlink( rel * "/" * v , dst * "/" * v )
	end
end
