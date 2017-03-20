#!/bin/bash
# Parses OUTCAR extracting atomic forces and displacements from each timetep.
# Antonio Mei Nov/2014
# Antonio Mei Jan/2017
usage_ () {
    echo "Creates infile.force_position based on supplied outcar files."
    echo ""
    echo "Usage: $0 [-h] [-t] [-f] -o <outcar_list> -n <natoms> [-c <compress_name>]"
    echo ""
    echo "Example: $0 -f -t -o \"\$(find . -name OUTCAR | grep 4.00-300)\" -n 250"
    echo ""
    echo "-h : prints this message"
    echo "-n : [REQUIRED] number of atoms"
    echo "-o : [REQUIRED] list of outcar files to parse"
    echo "-t : trims the last md run (useful for removing runs which have not completed)"
    echo "-f : overwrites existing infile.force_position"
    echo "-c : compresses infile.force_position to a tar.gz file"
    echo ""
    echo "infile.force_position file contents:"
    echo "   x position   y position   z position     x force      y force      z force"
    exit 1
}

main_ () {
    # trim the last md run which may not have completed
    trim_ () { tac $1 | awk '!found && /POSITION/{found=1;next}1' | tac ; }
    # get position and forces
    get_  () { cat $2 | grep -h -A $(($1+1)) POSITION  ; }
    # cut header lines
    cut_  () { cat $1 | sed '/^--$/d' | sed '/--/d' | sed '/POSITION/d' ; }
    # compress produced infile.force_position
    compress_ () { tar -zcvf infile.force_position.tar.gz infile.force_position ; }
    #
    if ${ISFORCE}; then
        if [ -f "./infile.force_position" ]; then
            rm ./infile.force_position
            printf " ... ./infile.force_position overwritten\n"
        fi
    fi
    #
    if ${ISTRIM}; then
        printf " ... trim:\n"
        for F in "${FLIST}"; do
            printf " ...     %-100s\n" "${F}"
            trim_ ${F} | get_ ${NATOMS} | cut_ >> infile.force_position
        done
    else
        printf " ... batch parsing without trim\n"
        get_ ${NATOMS} "${FLIST}" | cut_ >> infile.force_position
    fi
    #
    printf " ... infile.force_position created\n"
    #
    if ${ISCOMPRESS}; then
        printf " ... infile.force_position.tar.gz compressed\n"
        compress_ 
    fi
}

ISCOMPRESS=false; ISTRIM=false; ISFORCE=false;
if (($# == 0)); then usage_; exit 1; fi
while getopts "n:o:htfc" o; do
    case "${o}" in
        o)  FLIST=${OPTARG} ;;
        n)  NATOMS=${OPTARG} ;;
        c)  ISCOMPRESS=true ;;
        t)  ISTRIM=true ;;
        f)  ISFORCE=true ;;
        h)  usage_; exit 0 ;;
        *)  usage_; exit 1 ;;
    esac
done

main_
