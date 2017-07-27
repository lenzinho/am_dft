#!/bin/bash
# Parses OUTCAR extracting electronic band energies from each timetep.
# Antonio Mei Nov/2014
# Antonio Mei Jan/2017
usage_ () {
    echo "Creates infile.electron_energies based on supplied outcar files."
    echo ""
    echo "Usage: $0 [-h] [-t] [-f] -o <outcar_list> -n <nbands> [-c <compress_name>]"
    echo ""
    echo "Example: $0 -f -t -o \"\$(find . -name OUTCAR | grep 4.00-300)\" -n 751"
    echo ""
    echo "-h : prints this message"
    echo "-n : [REQUIRED] number of bands"
    echo "-o : [REQUIRED] list of outcar files to parse"
    echo "-t : trims the last md run (useful for removing runs which have not completed)"
    echo "-f : overwrites existing infile.electron_energies"
    echo "-c : compresses infile.electron_energies to a tar.gz file"
    echo ""
    echo "infile.electron_energies file contents:"
    echo "   En energy"
    exit 1
}
main_ () {
    # trim the last md run which may not have completed
    trim_ () { tac $1 | awk '!found && /POSITION/{found=1;next}1' | tac ; }
    # get energies
    get_  () { cat $2 | grep -h -A ${1} occupation  ; }
    # cut header lines
    cut_  () { cat $1 | sed '/^--$/d' | sed '/--/d' | sed '/occupation/d' ; }
    # compress produced infile.electron_energies
    compress_ () { tar -zcvf infile.electron_energies.tar.gz infile.electron_energies ; }
    #
    if ${ISFORCE}; then
        if [ -f "./infile.electron_energies" ]; then
            rm ./infile.electron_energies
            printf " ... ./infile.electron_energies overwritten\n"
        fi
    fi
    # 
    if ${ISTRIM}; then
        printf " ... trim:\n"
        for F in "${FLIST}"; do
            printf " ...     %-100s\n" "${F}"
            trim_ ${F} | get_ ${NBANDS} | cut_ >> infile.electron_energies
        done
    else
        printf " ... batch parsing without trim\n"
        get_ ${NBANDS} "${FLIST}" | cut_ >> infile.electron_energies
    fi
    #
    awk '{ print $2 }' infile.electron_energies > infile.electron_energies.tmp && mv infile.electron_energies.tmp infile.electron_energies
    #
    printf " ... infile.electron_energies created\n"
    #
    if ${ISCOMPRESS}; then
        printf " ... infile.electron_energies.tar.gz compressed\n"
        compress_ 
    fi
}
ISCOMPRESS=false; ISTRIM=false; ISFORCE=false;
if (($# == 0)); then usage_; exit 1; fi
while getopts "n:o:htfc" o; do
    case "${o}" in
        o)  FLIST=${OPTARG} ;;
        n)  NBANDS=${OPTARG} ;;
        c)  ISCOMPRESS=true ;;
        t)  ISTRIM=true ;;
        f)  ISFORCE=true ;;
        h)  usage_; exit 0 ;;
        *)  usage_; exit 1 ;;
    esac
done
main_