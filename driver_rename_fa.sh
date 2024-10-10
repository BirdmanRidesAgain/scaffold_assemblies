#!/bin/bash

# If no arguments are given, give usage and exit
USAGE="Usage: driver_rename_fa.sh <assembly.fa>"
OUTPUTS="Outputs: a fasta file"
if [ -z $1 ]
then
    echo -e "$USAGE\n${OUTPUTS}"
    exit 1
fi

# Get the prefix for your fasta file
PREFIX=$(basename -s .fa $1)
echo $PREFIX
# Write a fasta file to the command line and check the FA
./rename_fa.py -f $1 -o "${PREFIX}_rename"


# Testing output
echo ""
echo "# The first line of your output should equal '>contig_1'"
cat "${PREFIX}_rename.fa" | grep '>' | head -n 1 
exit


#COUNT=1
#while read -r i
#do
#    LINE=$(echo $i | awk -v VAR=$COUNT '$2=VAR')
#    echo $LINE | awk '{print $1, "\t", "contig_"$2}' >> ${PREFIX}.tsv
#    COUNT=$((${COUNT}+1))
#done < sort2.tsv
#
## rerun rename_fa.py to generate the new .fa
#./rename_fa.py -f $1 -n ${PREFIX}.tsv
#
##rm sort.tsv #${PREFIX}.tsv
