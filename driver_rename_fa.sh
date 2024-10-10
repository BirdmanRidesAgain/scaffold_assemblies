#!/bin/bash

# If no arguments are given, give usage and exit
USAGE="Usage: driver_rename_fa.sh <assembly.fa> <refseq.fa>"
OUTPUTS="Outputs: a fasta file"
if [ $# -ne 2 ]
then
    echo -e "$USAGE\n${OUTPUTS}"
    exit 1
fi

# Get the prefix for your fasta file
PREFIX_Q=$(basename -s .fa $1)
PREFIX_REF=$(basename -s .fa $2)



##### RENAME FILES #####
# Write a fasta file to the command line and check the FA with no prefix argument
./rename_fa.py -f $1 -o "${PREFIX_Q}_rename"

# Check to see if your reference sequence has been renamed already. If not, rename it again.
OLD_REF_CHROMNAME=$(cat $2 | head -n 1)
if [ "$OLD_REF_CONTIG_NAME" != ">chrom_1" ]
then
    echo "Renaming the reference sequence. Chromosomal assembly assumed."
    ./rename_fa.py -f $2 -o "${PREFIX_REF}_rename_refseq" -p chrom
fi

##### RUN RAGTAG #####
# Ragtag your genome against references
    # ragtag.py expected to be in your path. NOT activated via conda
ragtag.py scaffold "${PREFIX_REF}_rename_refseq.fa" "${PREFIX_Q}_rename.fa" -o ${PREFIX_Q}_${PREFIX_REF}_ragtag -C # '-C' cats all the unplaced contigs and renames them 'chr0'
exit

#### Testing output
#echo ""
#echo "# The first line of your output should equal '>contig_1'"
#cat "${PREFIX}_rename.fa" | grep '>' | head -n 1 


# We now need to scaffold the renamed assembly to 
