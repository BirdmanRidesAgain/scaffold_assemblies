#!/bin/bash
################################
# INTRO
################################
# This is a script to run scaffolding and genome visualization
# It renames files and then runs ragtag on them

# If no arguments are given, give usage and exit
USAGE="Usage: driver_rename_fa.sh <assembly.fa> <refseq.fa>"
OUTPUTS="Outputs: a fasta file"
if [ $# -ne 2 ]
then
    echo -e "$USAGE\n${OUTPUTS}"
    exit 1
fi

################################
# HELPER FUNCTIONS
################################
function check_if_file_exists {
  # usage: pass it a name of a file. If it exists, return 0
    if [ -f $1 ]
    then
        echo "$1 already exists. Skipping."
        return 0
    fi
    echo "Output $1 does not exist. Creating."
    return 1
}

################################
# BODY
################################
# Get the prefix for your fasta file
PREFIX_Q=$(basename -s .fa $1)
PREFIX_REF=$(basename -s .fa $2)


##### RENAME FILES #####
# Write a fasta file to the command line and check the FA with no prefix argument
check_if_file_exists "${PREFIX_Q}_rename.fa"
if [ $? -ne 0 ] 
then
    ./rename_fa.py -f $1 -o "${PREFIX_Q}_rename"
fi

# Check to see if your reference sequence has been renamed already. If not, rename it again.
check_if_file_exists "${PREFIX_REF}_rename_refseq.fa"
if [ $? -ne 0 ]
then
    ./rename_fa.py -f $2 -o "${PREFIX_REF}_rename_refseq" -p chrom
fi

##### RUN RAGTAG #####
# Ragtag your genome against references
PAF_LOC="${PREFIX_Q}_${PREFIX_REF}_ragtag/ragtag.scaffold.asm.paf"
check_if_file_exists $PAF_LOC
if [ $? -ne 0 ]
then
    # ragtag.py expected to be in your path. NOT activated via conda
    RAGTAG_PATH=./RagTag/ragtag.py
    $RAGTAG_PATH scaffold "${PREFIX_REF}_rename_refseq.fa" "${PREFIX_Q}_rename.fa" -o ${PREFIX_Q}_${PREFIX_REF}_ragtag -C # '-C' cats all the unplaced contigs and renames them 'chr0'
fi


##### VISUALIZE ALIGNMENT #####
# Requires R as a dependency
./get_contact_map.R  --paf "${PREFIX_Q}_${PREFIX_REF}_ragtag/ragtag.scaffold.asm.paf" --prefix "${PREFIX_Q}_${PREFIX_REF}_ragtag_dotplot"
