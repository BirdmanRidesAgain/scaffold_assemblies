#!/bin/bash

# If no arguments are given, give usage and exit
USAGE="Usage: test_rename_fa.sh <assembly.fa>"
OUTPUTS="Outputs: Pass/fail information for various tests"
if [ -z $1 ]
then
    echo -e "$USAGE\n${OUTPUTS}"
    exit 1
fi

##############
# TEST BENCH #
##############
# TEST1 - Output renamed fa to stdout
TEST1=$(./rename_fa.py -f $1 | grep '>' | head -n 1) # ignore the broken pipe error
if [ "$TEST1" != ">contig_1" ]
then
    echo 'TEST1 failed.'
    exit
fi

# TEST2 - Output renamed fa to file
./rename_fa.py -f $1 -o 'TESTFA'
TEST2=$(cat TESTFA.fa | grep '>' | head -n 1)

if [ "$TEST2" != ">contig_1" ]
then
    echo 'TEST2 failed.'
    exit
fi
rm TESTFA.fa

# TEST3 - Output renamed fa to file with a nonstandard contig prefix ('-p' option)
./rename_fa.py -f $1 -o 'TESTFA' -p chrom
TEST3=$(cat TESTFA.fa | grep '>' | head -n 1)

if [ "$TEST3" != ">chrom_1" ]
then
    echo 'TEST3 failed.'
    exit
fi
rm TESTFA.fa

###################
# ENDING BEHAVIOR #
###################
echo "All tests passed"
exit

# Write a fasta file to the command line and check the FA with no prefix argument
#./rename_fa.py -f $1 -o "${PREFIX}_rename"
#echo ""
#echo "# The first line of your output should equal '>contig_1'"
#cat "${PREFIX}_rename.fa" | grep '>' | head -n 1 
#
## Write a fasta file to the command line and check the FA with a prefix argument
#./rename_fa.py -f $1 -o "${PREFIX}_rename_refseq" -p chrom
#echo ""
#echo "# The first line of your output should equal '>chrom_1'"
#cat "${PREFIX}_rename_refseq.fa" | grep '>' | head -n 1 
#
exit

# Ragtag your genome against references
    # ragtag.py expected to be in your path. NOT installed via conda


#### Testing output
#echo ""
#echo "# The first line of your output should equal '>contig_1'"
#cat "${PREFIX}_rename.fa" | grep '>' | head -n 1 


# We now need to scaffold the renamed assembly to 
