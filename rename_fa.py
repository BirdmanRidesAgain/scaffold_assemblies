#!/usr/bin/env python3.11
"""
    reverse_compliment.py
    12 Jun 2024

    Copyright (C) 2024-2025 Keiler Collier
    This is a script which takes a fasta file (multiline or oneline) and a TSV.
    The TSV has old contig names on the left and new ones on the right, separated by a tab character.
    This script goes in and renames the contigs in the fasta to fit the new names.
    If a contig is present in the fasta but not in the TSV, it is ignored by default.
        However, if you enable the -r (--remove) flag, it removes those contigs from the assembly.
    
    It will also support pattern renaming - this is mostly useful for our uncurated assemblies.
    

"""
#--------------------------------------------
# Import statements
#--------------------------------------------
import argparse
import os
import pandas as pd

#--------------------------------------------
# Classes and helper functions
#--------------------------------------------
from fa_classes import *

#--------------------------------------------
# MAIN:
#--------------------------------------------    
def main():
    # PARSE ARGS:
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help = "Supply a .fa file")
    parser.add_argument("-p", "--prefix", default = "contig", help = "Supply a prefix to rename contigs. Default is 'contig'")
    parser.add_argument("-n", "--names", default = None, help = "Supply a TSV file of contig names. Should be formatted as <old_name>\t<new_name>. Unnamed seqs will be renamed as per '-s'")
    parser.add_argument("-o", "--out_pre", default = "stdout", help = "Writes output to path provided. If unset, output written to stdout")
    args = parser.parse_args()
    # set default of 'sequential' to none. if a person sets the flag with a prefix, it sets a bool in the program to turn on seq. mode
    
    #--------------------------------------------
    # Sanity-check for help and fasta arguments
    #--------------------------------------------
    
    if (args.fasta == None):
        print("Required inputs not found:")
        print("\tPlease provide a fasta (-f or --fasta)")
        return 1
    if (not os.path.exists(args.fasta)):
        print("Required inputs not found:")
        if (not os.path.exists(args.fasta)):
            print("\t'" + args.fasta + "' does not exist")
        return 2
    
    
    #--------------------------------------------
    # parse a fasta and give us a list of Contigs.
    #--------------------------------------------
    ctg_lst = parse_fasta(args.fasta)
    ctg_lst.sort(reverse=True, key=sortSeq_lst) # sort the contig list by size


    #--------------------------------------------
    # Scan for 'unplaced' contigs/chroms and remove from the list
    #--------------------------------------------
    unplaced_keywords = ['unanchor', 'unplace', 'unscaffold']
    unplaced_ctgs = []
    for i in ctg_lst:
        for j in unplaced_keywords:
            if j in i.name:
                unplaced_ctgs.append(i) # gets a list of unplaced scaffolds probably
    
    # Use list comprehension to remove contigs
    len(ctg_lst)
    ctg_lst = [i for i in ctg_lst if i not in unplaced_ctgs]
    len(ctg_lst) 
    
    #--------------------------------------------
    # If no TSV argument is given, rename all contigs according to 'contig_<num>' and return a fasta
    #--------------------------------------------
    if (args.names == None):
        
        count = 1
        if args.out_pre == 'stdout':
            for i in ctg_lst:       # fasta output
                i.name = args.prefix + "_" + str(count)
                print(i.dump())
                count = count + 1
        else:
            filename = args.out_pre + '.fa'
            with open(filename, 'w') as f:
                for i in ctg_lst:
                    i.name = args.prefix + "_" + str(count)
                    print(i.dump(), file=f)
                    count = count + 1
                f.close()    
                
    return 0 # current end of the main function. We're not screwing with the TSV code.
    print("You should not see this")             
                


######################################################
##### BAD, CURRENTLY-UNREACHABLE CODE BELOW HERE #####
    #--------------------------------------------
    # If a TSV is present, rename the reads according to the TSV.
    #--------------------------------------------


    #--------------------------------------------
    # parse a TSV and give us a list of old/new seqnames
    #--------------------------------------------
    names_df = pd.read_csv(args.names, sep = '\t', names = ["oldname", "newname"])
    names_df["oldname"]
    names_df["newname"]
    
    #--------------------------------------------
    # Replace contigs found in 'oldname' with 'newname' values
    #--------------------------------------------
    if args.sequential != None:
        unnamed_ctr = 0 # this is a counter for any unnamed contigs. Used when -s is enabled

    for i in ctg_lst:
        if i.name in list(names_df["oldname"]):
            # if it's in oldnames, take the corresponding 'oldname' and get its rownum
            idx = (list(names_df["oldname"]).index(i.name))
            newname = names_df["newname"].iloc[idx] # use the rownum to get the 'newname' value
            i.name = newname # rename with newname

        else:
            if (args.sequential != None):
                unnamed_ctr = unnamed_ctr + 1
                i.name = args.sequential + "_" + str(unnamed_ctr)
            elif (args.remove): #FIXME - find a way to get rid of contigs
                i.name = None

    #--------------------------------------------
    # Write the fasta object to a file (or stdout)
    #--------------------------------------------
    # Print to stdout
    if args.out_pre == 'stdout':
        for i in ctg_lst:       # fasta output
            if i.name == None and args.remove == True:
                del(i)
                continue
            print(i.dump())  
    # Print to files
    else:
        filename = args.out_pre + ".fa" # fasta output
        with open(filename, 'w') as f:
            for i in ctg_lst:
                if i.name == None and args.remove == True:
                    del(i)
                    continue
                print(i.dump(), file=f)
            f.close()
    return 0 # End of the main function
    
if __name__ == "__main__":
    main()
