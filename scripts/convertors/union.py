'''
Created on Feb 24, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import sys
from optparse import OptionParser
from genepy import util

### ---------------------------------------------------------------------------
def union(file, common_name, locus_tag, organism_name, strain):
    """
    Merge scaffolds into one sequence file
    Run EMBOSS union if more than on '>' is found
    """
    util.checkFile(file)
    cmd = "grep '>' %s | wc -l" % file
    result = util.runProcess(cmd)
    if int(result) > 1:
        new_file = "%s.fsa" % common_name
        util.checkSoft("union")
        util.checkSoft("descseq")
        name = "%s [organism=%s] [strain=%s] [gcode=11]" % (locus_tag, organism_name, strain)
        cmd_union = "union -sequence %s -stdout Yes -auto Yes | descseq -filter Yes -name '%s' -auto Yes > %s" % (file, name, new_file)
        util.runProcess(cmd_union)
        return new_file
    else:
        return file

### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", metavar="FILE", help="FILE containing the list of all organism common names and its associated sequence file", action="store", type="string", dest="list")
    
    (options, args) = parser.parse_args()

    # Print help if no argument given
    if util.printHelp(options):
        parser.print_help()
        sys.exit()
        
    # Get and check input arguments
    if options.list:
        # Read organism common name and related fasta sequence file
        list_file = options.list
        util.checkFile(list_file)
        for line in open(list_file, "r"):
            if line[0] == '!':
                continue
            if line.count('||') < 1:
                continue
            # ! common_name||organim_name||strain||locus_tag||fasta_file
            line = line.strip()
            values = line.split('||')
            print "Processing %s" % values[0]
            union(file=values[4], common_name=values[0], locus_tag=values[3], organism_name=values[1], strain=values[2])           

if __name__ == '__main__':
    main()
