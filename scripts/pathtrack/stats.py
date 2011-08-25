'''
Created on May 18, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.

Statistics report on either a fasta file to display total size and N50 style information
'''

from genepy import logsetup
from optparse import OptionParser
from genepy import util
import sys, os

### ---------------------------------------------------------------------------
### Logging configuration
### ---------------------------------------------------------------------------
import logging
log = logging.getLogger('genepy.pathtrack')

### ---------------------------------------------------------------------------
def infoseq(file):
    """
    Run EMBOSS infoseq to 
    Display basic information about sequences
    """
    util.checkFile(file)
    util.checkSoft("infoseq")
    cmd = "infoseq -only -length -noheading %s -outfile %s.infoseq" % (file, file)
    util.runProcess(cmd)

### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--file", metavar="FILE", help="Fasta FILE to analyse", action="store", type="string", dest="file")
    
    (options, args) = parser.parse_args()

    if not (options.file):
        parser.print_help()
        sys.exit()
        
    infoseq(options.file)
    infoseq_file = open("%s.infoseq" % options.file, 'r').readlines()
    total_nb_residues = 0
    number_of_sequences = 0
    stat_list = []
    for line in infoseq_file:
        line = line.strip()
        number_of_sequences = number_of_sequences + 1
        total_nb_residues = total_nb_residues + int(line)
        stat_list.append(int(line))
    average_length = total_nb_residues / number_of_sequences 
    stat_list.sort()
    smallest = stat_list[0]
    largest = stat_list[-1]

    stats_file = open("%s.stats" % options.file, 'w')

    # tab delimited output
    stats_file.write("#seq\t#bases\tsmallest\tlargest\tavg\tN50_size\tN50_#seq\n")
    stats_file.write("%s\t%s\t%s\t%s\t%s\t" % (number_of_sequences, total_nb_residues, smallest, largest, average_length))
    
    # N50
    stat_list.reverse()
    n50_sum = 0
    n50_size = 0
    n50_number_of_sequences = 0
    for x in stat_list:
        n50_sum = n50_sum + x
        n50_number_of_sequences = n50_number_of_sequences + 1
        n50_size = x
        if (n50_sum > (total_nb_residues / 2)):
            stats_file.write("%s\t%s\n" % (n50_size, n50_number_of_sequences))
            break

    # clean tmp file
    util.rmFile("%s.infoseq" % options.file)

    log.info("Results in %s.stats" % options.file)
        

### ---------------------------------------------------------------------------
if __name__ == '__main__':
    main()

