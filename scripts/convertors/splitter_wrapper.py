'''
Created on Jul 01, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.

Splits one input sequence into smaller subsequences using EMBOSS splitter. 
'''

import os, sys
from optparse import OptionParser
import subprocess

SOFTNAME = 'splitter'

### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------
def splitter(sequence, size, outdir):
    # create outdir if it does not exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # check splitter installed
    check_process =  subprocess.Popen(['which %s' % SOFTNAME], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    check_process_out, check_process_err = check_process.communicate()
    check_retval = check_process.poll()
    if not check_retval == 0:
        print "Error: software %s is not installed, please install." % SOFTNAME
        sys.exit(1)
    
    # run splitter
    cmd = "%s -sequence fasta::%s -outseq sequence -osformat2 fasta -osdirectory2 %s -ossingle2 Yes -size %s" % (SOFTNAME, sequence, outdir, size)
    process = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process_out, process_err = process.communicate()
    retval = process.poll()
    if not retval == 0:
        print "Error executing %s command: %s" % (cmd.split()[0], process_err)
        sys.exit(1)
    else:
        print process_out


### ---------------------------------------------------------------------------
def renameFiles(size, outdir):
    # Rename new sequence files for job array
    # from seq_name_1-50000.fasta to be seq_name_1.faa
    # from seq_name_50001-10000.fasta to be seq_name_2.faa
    # ...
    for seq_file in os.listdir(outdir):
        if not '.fasta' in seq_file:
            continue
        seq_range = seq_file.split('.')[0].split('_')[-1]
        seq_num = ((int(seq_range.split('-')[0]) - 1) / size) + 1
        new_seq_file = seq_file.replace("%s" % seq_range, "%s" % seq_num)
        os.rename("%s/%s" % (outdir, seq_file), "%s/%s" % (outdir, new_seq_file))

    
### ---------------------------------------------------------------------------
def doRun():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("--sequence", metavar="FILE", help="Input Fasta sequence FILE", action="store", type="string", dest="sequence")
    parser.add_option("--size", metavar="VALUE", help="Size to split at [default=50000]", action="store", type="int", dest="size", default=50000)
    parser.add_option("--outdir", metavar="DIR", help="Output directory", action="store", type="string", dest="outdir")
    (options, args) = parser.parse_args()

    if not (options.sequence and options.outdir):
        parser.print_help()
        sys.exit(0)


    # Run the conversion 
    if os.path.exists(options.sequence):
        # Split sequence into chunks in different files using EMBOSS splitter
        splitter(options.sequence, options.size, options.outdir)
        # Rename files to be ready for use by LSF job array
        renameFiles(options.size, options.outdir)
    else:
        print "Sequence file %s does not exist" % options.sequence
        
if __name__ == '__main__':
    doRun()
