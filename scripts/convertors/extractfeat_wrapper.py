'''
Created on Jul 26, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.

Splits input EMBL sequence based on its features using EMBOSS extractfeat,
translate the DNA sequence into protein using EMBOSS transeq,
and rename the files to enable job array submissions. 
'''

import os, sys
from optparse import OptionParser
import subprocess

SOFTNAME = 'extractfeat'
SOFTNAME_TRANSEQ = 'transeq'

### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------
def extract(sequence, type, outdir):
    # create outdir if it does not exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # check extractfeat installed
    check_process =  subprocess.Popen(['which %s' % SOFTNAME], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    check_process_out, check_process_err = check_process.communicate()
    check_retval = check_process.poll()
    if not check_retval == 0:
        print "Error: software %s is not installed, please install." % SOFTNAME
        sys.exit(1)
    
    # check transeq installed
    check_process =  subprocess.Popen(['which %s' % SOFTNAME_TRANSEQ], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    check_process_out, check_process_err = check_process.communicate()
    check_retval = check_process.poll()
    if not check_retval == 0:
        print "Error: software %s is not installed, please install." % SOFTNAME_TRANSEQ
        sys.exit(1)

    # run extractfeat & transeq
    cmd = "%s -sequence embl::%s -type %s -featinname YES -stdout Yes -auto Yes | %s  -filter Yes -outseq fasta:: -ossingle2 Yes -osdirectory2 %s" % (SOFTNAME, sequence, type, SOFTNAME_TRANSEQ, outdir)
    process = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process_out, process_err = process.communicate()
    retval = process.poll()
    if not retval == 0:
        print "Error executing %s command: %s" % (cmd.split()[0], process_err)
        sys.exit(1)
    else:
        print process_out


### ---------------------------------------------------------------------------
def renameFiles(outdir, name):
    # Rename sequence files for job array and translate
    # from seq_name_481_3516_cds_1.pep to be seq_name_cds_1.pep
    # from seq_name_3497_4207_cds_1.pep to be seq_name_cds_2.pep
    # ...

    seq_num = 0
    for seq_file in os.listdir(outdir):
        if not '.pep' in seq_file:
            continue
        #seq_range_begin = seq_file.split('.')[0].split('_')[-4]
        #seq_range_end = seq_file.split('.')[0].split('_')[-3]
        seq_num += 1
        #new_seq_file = seq_file.replace("_%s_%s" % (seq_range_begin, seq_range_end), "")
        #new_seq_file = new_seq_file.replace("_1.", "_%s." % seq_num)
        new_seq_file = "%s_%s.pep" % (name, seq_num)
        # rename
        os.rename("%s/%s" % (outdir, seq_file), "%s/%s" % (outdir, new_seq_file))

    
### ---------------------------------------------------------------------------
def doRun():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("--sequence", metavar="FILE", help="Input EMBL sequence FILE", action="store", type="string", dest="sequence")
    parser.add_option("--name", metavar="NAME", help="Sequence NAME to be used for naming files", action="store", type="string", dest="name")
    parser.add_option("--type", metavar="NAME", help="Feature type you wish to extract [default='CDS']", action="store", type="string", dest="type", default="CDS")
    parser.add_option("--outdir", metavar="DIR", help="Output directory", action="store", type="string", dest="outdir")
    (options, args) = parser.parse_args()

    if not (options.sequence and options.name and options.outdir):
        parser.print_help()
        sys.exit(0)

    # Run the conversion 
    if os.path.exists(options.sequence):
        # extract sequence
        extract(options.sequence, options.type, options.outdir)
        # Rename files to be ready for use by LSF job array
        renameFiles(options.outdir, options.name)
    else:
        print "Sequence file %s does not exist" % options.sequence
        
if __name__ == '__main__':
    doRun()
