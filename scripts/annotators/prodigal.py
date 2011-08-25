'''
Created on Jan 26, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import os
from optparse import OptionParser
from genepy import util


### ---------------------------------------------------------------------------
### Logging configuration
### ---------------------------------------------------------------------------
import logging
import logging.config
logging.config.fileConfig('%s/conf/logging.conf' % os.path.realpath(os.path.dirname(__file__)).replace('/annotators', ''))
log = logging.getLogger('genepy.annotators.prodigal')


### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------
def checkValidInput(input_file, common_name):
    """
    Check if the input fasta sequence file is not splitted by scaffolds.
    Run EMBOSS union before if more than on '>' is found
    """
    try:
        util.checkFile(input_file)
        cmd = "grep '>' %s | wc -l" % input_file
        result = util.runProcess(cmd)
        if int(result) > 1:
            new_input_file = "%s.fna" % common_name
            util.checkSoft("union")
            util.checkSoft("descseq")
            cmd_union = "union -sequence %s -stdout Yes -auto Yes | descseq -filter Yes -name '%s' -auto Yes > %s" % (input_file, common_name, new_input_file)
            util.runProcess(cmd_union)
            return new_input_file
        else:
            return input_file
    except util.UtilException, ue:
        raise ue
    except Exception, e:
        raise e

### ---------------------------------------------------------------------------
def convertToTab(result_file, common_name):
    try:
        tab_file = "%s.prodigal.tab" % common_name
        util.checkFile(result_file)
        f_input = open (result_file, 'r')
        f_output = open (tab_file, 'w')
        for line in f_input:
            line = line.strip()
            #      CDS             complement(14682..18617)
            values = line.split()
            location = values[1]
            
            f_output.write("FT   CDS             %s\n" % location)
            f_output.write("FT                   /colour=4\n")
            f_output.write("FT                   /method=\"PRODIGAL\"\n")
        f_input.close()
        f_output.close()
        return tab_file
    except util.UtilException, ue:
        raise ue
    except Exception, e:
        raise e

### ---------------------------------------------------------------------------
def doRun():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-o", metavar="NAME", help="organism common name", action="store", type="string", dest="name")
    parser.add_option("-i", metavar="FILE", help="input organism sequence file in FASTA format", action="store", type="string", dest="input")
    (options, args) = parser.parse_args()

    try:
        common_name = options.name
        input_file = checkValidInput(options.input, common_name)
        output_file = "%s.prodigal" % common_name
    
        # Print info
        log.info("Running prodigal on %s\n" % common_name)
        log.info("Getting sequence from %s\n" % input_file)
    
        # Run prodigal
        softname = "prodigal"
        util.checkSoft(softname)
        cmd = "%s < %s > %s" % (softname, input_file, output_file)
        util.runProcess(cmd)
    
        # Run the conversion only if successful 
        if os.path.exists(output_file):
            # Convert output results into a feature table EMBL file.
            tab_file = convertToTab(output_file, common_name)
        
            # Tidy up
            util.rmFile(common_name + ".fna")
            util.rmFile(output_file)
    
            log.info("%s is the final feature table Prodigal predictions\n" % tab_file)
        else:
            log.info("%s file does not exists\n" % output_file)
    except Exception, e:
        log.error(e)
        raise e

if __name__ == '__main__':
    doRun()
