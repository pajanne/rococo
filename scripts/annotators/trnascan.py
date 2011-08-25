'''
Created on Jan 27, 2010
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
log = logging.getLogger('genepy.annotators.trnascan')


### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------
def convertToTab(result_file, common_name):
    try:
        tab_file = "%s.trna.tab" % common_name
        util.checkFile(result_file)
        f_input = open (result_file, 'r')
        f_output = open (tab_file, 'w')
        for line in f_input:
            line = line.strip()
            # name=Alistipes_shahii 	1	89017	88933	Undet	???	0	0	61.20
            values = line.split()
            start = int(values[2])
            end = int(values [3])
            trna_type = values[4]
            anti_codon = values[5]
            score = values [8]

            if start <= end:
                location = "%s..%s" % (start, end)
            else:
                location = "complement(%s..%s)" % (end, start)

            
            f_output.write("FT  tRNA             %s\n" % location)
            f_output.write("FT                   /note=\"tRNA-%s(%s) Cove Score %s\"\n" % (trna_type, anti_codon, score))
            f_output.write("FT                   /colour=4\n")
            f_output.write("FT                   /method=\"tRNAscan-SE\"\n")
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
        input_file = options.input
        output_file = "%s.trna" % common_name
    
        # Print info
        log.info("Running tRNAscan on %s\n" % common_name)
        log.info("Getting sequence from %s\n" % input_file)
    
        # Run 
        softname = "tRNAscan-SE"
        util.checkSoft(softname)
        cmd = "%s -P -o %s -q -b -Q -C %s" % (softname, output_file, input_file)
        util.runProcess(cmd)
    
        # Run the conversion only if successful 
        if os.path.exists(output_file):
            # Convert output results into a feature table EMBL file.
            tab_file = convertToTab(output_file, common_name)
        
            # Tidy up
            #util.rmFile(output_file)
    
            log.info("%s is the final feature table tRNAscan-SE predictions\n" % tab_file)
            pass
        else:
            log.info("%s file does not exists\n" % output_file)
    except Exception, e:
        log.error(e)
        raise e

if __name__ == '__main__':
    doRun()
