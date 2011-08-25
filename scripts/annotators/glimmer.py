'''
Created on Jan 14, 2010
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
log = logging.getLogger('genepy.annotators.glimmer')


### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------
def checkValidInput(input_file, common_name):
    """
    Check if the input fasta sequence file is of correct format.
    Segmentation fault while running glimmer on splitted sequences with a fasta file
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
def convertToTab(g3_predict_file, common_name):
    try:
        g3_tab = "%s.g3.tab" % common_name
        util.checkFile(g3_predict_file)
        f_input = open (g3_predict_file, 'r')
        f_output = open (g3_tab, 'w')
        for line in f_input:
            if line[0] == '>':
                continue
            line = line.strip()
            # id    start    end    direction    score
            
            values = line.split()
            id = values[0]
            start = int(values[1])
            end = int(values[2])
            direction = values[3]
            score = values[4]
            
            if direction[0] == "+":
                location = "%s..%s" % (start, end)
            else:
                location = "complement(%s..%s)" % (end, start)
    
            if not ((direction[0] == '+' and start > end) or (direction[0] == '-' and start < end)):
    
                f_output.write("FT   CDS             %s\n" % location)
                f_output.write("FT                   /note=\"Raw score %s\"\n" % score)
                f_output.write("FT                   /label=%s\n" % id)
                f_output.write("FT                   /colour=4\n")
                f_output.write("FT                   /method=\"GLIMMER\"\n")
        f_input.close()
        f_output.close()
        return g3_tab
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
    
        # Print info
        log.info("Running Glimmer3 on %s\n" % common_name)
        log.info("Getting sequence from %s\n" % input_file)
    
        # Run glimmer3 iterated
        script = "/software/pathogen/external/applications/glimmer/glimmer/scripts/g3-iterated.csh"
        util.checkFile(script)
        cmd = "%s %s %s" % (script, input_file, common_name)
        util.runProcess(cmd)
    
        # Run the conversion only if g3 successful 
        g3_predict_file = "%s.predict" % common_name
        if os.path.exists(g3_predict_file):
            # Convert output results into a feature table EMBL file.
            g3_tab = convertToTab(g3_predict_file, common_name)
        
            # Tidy up
            util.rmFile(common_name + ".longorfs")
            util.rmFile(common_name + ".train")
            util.rmFile(common_name + ".icm")
            util.rmFile(common_name + ".run1.detail")
            util.rmFile(common_name + ".run1.predict")
            util.rmFile(common_name + ".coords")
            util.rmFile(common_name + ".upstream")
            util.rmFile(common_name + ".motif")
            util.rmFile(common_name + ".detail")
            util.rmFile(g3_predict_file)
    
            log.info("%s is the final feature table Glimmer3 predictions\n" % g3_tab)
        else:
            log.info("%s file does not exists\n" % g3_predict_file)
    except Exception, e:
        log.error(e)
        raise e
        
if __name__ == '__main__':
    doRun()
