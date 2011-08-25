'''
Created on Jun 24, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import os, sys
from optparse import OptionParser

### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------
def convertToTab(g3_file, tab_file):

    f_input = open (g3_file, 'r')
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


### ---------------------------------------------------------------------------
def doRun():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", metavar="FILE", help="Prodigal input FILE", action="store", type="string", dest="input")
    parser.add_option("-o", metavar="FILE", help="Prodigal output FILE converted into EMBL feature tab", action="store", type="string", dest="output")
    (options, args) = parser.parse_args()

    if not (options.input and options.output):
        parser.print_help()
        sys.exit()


    # Run the conversion 
    if os.path.exists(options.input):
        # Convert output results into a feature table EMBL file.
        convertToTab(options.input, options.output)
        
if __name__ == '__main__':
    doRun()
