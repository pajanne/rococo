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


### ---------------------------------------------------------------------------
def doRun():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", metavar="FILE", help="Glimmer3 input FILE", action="store", type="string", dest="input")
    parser.add_option("-o", metavar="FILE", help="Glimmer3 output FILE converted into EMBL feature tab", action="store", type="string", dest="output")
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
