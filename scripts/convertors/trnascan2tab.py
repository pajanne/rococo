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


### ---------------------------------------------------------------------------
def doRun():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", metavar="FILE", help="tRNAscan-SE input FILE", action="store", type="string", dest="input")
    parser.add_option("-o", metavar="FILE", help="tRNAscan-SE output FILE converted into EMBL feature tab", action="store", type="string", dest="output")
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
