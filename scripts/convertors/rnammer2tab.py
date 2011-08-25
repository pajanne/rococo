'''
Created on Jul 15, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import os, sys
from optparse import OptionParser

### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------
"""
##gff-version2
##source-version RNAmmer-1.2
##date 2010-07-15
##Type DNA
# seqname           source                      feature     start      end   score   +/-  frame  attribute
# ---------------------------------------------------------------------------------------------------------
Alistipes_shahii_WAL8301	RNAmmer-1.2	rRNA	883349	886219	2795.1	+	.	23s_rRNA	
Alistipes_shahii_WAL8301	RNAmmer-1.2	rRNA	886334	886445	83.1	+	.	5s_rRNA	
Alistipes_shahii_WAL8301	RNAmmer-1.2	rRNA	881290	882801	1566.8	+	.	16s_rRNA	
# ---------------------------------------------------------------------------------------------------------

"""
def convertToTab(file, tab_file):

    f_input = open (file, 'r')
    f_output = open (tab_file, 'w')
    for line in f_input:
        line = line.strip()
        values = line.split()
        if line[0] == '#':
            continue
        if len(values) == 0:
            continue
        source = values[1]
        start = int(values[3])
        end = int(values[4])
        strand = values[6]
        feature = values[8]

        if strand == '+':
            location = "%s..%s" % (start, end)
        else:
            location = "complement(%s..%s)" % (end, start)


        f_output.write("FT  rRNA             %s\n" % location)
        f_output.write("FT                   /note=\"%s\"\n" % feature)
        f_output.write("FT                   /method=\"%s\"\n" % source)
    f_input.close()
    f_output.close()


### ---------------------------------------------------------------------------
def doRun():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", metavar="FILE", help="RNAmmer input FILE", action="store", type="string", dest="input")
    parser.add_option("-o", metavar="FILE", help="RNAmmer output FILE converted into EMBL feature tab", action="store", type="string", dest="output")
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
