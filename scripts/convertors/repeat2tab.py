'''
Created on Jul 12, 2010
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
   SW   perc perc perc  query                     position in query           matching  repeat           position in repeat
score   div. del. ins.  sequence                  begin   end        (left)   repeat    class/family   begin  end    (left)  ID

  997    0.8  0.0  0.0  Alistipes_shahii_WAL8301      215     339 (3762978) + R=19      Unknown             1    125    (0)   1  
   22    0.0  0.0  0.0  Alistipes_shahii_WAL8301    26000   26021 (3737296) C GC_rich   Low_complexity  (278)     22      1   2  
  516   14.8  0.0  0.0  Alistipes_shahii_WAL8301    31405   31512 (3731805) C R=0       Unknown           (2)    112      5   3  
  499   21.0  0.0  0.0  Alistipes_shahii_WAL8301    41588   41687 (3721630) + R=0       Unknown            13    112    (2)   4  
  684    0.0  0.0  0.0  Alistipes_shahii_WAL8301    63180   63260 (3700057) + R=12      Unknown             1     81    (0)   5  

"""
def convertToTab(file, tab_file):

    f_input = open (file, 'r')
    f_output = open (tab_file, 'w')
    for line in f_input:
        line = line.strip()
        values = line.split()
        if len(values) == 0:
            continue
        if values[0] =='SW':
            continue
        elif values[0] == 'score':
            continue
        start = int(values[5])
        end = int(values[6])
        strand = values[8]
        repeat_type = values[10]

        if strand == '+':
            location = "%s..%s" % (start, end)
        else:
            location = "complement(%s..%s)" % (end, start)


        f_output.write("FT  repeat_region    %s\n" % location)
        f_output.write("FT                   /note=\"Repeat family: %s\"\n" % repeat_type)
        f_output.write("FT                   /method=\"RepeatScout library\"\n")
    f_input.close()
    f_output.close()


### ---------------------------------------------------------------------------
def doRun():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", metavar="FILE", help="RepeatMasker input FILE", action="store", type="string", dest="input")
    parser.add_option("-o", metavar="FILE", help="RepeatMasker output FILE converted into EMBL feature tab", action="store", type="string", dest="output")
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
