'''
Created on Feb 19, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import os, sys
from optparse import OptionParser
from genepy import util

### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-o", metavar="EXTENSION", help="Old EXTENSION", action="store", type="string", dest="old")
    parser.add_option("-n", metavar="EXTENSION", help="New EXTENSION", action="store", type="string", dest="new")
    parser.add_option("--convert", help="Do convert genbank file into embl", action="store_true", dest="convert")

    (options, args) = parser.parse_args()
    
    # Print help if no argument given
    if util.printHelp(options):
        parser.print_help()
        sys.exit()
    
    for file in os.listdir('.'):
        if options.old in file:
            oldfile = file
            newfile = "%s.%s" % (oldfile.split(".")[0], options.new)
            print "Convert file %s into %s" % (oldfile, newfile)
            if options.convert:
                cmd = "seqret -sequence gb::%s -feature Yes -outseq embl::%s" % (oldfile, newfile)
                try:
                    util.runProcess(cmd)
                except Exception, e:
                    print "Error to convert %s" % oldfile
                    print e
    if not options.convert:
        print "To perform the action, please use --convert"
            
### ---------------------------------------------------------------------------
if __name__ == '__main__':
    main()