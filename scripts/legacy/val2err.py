'''
Created on Feb 24, 2010
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
    parser.add_option("-o", metavar="EXTENSION", help="val EXTENSION", action="store", type="string", dest="old")
    parser.add_option("-n", metavar="EXTENSION", help="err EXTENSION", action="store", type="string", dest="new")
    parser.add_option("--extract", help="Extract ERRORs only", action="store_true", dest="extract")

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
            if options.extract:
                cmd = "grep ERROR %s > %s" % (oldfile, newfile)
                try:
                    util.runProcess(cmd)
                except Exception, e:
                    print "Error to extract %s" % oldfile
                    print e
    if not options.extract:
        print "To perform the action, please use --extract"
            
### ---------------------------------------------------------------------------
if __name__ == '__main__':
    main()