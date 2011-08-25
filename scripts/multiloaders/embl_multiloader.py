#!/usr/bin/env python
# encoding: utf-8
'''
Created on Jan 5, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import os
from genepy import util
import sys
import loaders.db as db
import ropy.util
import ropy.query
from optparse import OptionParser
from ropy.log import LogSetup

### ---------------------------------------------------------------------------
### Logging setup
### ---------------------------------------------------------------------------
logsetup = LogSetup()
logsetup.logname = "genepy-multiloader"
logsetup.logpath = "%s/logs.txt" % os.path.realpath(os.path.dirname(__file__))
logsetup.setupLogging()

logger = logsetup.logger

### ---------------------------------------------------------------------------
### main
### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", "--list", metavar="FILE", help="FILE containing the list of all organism common names and its associated file to load", action="store", type="string", dest="list")
    parser.add_option("-D", action="store", dest="dbhost")

    (options, args) = parser.parse_args()
    
    # Print help if no argument given
    if util.printHelp(options):
        parser.print_help()
        sys.exit()
    
    # Print command line
    cmdline = "$ python "
    for argv in sys.argv:
        cmdline += argv + " " 
    logger.info(cmdline)
    
    # Print logger file info
    logger.info(logsetup.logpath)
    
    # Setup database connection
    host = ropy.util.getDArg("dbhost", raiseOnEmpty = True)
    database = ropy.util.getDArg("dbname", raiseOnEmpty = True)
    port = ropy.util.getDArg("dbport", raiseOnEmpty = True)
    user = ropy.util.getDArg("dbuser", raiseOnEmpty = True)
    #password = ropy.util.getDArg("dbpassword", raiseOnEmpty = True)
    
    # Check if chado_load is installed
    util.isSoftInstalled("chado_load")

    # Read organism common name and load related embl file into the database
    data_path = options.list
    for line in open(data_path, "r"):
        if line[0] == '!':
            continue
        if line.count('||') < 1:
            continue
        # ! common_name||taxon_id||filename
        line = line.strip()
        list = line.split('||')
        common_name = list[0]
        filename = list[2]
        util.checkFile(filename)
        # Loader command
        cmd = "chado_load embl -o %s -t contig -D %s:%s/%s -U %s %s" % (common_name, host, port, database, user, filename)
        # Run command
        util.runProcess(cmd)

### ---------------------------------------------------------------------------


### ---------------------------------------------------------------------------
if __name__ == '__main__':
    main()
