'''
Created on Mar 10, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import sys, os
from optparse import OptionParser
import util
import logsetup

### ---------------------------------------------------------------------------
### Logging configuration
### ---------------------------------------------------------------------------
import logging
log = logging.getLogger('genepy.submit')

### ---------------------------------------------------------------------------
### Submit constants
### ---------------------------------------------------------------------------
SUBMIT_SCRIPTS = ['annotated_genome', 'genome_project']

### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", "--list", metavar="FILE", help="FILE containing the list of all organism common names and its associated parameters depending on submitter type", action="store", type="string", dest="list")
    parser.add_option("-r", "--run", metavar="SCRIPT", help="name of the script to run from %s against each genome of the list" % SUBMIT_SCRIPTS, action="store", choices=SUBMIT_SCRIPTS, dest="run")
    parser.add_option("--submit", help="To submit data, not only for checking", action="store_true", dest="submit")
    
    (options, args) = parser.parse_args()

    # Print help if no argument given
    if util.printHelp(options):
        parser.print_help()
        sys.exit()

    # Print command line
    cmdline = "$ python "
    for argv in sys.argv:
        cmdline += argv + " " 
    log.info(cmdline)
    
    script = options.run
    try:
        if options.list:
            util.checkFile(options.list)
            if script == 'genome_project':
                import submitters.genome_project as genome_project
                genome_project.doRun(options.list, options.submit)
            elif script == 'annotated_genome':
                import submitters.annotated_genome as annotated_genome
                annotated_genome.doRun(options.list, options.submit)
        else:
            log.info("Organism list file not provided! Please provide one using -l FILE or --list=FILE")

    except Exception, e:
        import traceback
        log.error(traceback.extract_stack())
        log.error(e)
        
if __name__ == '__main__':
    main()

