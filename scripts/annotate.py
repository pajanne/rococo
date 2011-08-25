'''
Created on Jan 27, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import sys, os
from optparse import OptionParser
import util

### ---------------------------------------------------------------------------
### Logging configuration
### ---------------------------------------------------------------------------
import logging
import logging.config
logging.config.fileConfig('%s/conf/logging.conf' % os.path.realpath(os.path.dirname(__file__)))
log = logging.getLogger('genepy.annotate')

### ---------------------------------------------------------------------------
### Annotate constants
### ---------------------------------------------------------------------------
ANNOTATOR_SCRIPTS = ['glimmer', 'prodigal', 'trnascan']
ANNOTATOR_EXTENSION = {'glimmer':'g3', 
                       'prodigal':'prodigal',
                       'trnascan':'trna'}
ANNOTATOR_QUEUE = {'glimmer':'small', 
                   'prodigal':'small',
                   'trnascan':'normal'}

### ---------------------------------------------------------------------------
### Main methods
### ---------------------------------------------------------------------------
def checkDoRun(script):
    file = open("%s/annotators/%s.py" % (os.path.realpath(os.path.dirname(__file__)), script), 'r')
    if 'def doRun()' not in file.read():
        log.error("doRun() method not found in annotators.%s.py. Please implement it!" % script)
        sys.exit(1)
    
### ---------------------------------------------------------------------------
def mainUsage():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", "--list", metavar="FILE", help="FILE containing the list of all organism common names and its associated FASTA sequence file", action="store", type="string", dest="list")
    parser.add_option("-r", "--run", metavar="SCRIPT", help="name of the script to run from %s against each genome of the list" % ANNOTATOR_SCRIPTS, action="store", choices=ANNOTATOR_SCRIPTS, dest="run")
    
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
    
    return options

### ---------------------------------------------------------------------------
def main():
    options = mainUsage()
    
    script = options.run
    checkDoRun(script)

    try:
        # Read organism common name and related fasta sequence file
        list_file = options.list
        util.checkFile(list_file)
        for line in open(list_file, "r"):
            if line[0] == '!':
                continue
            if line.count('||') < 1:
                continue
            # ! common_name||sequence_file
            common_name, input_file = util.splitLine(line)
            util.checkFile(input_file)
            # Run command
            cmd = 'python -c "from genepy.annotators.%s import doRun; doRun()" -o %s -i %s' % (script, common_name, input_file)
            if util.isLsf():
                job_name = "%s.%s" % (common_name, ANNOTATOR_EXTENSION[script])
                util.submitJob(job_name, cmd, ANNOTATOR_QUEUE[script])
            else:
                util.runProcess(cmd)
    except Exception, e:
        log.error(e)
        
if __name__ == '__main__':
    main()
