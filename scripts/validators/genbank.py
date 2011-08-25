'''
Created on Mar 16, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''
import os, sys
import traceback
from optparse import OptionParser
from genepy import util
from Bio import SeqIO
from genepy import logsetup

### ---------------------------------------------------------------------------
### Logging configuration
### ---------------------------------------------------------------------------
import logging
log = logging.getLogger('genepy.validate')

### ---------------------------------------------------------------------------
def doConvert(embl_file, tbl_file, locus_tag):
    # convert .4dep.embl into .4val.embl to get a std ID line
    embl_file_val = embl_file.replace('4dep','4val')
    cmd = "seqret -sequence embl::%s -feature Yes -outseq embl::%s" % (embl_file, embl_file_val)
    util.runProcess(cmd)
    # convert .4val.embl into tbl file
    record = SeqIO.read(open(embl_file_val), "embl")
    table = open(tbl_file, 'w')
    table.write('>Feature %s\n' % locus_tag)
    for feature in record.features:
        if not ('source' in feature.type or 'gap' in feature.type):
            table.write('%s\t%s\n' % (getLocation(feature), feature.type))
            for qualifier in feature.qualifiers:
                if not 'translation' in qualifier:
                    table.write('\t\t%s\t%s\n' % (qualifier, feature.qualifiers[qualifier][0]))
            table.write('\t\tinference\tab initio prediction:IMG/ER\n')
    table.close()

### ---------------------------------------------------------------------------
def getLocation(feature):
    if feature.strand == -1:
        if str(feature.location.start)[0] == '<' or str(feature.location.start)[0] == '>':
            return "%s\t%s%s" % (feature.location.end.position, str(feature.location.start)[0], int(feature.location.start.position) + 1)
        elif str(feature.location.end)[0] == '<' or str(feature.location.end)[0] == '>':
            return "%s%s\t%s" % (str(feature.location.end)[0], feature.location.end.position, int(feature.location.start.position) + 1)
        else:
            return "%s\t%s" % (feature.location.end.position, int(feature.location.start.position) + 1)
    else:
        if str(feature.location.start)[0] == '<' or str(feature.location.start)[0] == '>':
            return "%s%s\t%s" % (str(feature.location.start)[0], int(feature.location.start.position) + 1, feature.location.end.position)
        elif str(feature.location.end)[0] == '<' or str(feature.location.end)[0] == '>':
            return "%s\t%s%s" % (int(feature.location.start.position) + 1, str(feature.location.end)[0], feature.location.end.position)
        else:
            return "%s\t%s" % (int(feature.location.start.position) + 1, feature.location.end.position)
    
### ---------------------------------------------------------------------------
def doValidate():
    cmd = "%s/tbl2asn -p /Users/ap12/Documents/metahit_data/EMBLValidation -t /Users/ap12/Documents/metahit_data/EMBLValidation/template -V v" % os.path.realpath(os.path.dirname(__file__))
    try:
        util.runProcess(cmd)
    except Exception, e:
        log.error(traceback.extract_stack())
        log.error(e)
    # extract errors only
    for file in os.listdir('.'):
        if 'val' in file:
            oldfile = file
            newfile = "%s.err" % (oldfile.split(".")[0])
            cmd = "grep ERROR %s > %s" % (oldfile, newfile)
            try:
                util.runProcess(cmd)
            except Exception, e:
                log.error(traceback.extract_stack())
                log.error(e)

### ---------------------------------------------------------------------------
def doClean():
    pass

### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", metavar="FILE", help="FILE containing the list of all organism common names and its associated information (common_name||organim_name||strain||locus_tag||genome_project_id||coverage)", action="store", type="string", dest="list")
    parser.add_option("--convert", help="Do convert embl file into tbl", action="store_true", dest="convert")

    (options, args) = parser.parse_args()
    
    # Print help if no argument given
    if util.printHelp(options):
        parser.print_help()
        sys.exit()
    
    # Get and check input arguments
    if options.list:
        # Read organism common name and related locus tag
        list_file = options.list
        util.checkFile(list_file)
        for line in open(list_file, "r"):
            if line[0] == '!':
                continue
            if not line.count('||') == 6:
                continue
            # ! common_name||organim_name||strain||locus_tag||genome_project_id||coverage||source
            line = line.strip()
            values = line.split('||')
            common_name=values[0]
            locus_tag=values[3]

            embl_file = "../IMG/%s.4dep.embl" % common_name
            util.checkFile(embl_file)
            tbl_file = "%s.tbl" % common_name
            log.info("Convert file %s into %s" % (embl_file, tbl_file))
            if options.convert:
                try:
                    doConvert(embl_file, tbl_file, locus_tag)
                except Exception, e:
                    log.error("Converting %s" % embl_file)
                    log.error(traceback.extract_stack())
                    log.error(e)
    if options.convert:
        doValidate()
        doClean()
                    
    if not options.convert:
        log.info("To perform the action, please use --convert")
            
### ---------------------------------------------------------------------------
if __name__ == '__main__':
    main()