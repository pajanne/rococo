'''
Created on Feb 24, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import os, sys
from optparse import OptionParser
from genepy import util
from Bio import SeqIO


### ---------------------------------------------------------------------------
def doConvert(gbk_file, tbl_file, locus_tag):
    record = SeqIO.read(open(gbk_file), "embl")
    table = open(tbl_file, 'w')
    table.write('>Feature %s\n' % locus_tag)
    # remove features on partial genes and features with illegal start codons
    features_to_remove = []
    for feature in record.features:
        for qualifier in feature.qualifiers:
            if 'translation' in qualifier:
                if 'X' in feature.qualifiers[qualifier][0]:
                    features_to_remove.append(getLocation(feature))
                #if not 'M' == feature.qualifiers[qualifier][0][1]:
                #    features_to_remove.append(getLocation(feature))
    
    print "Number of features removed %s" % len(features_to_remove)
    for feature in record.features:
        if not ('source' in feature.type):
            if getLocation(feature) not in features_to_remove:
                table.write('%s\t%s\n' % (getLocation(feature), feature.type))
                for qualifier in feature.qualifiers:
                    if not ('translation' in qualifier or 'note' in qualifier):
                        if 'locus_tag' in qualifier:
                            table.write('\t\t%s\t%s\n' % (qualifier, getLocusTag(feature.qualifiers[qualifier][0], locus_tag)))
                        elif 'tRNA' in feature.type and 'product' in qualifier:
                            continue
                        else:
                            table.write('\t\t%s\t%s\n' % (qualifier, feature.qualifiers[qualifier][0]))
                table.write('\t\tinference\tab initio prediction:IMG/ER\n')
    table.close()

### ---------------------------------------------------------------------------
def getLocusTag(value, locus_tag):
    return '%s_%5s' % (locus_tag, value[len(value)-5:len(value)])
        
### ---------------------------------------------------------------------------
def getLocation(feature):
    if feature.strand == -1:
        return "%s\t%s" % (feature.location.end.position, int(feature.location.start.position) + 1)
    else:
        return "%s\t%s" % (int(feature.location.start.position) + 1, feature.location.end.position)
    
### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", metavar="FILE", help="FILE containing the list of all organism common names and its associated locus tag", action="store", type="string", dest="list")
    parser.add_option("--convert", help="Do convert genbank file into embl", action="store_true", dest="convert")

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
            if line.count('||') < 1:
                continue
            # ! common_name||organim_name||strain||locus_tag||fasta_file
            line = line.strip()
            values = line.split('||')
            common_name=values[0]
            locus_tag=values[3]

            gbk_file = "%s.img.embl" % common_name
            util.checkFile(gbk_file)
            tbl_file = "%s.tbl" % common_name
            print "Convert file %s into %s" % (gbk_file, tbl_file)
            if options.convert:
                try:
                    doConvert(gbk_file, tbl_file, locus_tag)
                except Exception, e:
                    print "ERROR to convert %s" % gbk_file
                    print e
                    
    if not options.convert:
        print "To perform the action, please use --convert"
            
### ---------------------------------------------------------------------------
if __name__ == '__main__':
    main()