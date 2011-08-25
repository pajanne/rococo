'''
Created on Dec 15, 2009
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

from Bio import SeqIO
import sys, os
from optparse import OptionParser
import util

### ---------------------------------------------------------------------------
### Main methods
### ---------------------------------------------------------------------------
def stat(record):
    feature_types = []

    print "Sequence length: %s" % len(record.seq)
    
    for feature in record.features:
        if not feature.type in feature_types:
            feature_types.append(feature.type)
    #print feature_types
    
    total_nb_features = 0
    for type in feature_types:
        count = 0
        if type == 'source':
            continue
        for feature in record.features:
            if feature.type == type:
                count = count + 1
        total_nb_features = total_nb_features + count
        print "Number of %s: %i" % (type, count)
    print "Total number of features: %i" % total_nb_features

### ---------------------------------------------------------------------------
def compare(record_a, record_b):
    count_exact_match = 0
    for feature_a in record_a.features:
        exact_match = False
        for feature_b in record_b.features:
            if feature_a.location.start == feature_b.location.start and feature_a.location.end == feature_b.location.end:
                exact_match = True
                #print "A location is %s, B location is %s" % (feature_a.location, feature_b.location)
            if exact_match:
                break
        if not exact_match:
            count_exact_match = count_exact_match + 1
            print "A feature type %s location %s not found in B" % (feature_a.type, feature_a.location)
    print "%s features not found" % (count_exact_match)
    
### ---------------------------------------------------------------------------
def transfer(record_a, record_b):
    count_exact_match = 0
    count_product_transfer = 0
    for feature_a in record_a.features:
        exact_match = False
        for feature_b in record_b.features:
            if feature_a.location.start == feature_b.location.start and feature_a.location.end == feature_b.location.end:
                exact_match = True
                if 'product' not in feature_b.qualifiers.keys():
                    # add it into A with a note 
                    if 'product' in feature_a.qualifiers.keys():
                        feature_b.qualifiers['product'] = feature_a.qualifiers['product']
                        feature_b.qualifiers['note'] = 'product obtained from automatic annotation'
                        count_product_transfer = count_product_transfer + 1
                        print feature_a.qualifiers['product']
                #print "A location is %s, B location is %s" % (feature_a.location, feature_b.location)
            if exact_match:
                break
        if not exact_match:
            count_exact_match = count_exact_match + 1
            print "A feature type %s location %s not found in B" % (feature_a.type, feature_a.location)
    print "%s features not found" % (count_exact_match)
    print "%s products has been added" % (count_product_transfer)
    return record_b
    
### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-a", metavar="EMBL-a", help="First EMBL file", action="store", type="string", dest="first_embl")
    parser.add_option("-b", metavar="EMBL-b", help="Second EMBL file to compare", action="store", dest="second_embl")
    parser.add_option("--merge", help="To transfer /product of identical annotations into a merged file", action="store_true", dest="merge")
    
    (options, args) = parser.parse_args()

    # Print help if no argument given
    if util.printHelp(options):
        parser.print_help()
        sys.exit()
    
    first_record = SeqIO.read(open(options.first_embl), "embl")
    second_record = SeqIO.read(open(options.second_embl), "embl")

    print "Analysis of EMBL features A from %s" % options.first_embl
    print "Analysis of EMBL features B from %s" % options.second_embl

    stat(first_record)
    
    if options.merge:
        merged_record = transfer(first_record, second_record)
        # Write out genbank file
        SeqIO.write([merged_record], open("merged.embl", "w"), "embl")
    

if __name__ == '__main__':
    main()
