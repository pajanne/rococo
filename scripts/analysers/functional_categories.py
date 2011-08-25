'''
Created on Jul 19, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.

'''

import sys, os
from Bio import SeqIO
from Bio.SeqFeature import Reference, SeqFeature, FeatureLocation, ExactPosition, BeforePosition, AfterPosition
from datetime import date
import re
import traceback
import ftplib
import subprocess

### ---------------------------------------------------------------------------
### constants
### ---------------------------------------------------------------------------
categories = {}
organisms = []

### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------

### ---------------------------------------------------------------------------
def doAnalyse(embl_file):
    record = SeqIO.read(open(embl_file), "embl")
    organism = record.annotations["organism"]
    organisms.append(organism)
    
    for i in range(len(record.features)):
        feature = record.features[i]

        if feature.type == 'CDS':
            if 'product' in feature.qualifiers.keys():
                product_name = feature.qualifiers['product'][0]
                if not categories.has_key(product_name):
                    categories[product_name] = {organism: 1}
                else:
                    if not categories[product_name].has_key(organism):
                        categories[product_name] = {organism: 1}
                    else:
                        categories[product_name][organism] += 1

    #cat_nb = 0
    #for value, category in sorted([(value, category) for (category, value) in categories.items()]):
    #    if value > 1:
    #        print "%s\t%s" % (value, category)
    #        cat_nb += 1

    #print "Number of categories: %s" % len(categories)
    #print "Nb of categories > 1: %s" % cat_nb
        
### ---------------------------------------------------------------------------
def doRun(list):
    print "Reading input file %s" % list
    processed_lines = 0
    for line in open(list, "r"):
        if line[0] == '!':
            continue
        # embl_file
        processed_lines += 1 
        line = line.strip()
        values = line.split()
        embl_file=values[0]

        # Check if EMBL file exists
        if not os.path.exists(embl_file):
            print "ERROR: file %s does not exist" % embl_file
            sys.exit(1)
            
        print "Analysing file %s ..." % (embl_file)
        doAnalyse(embl_file)
            
    print "%s processed lines from input file" % processed_lines
    for key in sorted([key for key in categories.keys()]):
        value = ""
        for organism in organisms:
            if categories[key].has_key(organism):
                value = "%s\t%s" % (value, categories[key][organism])
            else:
                value = "%s\t0" % value

        print "%s\t%s" % (value, key)
    print organisms
    print len(categorie)
    

### ---------------------------------------------------------------------------
def main():
    from optparse import OptionParser

    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", "--list", metavar="FILE", help="FILE containing the list of EMBL files to analyse", action="store", type="string", dest="list")

    (options, args) = parser.parse_args()
    
    # Print help if no argument given
    if not (options.list):
        parser.print_help()
        sys.exit()
    
    # Check input file list exists
    list = options.list
    if not os.path.exists(list):
        print "ERROR: file %s does not exist" % list
        sys.exit(1)    

    doRun(list)
            
### ---------------------------------------------------------------------------
if __name__ == '__main__':
    main()
