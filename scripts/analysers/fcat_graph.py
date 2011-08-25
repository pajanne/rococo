'''
Created on Aug 24, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.

'''

import sys, os
from rpy import *

### ---------------------------------------------------------------------------
### constants
### ---------------------------------------------------------------------------
categories = {}
organisms = {}
goslim = {}
cat_set = set()

### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------

### ---------------------------------------------------------------------------
def doGraph():
    for organism_name in organisms.keys():
        for category in organisms[organism_name].keys():
            if category not in cat_set:
                cat_set.add(category)

    for cat in cat_set:
        frequencies = {}
        for organism_name in organisms.keys():
            if cat in organisms[organism_name].keys():
                frequencies[organism_name] = organisms[organism_name][cat]
            else:
                frequencies[organism_name] = 0
        categories[cat] = frequencies

    r.pdf('fcat_graph.pdf')
    r.par(mfrow=(2,2))
    for cat in categories.keys():
        if cat in goslim.keys():
            main_title = "%s - %s" % (cat, goslim[cat])
        else:
            "WARNING: GO slim %s does not exist!" % cat
            main_title = "%s" % (cat)
        coords = r.barplot(categories[cat].values(), main=main_title, ylab="Organisms", xlab="Frequency", horiz="TRUE")
        r.text(0, coords, categories[cat].keys(), xpd=True, cex=0.5)
    r.dev_off()


### ---------------------------------------------------------------------------
def doAnalyse(file, common_name):
    """
    pfamscan stats
GO:0003676: 7
GO:0043226: 2
GO:0005488: 61
UNK_PFAM2GO: 702
    """
    frequencies = {}
    # read pfam stats data
    f_input = open (file, 'r')
    for line in f_input:
        line = line.strip()
        values = line.split()
        if not len(values) == 2:
            continue
        key = values[0][:-1]
        value = int(values[1])
        frequencies[key] = value
    organisms[common_name] = frequencies

        
### ---------------------------------------------------------------------------
def doRun(list):
    print "Reading input file %s" % list
    processed_lines = 0
    for line in open(list, "r"):
        if line[0] == '!':
            continue
        # common_name||pfamstats_file
        line = line.strip()
        values = line.split('||')
        if not len(values) == 2:
            continue
        processed_lines += 1
        common_name = values[0]
        pfamstats_file = values[1]

        # Check if Pfam stats file exists
        if not os.path.exists(pfamstats_file):
            print "ERROR: file %s does not exist" % pfamstats_file
            sys.exit(1)
            
        print "Analysing file %s ..." % (pfamstats_file)
        doAnalyse(pfamstats_file, common_name)

    doGraph()

    print "%s line(s) processed from %s" % (processed_lines, list)


### ---------------------------------------------------------------------------
def doGoSlim(file):
    """
[Term]
id: GO:0000003
name: reproduction

    """
    print "Reading GO slim terms from input file %s" % file
    is_term = False
    for line in open(file, "r"):
        line = line.strip()
        values = line.split(' ')
        if line == '[Term]':
            is_term = True
        if is_term and values[0] == 'id:':
            go_id = values[1]
        if is_term and values[0] == 'name:':
            go_name = " ".join(values[1:])
            is_term = False
            goslim[go_id] = go_name

        
### ---------------------------------------------------------------------------
def main():
    from optparse import OptionParser

    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", "--list", metavar="FILE", help="FILE containing the list of Pfam stats files to analyse", action="store", type="string", dest="list")
    parser.add_option("-g", "--goslim", metavar="FILE", help="FILE containing GO slim terms", action="store", type="string", dest="goslim")

    (options, args) = parser.parse_args()
    
    # Print help if no argument given
    if not (options.list and options.goslim):
        parser.print_help()
        sys.exit()
    
    # Check input file list exists
    list = options.list
    if not os.path.exists(list):
        print "ERROR: file %s does not exist" % list
        sys.exit(1)

    goslim = options.goslim
    if not os.path.exists(goslim):
        print "ERROR: file %s does not exist" % go_slim
        sys.exit(1)

    doGoSlim(goslim)
    doRun(list)
            
### ---------------------------------------------------------------------------
if __name__ == '__main__':
    main()
