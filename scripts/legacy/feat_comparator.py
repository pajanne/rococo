'''
Created on Feb 19, 2010
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
    stat_dict = {}
    
    stat_dict['Seq length'] = len(record.seq)
    
    for feature in record.features:
        if not feature.type in feature_types:
            feature_types.append(feature.type)
    #print feature_types
    
    total = 0
    for type in feature_types:
        count = 0
        if type == 'source':
            continue
        for feature in record.features:
            if feature.type == type:
                count = count + 1
        total = total + count
        stat_dict[type] = count
    stat_dict['Total #features'] = total
    return stat_dict

### ---------------------------------------------------------------------------
def getExactMatch(record_a, record_b):
    exact_match_dict = {}
    count_non_exact_match = 0
    for feature_a in record_a.features:
        exact_match = False
        for feature_b in record_b.features:
            if feature_a.location.start == feature_b.location.start and feature_a.location.end == feature_b.location.end and feature_a.strand == feature_b.strand:
                exact_match = True
                key = getLocation(feature_a)
                exact_match_dict[key] = "%5s" % feature_a.type
                #print "A location is %s, B location is %s" % (feature_a.location, feature_b.location)
            if exact_match:
                break
        if not exact_match:
            count_non_exact_match = count_non_exact_match + 1
            #print "A feature type %s location %s not found in B" % (feature_a.type, feature_a.location)
    #print "%s features not found" % (count_non_exact_match)
    return exact_match_dict
    
### ---------------------------------------------------------------------------
def getExactMatch3(dict, record_c):
    exact_match_dict = {}
    count_non_exact_match = 0
    for feature_c in record_c.features:
        value = getLocation(feature_c)
        #print value
        if value in dict:
            exact_match_dict[value] = "%5s" % feature_c.type
        else:
            count_non_exact_match = count_non_exact_match + 1
            #print "A feature type %s location %s not found in B" % (feature_a.type, feature_a.location)
    #print "%s features not found" % (count_non_exact_match)
    return exact_match_dict

### ---------------------------------------------------------------------------
def getNoMatch(record_a, abc_matches, ab_matches, ac_matches):
    no_match_dict = {}
    for feature_a in record_a.features:
        key = getLocation(feature_a)
        if (key not in abc_matches) and (key not in ab_matches) and (key not in ac_matches):
            no_match_dict[key] = "%5s" % feature_a.type
    return no_match_dict

### ---------------------------------------------------------------------------
def getLocation(feature):
    if feature.strand == -1:
        return "complement(%s..%s)" % (feature.location.start, feature.location.end)
    else:
        return "(%s..%s)" % (feature.location.start, feature.location.end)

### ---------------------------------------------------------------------------
def doCompare(common_name, gendb_file, rast_file, img_file):
    gendb_record = SeqIO.read(open(gendb_file), "embl")
    rast_record = SeqIO.read(open(rast_file), "embl")
    img_record = SeqIO.read(open(img_file), "embl")

    print "Analysis of GenDB (A) features from %s" % gendb_file
    a_stat = stat(gendb_record)
    print "Analysis of RAST (B) features from %s" % rast_file
    b_stat = stat(rast_record)
    if not len(gendb_record.seq) == len(rast_record.seq):
        print "Sequences do not have the same length"
        
    print "Analysis of IMG (C) features from %s" % img_file
    c_stat = stat(img_record)
    ab_matches = getExactMatch(gendb_record, rast_record)
    abc_matches = getExactMatch3(ab_matches, img_record)
    # remove abc_matches in ab_matches
    for abc_match in abc_matches:
        if abc_match in ab_matches:
            del ab_matches[abc_match]
    bc_matches = getExactMatch(rast_record, img_record)
    # remove abc_matches in bc_matches
    for abc_match in abc_matches:
        if abc_match in bc_matches:
            del bc_matches[abc_match]
    ac_matches = getExactMatch(gendb_record, img_record)
    # remove abc_matches in ac_matches
    for abc_match in abc_matches:
        if abc_match in ac_matches:
            del ac_matches[abc_match]
    # get no match
    a_only = getNoMatch(gendb_record, abc_matches, ab_matches, ac_matches)
    b_only = getNoMatch(rast_record, abc_matches, ab_matches, bc_matches)
    c_only = getNoMatch(img_record, abc_matches, bc_matches, ac_matches)
    
    # write report
    res_file = open("%s.compare" % common_name, 'w')
    res_file.write("------------------------------\n")
    res_file.write("ABC: %s\n" % len(abc_matches))
    res_file.write("------------------------------\n")
    for match in abc_matches:
        res_file.write("     %s           %s\n" % (abc_matches[match],  match))
    res_file.write("------------------------------\n")
    res_file.write("AB: %s\n" % len(ab_matches))
    res_file.write("------------------------------\n")
    for match in ab_matches:
        res_file.write("     %s           %s\n" % (ab_matches[match],  match))
    res_file.write("------------------------------\n")
    res_file.write("BC: %s\n" % len(bc_matches))
    res_file.write("------------------------------\n")
    for match in bc_matches:
        res_file.write("     %s           %s\n" % (bc_matches[match],  match))
    res_file.write("------------------------------\n")
    res_file.write("AC: %s\n" % len(ac_matches))
    res_file.write("------------------------------\n")
    for match in ac_matches:
        res_file.write("     %s           %s\n" % (ac_matches[match],  match))
    res_file.write("------------------------------\n")
    res_file.write("A: %s\n" % len(a_only))
    res_file.write("------------------------------\n")
    for match in a_only:
        res_file.write("     %s           %s\n" % (a_only[match],  match))
    res_file.write("------------------------------\n")
    res_file.write("B: %s\n" % len(b_only))
    res_file.write("------------------------------\n")
    for match in b_only:
        res_file.write("     %s           %s\n" % (b_only[match],  match))
    res_file.write("------------------------------\n")
    res_file.write("C: %s\n" % len(c_only))
    res_file.write("------------------------------\n")
    for match in c_only:
        res_file.write("     %s           %s\n" % (c_only[match],  match))
    res_file.write("------------------------------\n")
    res_file.write("Summary:\n")
    res_file.write("------------------------------\n")
    res_file.write("A: GenDB; B: RAST; C: IMG\n")
    res_file.write("------------------------------\n")
    for key in a_stat:
        res_file.write("%20s%10s\n" % (key, a_stat[key]))
    res_file.write("------------------------------\n")
    for key in b_stat:
        res_file.write("%20s%10s\n" % (key, b_stat[key]))
    res_file.write("------------------------------\n")
    for key in c_stat:
        res_file.write("%20s%10s\n" % (key, c_stat[key]))
    res_file.write("------------------------------\n")
    res_file.write("   ABC    AB    BC    AC     A     B     C\n")
    res_file.write("%6s%6s%6s%6s%6s%6s%6s\n" % (len(abc_matches), len(ab_matches), len(bc_matches), len(ac_matches), len(a_only), len(b_only), len(c_only)))
    res_file.write("------------------------------\n")

### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", "--list", metavar="FILE", help="FILE containing the list of all organism common names to compare", action="store", type="string", dest="list")
    
    (options, args) = parser.parse_args()

    # Print help if no argument given
    if util.printHelp(options):
        parser.print_help()
        sys.exit()
    
    if options.list:
        # Read organism common name and related fasta sequence file
        list_file = options.list
        util.checkFile(list_file)
        for line in open(list_file, "r"):
            if line[0] == '!':
                continue
            # ! common_name
            common_name = line.strip()
                        
            gendb_file = "GenDB/%s.gendb.embl" % common_name
            rast_file = "RAST/%s.rast.embl" % common_name
            img_file = "IMG/%s.img.embl" % common_name
            if not os.path.exists(gendb_file) or not os.path.exists(rast_file) or not os.path.exists(img_file):
                print "No three results for %s" % common_name
                continue
            
            print "Processing %s" % common_name
            doCompare(common_name, gendb_file, rast_file, img_file)

if __name__ == '__main__':
    main()
