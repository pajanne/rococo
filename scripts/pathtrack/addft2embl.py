'''
Created on Oct 22, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import os, sys
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import Reference, SeqFeature, FeatureLocation, ExactPosition, BeforePosition, AfterPosition
from optparse import OptionParser

### ---------------------------------------------------------------------------
features = []

### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------
def populateFeatures(tab_file):
    # read tab FT file
    """
FT  rRNA             1436971..1437083
FT                   /note="5s_rRNA"
FT                   /method="RNAmmer-1.2"

FT  tRNA             1614340..1614416
FT                   /note="tRNA-Arg(TCT) Cove Score 86.86"
FT                   /colour=4
FT                   /method="tRNAscan-SE"

FT   misc_feature    695000..700000
FT                   /colour=255 244 244
FT                   /algorithm="alien_hunter"
FT                   /note="threshold: 12.888"
FT                   /score=14.794
"""

    f_input = open (tab_file, 'r')
    ft_count = 0
    for line in f_input:
        values = line.strip().split()
        print values
        if values[0] == 'FT':
            if len(values) == 3:
                ft_type = values[1]
                print values[2]
                location = values[2].split('..')
                start = location[0]
                end = location[1]
                strand = 1
                if start[0:10] == 'complement':
                    start = start.split('(')[1]
                    end = end[:-1]
                    strand = -1
                if start[0] == '>':
                    start_position = AfterPosition(int(start[1:])-1)
                elif start[0] == '<':
                    start_position = BeforePosition(int(start[1:])-1)
                else:
                    start_position = ExactPosition(int(start)-1)
                if end[0] == '>':
                    end_position = AfterPosition(int(end[1:]))
                elif end[0] == '<':
                    end_position = BeforePosition(int(end[1:]))
                else:
                    end_position = ExactPosition(int(end))
                ft_count += 1
                feature = SeqFeature(FeatureLocation(start_position, end_position), strand=strand, type=ft_type)
                features.append(feature)
            if len(values) == 2:
                if values[1].startswith('/note'):
                    note = values[1].split('"')[1]
                    feature.qualifiers['note'] = [note]
                if values[1].startswith('/method'):
                    method = values[1].split('"')[1]
                    feature.qualifiers['method'] = [method]
                if values[1].startswith('/algorithm'):
                    method = values[1].split('"')[1]
                    feature.qualifiers['method'] = [method]
    f_input.close()
    print "INFO: %s features added from %s" % (ft_count, method)


### ---------------------------------------------------------------------------
def writeEmblFile(record, embl_file):
    features.extend(record.features)
    features.sort(compareFeatures)
    record.features = features
    print record
    SeqIO.write([record], open(embl_file, "w"), "embl")


### ---------------------------------------------------------------------------
def compareFeatures(fx, fy):
    if int(fx.location.start.position) > int(fy.location.start.position):
        return 1
    elif int(fx.location.start.position) == int(fy.location.start.position):
        return 0
    else: 
        return -1


### ---------------------------------------------------------------------------
def doRun():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", metavar="FILE", help="EMBL sequence input FILE", action="store", type="string", dest="in_embl")
    parser.add_option("-t", metavar="FILE", help="tab FILE", action="store", type="string", dest="tab")
    parser.add_option("-o", metavar="FILE", help="Merged EMBL output FILE", action="store", type="string", dest="out_embl")
    (options, args) = parser.parse_args()

    if not (options.in_embl and options.tab and options.out_embl):
        parser.print_help()
        sys.exit()

    # Create record from EMBL file
    if os.path.exists(options.in_embl):
        record = SeqIO.read(open(options.in_embl), "embl")
    # Read feature table
    if os.path.exists(options.tab):
        populateFeatures(options.tab)
    # Write EMBL file
    writeEmblFile(record, options.out_embl)
    
        
if __name__ == '__main__':
    doRun()

