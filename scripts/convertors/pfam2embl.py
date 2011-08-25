'''
Created on Aug 18, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import os, sys, copy
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import Reference, SeqFeature, FeatureLocation, ExactPosition, BeforePosition, AfterPosition
from optparse import OptionParser

### ---------------------------------------------------------------------------
go_regex = re.compile('GO:[0-9]{7}')
PFAM2GO = {}
GO2SLIM = {}
SEQIDS = {}
GENES = set()

### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------
"""
# <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan>

Alistipes_shahii_WAL8301_481_3516_CDS_1           172    313    172    314 PF07703.7   A2M_N_2           Family     1   135   136     77.7   5.4e-22   1 No_clan
  
Alistipes_shahii_WAL8301_481_3516_CDS_1           389    458    380    468 PF00207.15  A2M               Family    13    84    92     28.6   6.1e-07   1 CL0011 
  
Alistipes_shahii_WAL8301_481_3516_CDS_1           574    598    573    600 PF10569.2   Thiol-ester_cl    Domain     2    26    31     23.2   2.2e-05   1 No_clan
  

"""
def convertToEmbl(file, embl_file, tab_file, stats_file):

    new_features = []

    # read EMBL feature table
    record = SeqIO.read(open(embl_file), "embl")
    features = {}
    for i in range(len(record.features)):
        feature = record.features[i]
        if 'locus_tag' in feature.qualifiers.keys():
            features[feature.qualifiers['locus_tag'][0]] = feature
        if 'label' in feature.qualifiers.keys():
            features[feature.qualifiers['label'][0]] = feature

    # read pfam raw data
    f_input = open (file, 'r')
    cds_cat = {}
    for line in f_input:
        line = line.strip()
        values = re.split('\s+', line)
        if len(values) == 0:
            continue
        if values[0] == '#':
            continue
        if not len(values) == 15:
            continue
        cds_id = values[0] # XXX_599515_603324_CDS_1
        al_start = int(values[1]) # 172
        al_end = int(values[2]) # 313
        env_start = int(values[3]) # 172
        env_end = int(values[4]) # 314
        hmm_acc = values[5] # PF07703.7
        pfam_acc = hmm_acc.split('.')[0]
        hmm_name = values[6] # A2M_N_2
        hmm_type = values[7] # Family
        hmm_start = int(values[8]) # 1
        hmm_end = int(values[9]) # 135
        hmm_length = int(values[10]) # 136
        bit_score = values[11] # 77.7
        evalue = values[12] # 5.4e-22
        significance = values[13]
        clan = values[14]

        # find feature associated with this pfam result line
        feature = features[SEQIDS[cds_id]]
        
        # add new misc_feature
        new_feature = copy.copy(feature)
        new_feature.type = 'misc_feature'
        new_feature.qualifiers = {}
        new_feature.qualifiers['id'] = SEQIDS[cds_id]
        new_feature.qualifiers['note'] = "Pfam domain hit to %s (%s in clan %s), at HMM positions [%s-%s], bit score: %s, e-value:%s" % (hmm_acc, hmm_name, clan, hmm_start, hmm_end, bit_score, evalue)
        new_features.append(new_feature)
        
        # update EMBL feature table with Pfam accession code
        if feature.qualifiers.has_key('db_xref'):
            feature.qualifiers['db_xref'].append('PFAM:%s' % hmm_acc)
        else:
            feature.qualifiers['db_xref'] = ['PFAM:%s' % hmm_acc]
        
        # build a set of genes
        GENES.add(cds_id)

        # assign Pfam accession to GO then to GOslim
        if PFAM2GO.has_key(pfam_acc):
            go_id = PFAM2GO[pfam_acc]
            if GO2SLIM.has_key(go_id):
                goslim = GO2SLIM[go_id]
            else:
                goslim = "UNK_GO2SLIM"
        else:
            goslim = "UNK_PFAM2GO"
        if cds_cat.has_key(goslim):
            cds_cat[goslim].append(cds_id)
        else:
            cds_cat[goslim] = [cds_id]

        # add modified feature
        features[SEQIDS[cds_id]] = feature

    # write new EMBL feature table
    new_features.extend(features.values())
    new_features.sort(feature_compare)
    record.features = new_features
    SeqIO.write([record], open(tab_file, "w"), "embl")

    # write functional category stats
    f_stats = open (stats_file, 'w')
    cds_num = 0
    for key in cds_cat.keys():
        cds_num += len(cds_cat[key])
        f_stats.write("%s: %s\n" % (key, len(cds_cat[key])))
    f_stats.write("%s\n" % cds_num)
    f_stats.write("%s\n" % len(GENES))
    f_stats.write("%s\n" % len(cds_cat.keys()))

    # close files
    f_input.close()
    f_stats.close()


### ---------------------------------------------------------------------------
def feature_compare(fx, fy):
    if int(fx.location.start.position) > int(fy.location.start.position):
        return 1
    elif int(fx.location.start.position) == int(fy.location.start.position):
        return 0
    else: 
        return -1


### ---------------------------------------------------------------------------
def readPfam2Go(map):
    f_map = open(map, 'r')
    for line in f_map:
        if line[0] == '!':
            continue
        line.strip()
        values = line.split('>')
        key = values[0].split()[0][5:]
        value = values[1].split(';')[1].strip()
        if go_regex.search(value):
            PFAM2GO[key] = value
    f_map.close()


### ---------------------------------------------------------------------------
def readGoMap(map):
    f_map = open(map, 'r')
    for line in f_map:
        line.strip()
        values = line.split('=>')
        key = values[0].strip()
        value = values[1].strip().split(' ')[0].strip()
        if go_regex.search(key) and go_regex.search(value):
            GO2SLIM[key] = value
    f_map.close()


### ---------------------------------------------------------------------------
def readSeqIds(map):
    f_map = open(map, 'r')
    for line in f_map:
        line.strip()
        values = line.split()
        key = values[0][1:]
        value = values[2].split('"')[1]
        SEQIDS[key] = value
    f_map.close()
        

### ---------------------------------------------------------------------------
def doRun():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", metavar="FILE", help="Pfam-scan raw input FILE", action="store", type="string", dest="input")
    parser.add_option("-e", metavar="FILE", help="EMBL input FILE", action="store", type="string", dest="embl")
    parser.add_option("-o", metavar="FILE", help="Pfam-scan output FILE converted into EMBL", action="store", type="string", dest="output")
    parser.add_option("-s", metavar="FILE", help="Pfam-scan output stats FILE", action="store", type="string", dest="stats")
    parser.add_option("-p", metavar="FILE", help="mapping FILE between Pfam and GO terms", action="store", type="string", dest="pfam2go")
    parser.add_option("-g", metavar="FILE", help="mapping FILE between GO and GOslim terms generated with map2slim (e.g. map2slim goslim_generic.obo gene_ontology_ext.obo -outmap map)", action="store", type="string", dest="go")
    parser.add_option("-m", metavar="FILE", help="mapping FILE between EMBOSS and EMBL seq ids from description line of fasta file generated by extractfeat (e.g. extractfeat -sequence Ashahii_WAL8301.4dep.embl -type CDS -featinname Yes -describe locus_tag -stdout Yes -auto Yes | transeq -filter Yes -stdout Yes -auto Yes | grep '>' > Ashahii_WAL8301.seqids)", action="store", type="string", dest="map")
    (options, args) = parser.parse_args()

    if not (options.input and options.embl and options.output and options.stats and options.pfam2go and options.go and options.map):
        parser.print_help()
        sys.exit()

    # Run the conversion 
    if os.path.exists(options.input):
        # Read mapping between Pfam and GO - create a dictionary of mapping PFAM2GO
        readPfam2Go(options.pfam2go)
        # Read mapping between GO and GOslim - create a dictionary of mapping GO2SLIM
        readGoMap(options.go)
        # Read mapping between EMBOSS ans EMBL seq ids - create a dictionary of mapping SEQIDS
        readSeqIds(options.map)
        # Convert output results into a feature table EMBL file and build up functional categories based on GOslim
        convertToEmbl(options.input, options.embl, options.output, options.stats)
        
if __name__ == '__main__':
    doRun()
