'''
Created on Aug 8, 2010
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
GO2SLIM = {}
SEQIDS = {}
GENES = set()
COLOURS = {
    "ScanRegExp"  : "8",
    "Seg"         : "0",
    "ProfileScan" : "7",
    "TMHMM"       : "0",
    "Coil"        : "0",
    "HMMTigr"     : "5",
    "SignalPHMM"  : "0",
    "FPrintScan"  : "4",
    "superfamily" : "0",
    "BlastProDom" : "0",
    "HMMPfam"     : "9",
    "HMMSmart"    : "13",
    "HMMPanther"  : "0",
    "Gene3D"      : "0",
    "PatternScan" : "0",
    "HMMPIR"      : "0",
    "HAMAP"       : "0"
}

### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------
"""
XXX_599515_603324_CDS_1	877CD7BD496BE7E6	1270	HMMSmart	SM00481	no description	5	118	8.3e-22	T	02-Aug-2010	IPR003141	Polymerase/histidinol phosphatase, N-terminal	Molecular Function: DNA binding (GO:0003677), Molecular Function: DNA-directed DNA polymerase activity (GO:0003887), Biological Process: DNA replication (GO:0006260)
XXX_607511_609082_CDS_1	41B4CE154965B7EA	524	HMMSmart	SM00228	no description	86	157	0.0045	T	02-Aug-2010	IPR001478	PDZ/DHR/GLGF	Molecular Function: protein binding (GO:0005515)

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

    # read iprscan raw data
    f_input = open (file, 'r')
    cds_cat = {}
    for line in f_input:
        line = line.strip()
        values = line.split('\t')
        if len(values) == 0:
            continue
        if not len(values) == 14:
            continue
        cds_id = values[0] # XXX_599515_603324_CDS_1
        cds_length = int(values[2]) # 1270
        program = values[3] # HMMSmart
        program_acc = values[4] # SM00481
        short_description = values[5] # no description
        start = int(values[6]) # 5
        end = int(values[7]) # 118
        score = values[8] # 8.3e-22
        date = values[10] # 02-Aug-2010
        interpro_acc = values[11] # IPR003141
        description = values[12] # Polymerase/histidinol phosphatase, N-terminal
        go_description = values[13] # Molecular Function: DNA binding (GO:0003677), Molecular Function: DNA-directed DNA polymerase activity (GO:0003887), Biological Process: DNA replication (GO:0006260)

        # find feature associated with this iprscan result line
        feature = features[SEQIDS[cds_id]]
        
        # add new misc_feature
        new_feature = copy.copy(feature)
        new_feature.type = 'misc_feature'
        new_feature.qualifiers = {}
        new_feature.qualifiers['colour'] = COLOURS[program]
        new_feature.qualifiers['domain'] = "%s;%s;%s;%s;codon %s-%s" % (program, program_acc, short_description, score, start, end)
        new_feature.qualifiers['id'] = SEQIDS[cds_id]
        new_feature.qualifiers['label'] = program
        new_feature.qualifiers['note'] = "%s hit to %s (%s), score %s" % (program, program_acc, short_description, score)
        new_features.append(new_feature)
        
        # update EMBL feature table with InterPro ids
        feature.qualifiers['db_xref'] = ['InterPro:%s' % interpro_acc]
        
        # build a set of genes
        GENES.add(cds_id)

        # assign genes to categories using GOslim
        go_values = go_description.split(',')
        goslim = ''
        for go in go_values:
            go = go.strip()
            go_id = go_regex.search(go)
            go_key = go.split(':')[0]
            if go_key == 'Molecular Function':
                go_key = 'F'
            elif go_key == 'Cellular Component':
                go_key = 'C'
            else:
                go_key = 'P'
            if go_id:
                go_id = go_id.group(0)
                if GO2SLIM.has_key(go_id):
                    goslim = "%s:%s" % (go_key, GO2SLIM[go_id])
                else:
                    goslim = "%s:UNK" % go_key
                if cds_cat.has_key(goslim):
                    cds_cat[goslim].append(cds_id)
                else:
                    cds_cat[goslim] = [cds_id]
                # update EMBL feature table with GO ids
                feature.qualifiers['db_xref'].append(go_id)
        #f_stats.write("%s => %s\n" % (cds_id, goslim))

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
        if key[0] == 'F':
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
    parser.add_option("-i", metavar="FILE", help="iprscan raw input FILE", action="store", type="string", dest="input")
    parser.add_option("-e", metavar="FILE", help="EMBL input FILE", action="store", type="string", dest="embl")
    parser.add_option("-o", metavar="FILE", help="iprscan output FILE converted into EMBL feature tab", action="store", type="string", dest="output")
    parser.add_option("-s", metavar="FILE", help="iprscan output stats FILE", action="store", type="string", dest="stats")
    parser.add_option("-g", metavar="FILE", help="mapping FILE between GO and GOslim terms generated with map2slim (e.g. map2slim goslim_generic.obo gene_ontology_ext.obo -outmap map)", action="store", type="string", dest="go")
    parser.add_option("-m", metavar="FILE", help="mapping FILE between EMBOSS and EMBL seq ids from description line of fasta file generated by extractfeat (e.g. extractfeat -sequence Ashahii_WAL8301.4dep.embl -type CDS -featinname Yes -describe locus_tag -stdout Yes -auto Yes | transeq -filter Yes -stdout Yes -auto Yes | grep '>' > Ashahii_WAL8301.seqids)", action="store", type="string", dest="map")
    (options, args) = parser.parse_args()

    if not (options.input and options.embl and options.output and options.stats and options.go and options.map):
        parser.print_help()
        sys.exit()

    # Run the conversion 
    if os.path.exists(options.input):
        # Read mapping between GO and GOslim - create a dictionary of mapping GO2SLIM
        readGoMap(options.go)
        # Read mapping between EMBOSS ans EMBL seq ids - create a dictionary of mapping SEQIDS
        readSeqIds(options.map)
        # Convert output results into a feature table EMBL file and build up functional categories based on GOslim
        convertToEmbl(options.input, options.embl, options.output, options.stats)
        
if __name__ == '__main__':
    doRun()
