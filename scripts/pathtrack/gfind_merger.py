'''
Created on Sep 23, 2010
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
def getRecordFromSeq(fasta_file, organism_name):
    # read FASTA sequence file 
    records = SeqIO.parse(open(fasta_file), "fasta")
    new_seq = ''
    gap_count = 0
    for record in records:
        # Add 50 N between scaffolds, but not at the beginning
        if new_seq == '':
            new_seq = "%sNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" % (record.seq)
        else:
            new_seq = "%sNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN%s" % (new_seq, record.seq)
        # Add scaffold gaps
        start_N = len(new_seq)
        end_N = start_N + 50
        gap_feature = SeqFeature(FeatureLocation(start_N, end_N), strand=1, type="gap")
        gap_feature.qualifiers['estimated_length'] = ['50']
        gap_feature.qualifiers['note'] = ['scaffold gap']
        features.append(gap_feature)
        gap_count += 1
        
    new_record = SeqRecord(Seq(new_seq, generic_dna), id=organism_name)
    new_record.seq.alphabet = generic_dna
    print new_record
    print "INFO: %s scaffold gap features added" % gap_count
    return new_record


### ---------------------------------------------------------------------------
def populateFeatures(tab_file, method):
    # read tab FT file
    """
FT   CDS             complement(<1..174)
FT                   /colour=4
FT                   /method="PRODIGAL"
"""

    f_input = open (tab_file, 'r')
    cds_count = 0
    for line in f_input:
        values = line.strip().split()
        if len(values) == 3:
            if values[0] == 'FT' and values[1] == 'CDS':
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
                cds_count += 1
                cds_feature = SeqFeature(FeatureLocation(start_position, end_position), strand=strand, type="CDS")
                cds_feature.qualifiers['method'] = [method]
                features.append(cds_feature)
    f_input.close()
    print "INFO: %s CDS features added from %s" % (cds_count, method)

    
### ---------------------------------------------------------------------------
def addGaps(record):
    seq = record.seq
    in_N = False
    # TODO - Cope with a sequence which ends with N
    if seq[-1] == "N":
        print "WARNING: sequence ends with N"
    else:
        for i in range(len(seq)):
            if seq[i] == 'N' and not in_N:
                start_N = i
                in_N = True
            if in_N and not seq[i+1] == 'N':
                end_N = i + 1
                length = end_N - start_N
                assert length > 0
                assert str(seq[start_N:end_N]) == "N"*length
                # do not create FT for 1bp gap
                if length > 1:
                    gap_feature = SeqFeature(FeatureLocation(start_N, end_N), strand=1, type="gap")
                    gap_feature.qualifiers['estimated_length'] = [length]
                    gap_feature.qualifiers['note'] = ['contig gap']
                    features.append(gap_feature)
                in_N = False


### ---------------------------------------------------------------------------
def mergeFeatures(record):
    # Total number of features
    print "INFO: %s features in total" % len(features)
    # Remove features on N regions
    n_count = 0
    for feature in features:
        if feature.type == 'CDS':
            if record.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_end].count('N') > 1:
                n_count += 1
                features.remove(feature)
    print "INFO: %s CDS features removed located on N regions" % n_count
    # Group duplicate CDS features into one & remove duplicate scaffold gap features
    duplicate_count = 0
    i = 0
    while i < len(features):
        fx = features[i]
        i += 1
        j = i + 1
        while j < len(features):
            fy = features[j]
            j += 1
            if fx.type == 'CDS' and fy.type == 'CDS':
                if fx.location.nofuzzy_start == fy.location.nofuzzy_start and fx.location.nofuzzy_end == fy.location.nofuzzy_end:
                    duplicate_count += 1
                    features.remove(fy)
                    features[i].qualifiers['method'].append(fy.qualifiers['method'][0])
                    features[i].qualifiers['note'] = ['duplicate prediction']
            if fx.type == 'gap' and fy.type == 'gap':
                if fx.location.nofuzzy_start == fy.location.nofuzzy_start and fx.location.nofuzzy_end == fy.location.nofuzzy_end:
                    if fx.qualifiers['note'][0] == 'scaffold gap':
                        features.remove(fy)
                    else:
                        features.remove(fx)
    print "INFO: %s CDS features identical" % duplicate_count
    # Remove features where correlation score is below 50
    print "INFO: %s features" % len(features)


### ---------------------------------------------------------------------------
"""
Correlation score

 /**
  *  Return a String containing the correlation scores.
  **/
 protected String getScoresString(final Feature feature) {
   final int base_total = feature.getTranslationBases().length();

   final int c_total = feature.getBaseCount(Bases.getIndexOfBase('c'));
   final int g_total = feature.getBaseCount(Bases.getIndexOfBase('g'));

   final int g1_count = feature.getPositionalBaseCount(0, Bases.getIndexOfBase('g'));

   final int c3_count = feature.getPositionalBaseCount(2, Bases.getIndexOfBase('c'));
   final int g3_count = feature.getPositionalBaseCount(2, Bases.getIndexOfBase('g'));

   final double c3_score = 100.0 * (3 * c3_count - c_total) / c_total;
   final double g1_score = 100.0 * (3 * g1_count - g_total) / g_total;
   final double g3_score = 100.0 * (3 * g3_count - g_total) / g_total;

   final double cor1_2_score = feature.get12CorrelationScore();

   final NumberFormat number_format = NumberFormat.getNumberInstance();

   number_format.setMaximumFractionDigits(1);
   number_format.setMinimumFractionDigits(1);

   final String cor1_2_score_string = number_format.format(cor1_2_score);
   final String c3_score_string;
   final String g1_score_string;
   final String g3_score_string;


   if(c_total == 0)
     c3_score_string = "ALL";
   else 
     c3_score_string = number_format.format(c3_score);

   if(g_total == 0)
     g1_score_string = "ALL";
   else 
     g1_score_string = number_format.format(g1_score);

   if(g_total == 0)
     g3_score_string = "ALL";
   else
     g3_score_string = number_format.format(g3_score);

   String codon_usage_score_string = "";

   final CodonUsageAlgorithm codon_usage_alg = getBasePlotGroup().getCodonUsageAlgorithm();

   if(codon_usage_alg != null) {
     number_format.setMaximumFractionDigits(3);
     number_format.setMinimumFractionDigits(3);

     codon_usage_score_string = number_format.format(codon_usage_alg.getFeatureScore(feature)) + " ";
   }

   return codon_usage_score_string + padRightWithSpaces(cor1_2_score_string, 5) + " " +
          padRightWithSpaces(c3_score_string, 5) + " " +
          padRightWithSpaces(g1_score_string, 5) + " " +
          padRightWithSpaces(g3_score_string, 5);
 }

 /**
  *  Calculate the codon usage score for the given Feature.
  **/
 public float getFeatureScore (final Feature feature) {
   final String sequence = feature.getTranslationBases ();

   float total = 0F;

   for (int i = 0 ; i < sequence.length () ; i += 3) {

     final char base1 = sequence.charAt (i);
     final char base2 = sequence.charAt (i + 1);
     final char base3 = sequence.charAt (i + 2);

     final float this_weight =
       usage_data.getCodonValue (base1, base2, base3);

     total += Math.log (this_weight);
   }

   final int codon_count = sequence.length () / 3;

   return (float) Math.exp (total / codon_count);
 }

"""

### ---------------------------------------------------------------------------
def addLocusTag(prefix):
    # Sort the features
    features.sort(compareFeatures)
    # Add ids to all CDS features - prefix_000000 with 10 increment
    i = 0
    for feature in features:
        if feature.type == 'CDS':
            i += 1
            locus_tag = "%s_%06d" % (prefix, i*10)
            feature.qualifiers['locus_tag'] = [locus_tag]


### ---------------------------------------------------------------------------
def writeEmblFile(record, embl_file):
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
    parser.add_option("-s", metavar="FILE", help="FASTA sequence FILE", action="store", type="string", dest="seq")
    parser.add_option("-p", metavar="FILE", help="Prodigal EMBL tab FILE", action="store", type="string", dest="prodigal")
    parser.add_option("-g", metavar="FILE", help="Glimmer3 EMBL tab FILE", action="store", type="string", dest="glimmer")
    parser.add_option("-o", metavar="FILE", help="Merged EMBL output FILE", action="store", type="string", dest="embl")
    parser.add_option("-n", metavar="NAME", help="Organism NAME", action="store", type="string", dest="name")
    parser.add_option("-l", metavar="PREFIX", help="Locus tag PREFIX", action="store", type="string", dest="prefix")
    (options, args) = parser.parse_args()

    if not (options.seq and options.prodigal and options.glimmer and options.embl and options.name and options.prefix):
        parser.print_help()
        sys.exit()

    # Create record from FASTA sequence file
    if os.path.exists(options.seq):
        record = getRecordFromSeq(options.seq, options.name)
    # Read CDSs feature tables
    if os.path.exists(options.prodigal):
        populateFeatures(options.prodigal, 'PRODIGAL')
    if os.path.exists(options.glimmer):
        populateFeatures(options.glimmer, 'GLIMMER3')
    # Add gap features
    addGaps(record)
    # Modify features - remove duplicates
    mergeFeatures(record)
    # Add locus_tag
    addLocusTag(options.prefix)
    # Write EMBL file
    writeEmblFile(record, options.embl)
    
        
if __name__ == '__main__':
    doRun()

