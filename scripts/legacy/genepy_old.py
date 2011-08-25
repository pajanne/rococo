#!/usr/bin/env python
# encoding: utf-8
'''
Created on Oct 16, 2009
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

### ---------------------------------------------------------------------------
### GENPY Import
### ---------------------------------------------------------------------------
import sys, os
from optparse import OptionParser
import logging
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord 
from Bio.Alphabet import IUPAC
# ropy import
import ropy
# genepy import
import util
from setup import connectionFactory, logger
            

### ---------------------------------------------------------------------------
### GENEPY constants
### ---------------------------------------------------------------------------
# platform
IS_LSF = util.isLsf()
bsub_dir = "bsub"
# software
soft_lists = ['fasta35', 'pfscan']
# data
mygenome_dir = "mygenome" 
mygenome_fastafile_allcds = "mygenome_allcds.faa"
refgenomes_dir = "refgenomes"
# results
fasta_dir = "results_fasta"
reciprocalfasta_dir = "results_reciprocalfasta"
refgenomes_extractedseq_dir = "refgenomes_extractedseq"
hamap_dir = "results_hamap"


### ---------------------------------------------------------------------------
### GENEPY steps
### ---------------------------------------------------------------------------
def fasta2embl(infasta):
    """
    Transform sequence file format in fasta to embl using EMBOSS seqret
    Returns the name of created embl file
    """
    util.checkFile(infasta)
    outembl = infasta.split(".")[0] + ".embl"
    """
    Usage: seqret 
    Online documentation: http://emboss.open-bio.org/wiki/Appdoc:Seqret
    
      Standard (Mandatory) qualifiers:
      [-sequence]          seqall     (Gapped) sequence(s) filename and optional
                                      format, or reference (input USA)
      [-outseq]            seqoutall  [<sequence>.<format>] Sequence set(s)
                                      filename and optional format (output USA)
    
      The basic USA syntax is one of:
        "file"
        "file:entry"
        "format::file"
        "format::file:entry"
        "database:entry"
        "database"
        "@file"
    """
    # Create EMBOSS seqret command line
    cmd = "seqret -sequence fasta::%s -outseq embl::%s " % (infasta, outembl)
    # Call the subprocess using convenience method
    util.runProcess(cmd)
    logger.info("File", outembl, "created")
    return outembl

### ---------------------------------------------------------------------------
def concatFeatures(embl, features):
    """
    Concat CDS features in embl format into embl sequence file 
      - the first two lines of embl sequence containing ID & XX lines 
      - the CDS features file containing FT lines
      - the rest of embl sequence containing SQ lines
    Returns the name of created embl sequence file
    """
    util.checkFile(embl)
    util.checkFile(features)
    outembl = embl.split(".")[0] + "_with_cds.embl"
    # Create command line
    head_cmd = "head -2 %s > %s; cat %s >> %s;" % (embl, outembl, features, outembl)
    util.runProcess(head_cmd)
    tail_cmd = "tail +3 %s > tail; cat tail >> %s; rm tail;" % (embl, outembl)
    util.runProcess(tail_cmd)
    logger.info("File", outembl, "created")
    return outembl
    
### ---------------------------------------------------------------------------
def splitSeq(dir, embl, type):
    """
    Split sequence into separate file based on CDS features into dir/ directory
    based on EMBOSS extractfeat
    
    Usage: extractfeat
    Online documentation: http://emboss.open-bio.org/wiki/Appdoc:Extractfeat
    
      Standard (Mandatory) qualifiers:
      [-sequence]          seqall     Sequence(s) filename and optional format, or
                                      reference (input USA)
      [-outseq]            seqout     [.] Sequence filename and
                                      optional format (output USA)
   
      Additional (Optional) qualifiers:
       -type               string     [*] By default every feature in the feature
                                      table is extracted. You can set this to be
                                      any feature type you wish to extract.
                                      See http://www.ebi.ac.uk/Services/WebFeat/
                                      for a list of the EMBL feature types and see
                                      the Uniprot user manual in
                                      http://www.uniprot.org/manual/sequence_annotation
                                      for a list of the Uniprot feature types.
                                      The type may be wildcarded by using '*'.
                                      If you wish to extract more than one type,
                                      separate their names with the character '|',
                                      eg:
                                      *UTR | intron (Any string is accepted)
       -featinname         boolean    [N] To aid you in identifying the type of
                                      feature that has been output, the type of
                                      feature is added to the start of the
                                      description of the output sequence.
                                      Sometimes the description of a sequence is
                                      lost in subsequent processing of the
                                      sequences file, so it is useful for the type
                                      to be a part of the sequence ID name. If
                                      you set this to be TRUE then the name is
                                      added to the ID name of the output sequence.

       Associated qualifiers:
       "-outseq" associated qualifiers
       -ossingle2          boolean    Separate file for each entry
       -ofdirectory2       string     Output directory

      The basic USA syntax is one of:
        "file"
        "file:entry"
        "format::file"
        "format::file:entry"
        "database:entry"
        "database"
        "@file"
    """
    util.checkFile(embl)
    # Create directory
    util.createDir(dir)
    cmd = "extractfeat -sequence embl::%s -type %s -featinname YES -outseq fasta:: -osextension2 ffn -ossingle2 Yes -osdirectory2 %s" % (embl, type, dir)
    util.runProcess(cmd)
    logger.info("Sequences extracted into %s" % dir)
    
### ---------------------------------------------------------------------------
def splitSeqWithBiopython(embl, type):
    """
    Split sequence into separate file based on CDS features into sequences/ directory
    using Biopython
    
    """
    util.checkFile(embl)
    # Create directory sequences/
    dirname = "sequences/"
    util.createDir(dirname)
    record = SeqIO.read(open(embl, "rU"), "embl")
    if len(record.features) == 0:
        sys.exit("ERROR: EMBL file %s without features" % embl)
    for feature in record.features:
        if feature.type == 'CDS':
            seq = record.seq
            
            # Build up a list of (start,end) tuples that will be used to slice the sequence
            locations = []
            # If there are sub_features, then this gene is made up of multiple parts.  
            if len(feature.sub_features): 
                for sf in feature.sub_features:
                    locations.append((sf.location.start.position, sf.location.end.position))
            # This gene is made up of one part.  Store its start and end position.
            else:
                locations.append((feature.location.start.position, feature.location.end.position))

            # Store the joined sequence and nucleotide indices forming the CDS.
            seq_str = '' 
            for begin, end in locations:
                seq_str += seq[begin:end].tostring()

            # Reverse complement the sequence if the CDS is on the minus strand  
            if feature.strand == -1:  
                seq_obj = Seq(seq_str, IUPAC.ambiguous_dna)
                seq_str = seq_obj.reverse_complement().tostring()
            
            logger.debug(feature)
            logger.debug(SeqRecord(seq=Seq(seq_str), id=feature.qualifiers['systematic_id'][0], description=feature.type).format('fasta'))
              
    logger.info("Sequences extracted into %s" % dirname) 
    
### ---------------------------------------------------------------------------
def translateSeq(dir):
    """
    Translate nucleic acid sequence in fasta format into protein sequence using
    EMBOSS transeq
    
    Usage: transeq
    Online documentation: http://emboss.open-bio.org/wiki/Appdoc:Transeq
    
      Standard (Mandatory) qualifiers:
      [-sequence]          seqall     Nucleotide sequence(s) filename and optional
                                      format, or reference (input USA)
      [-outseq]            seqoutall  [.] Protein sequence
                                      set(s) filename and optional format (output
                                      USA)
      Additional (Optional) qualifiers:
       -table              menu       [0] Code to use (Values: 0 (Standard); 1
                                      (Standard (with alternative initiation
                                      codons)); 2 (Vertebrate Mitochondrial); 3
                                      (Yeast Mitochondrial); 4 (Mold, Protozoan,
                                      Coelenterate Mitochondrial and
                                      Mycoplasma/Spiroplasma); 5 (Invertebrate
                                      Mitochondrial); 6 (Ciliate Macronuclear and
                                      Dasycladacean); 9 (Echinoderm
                                      Mitochondrial); 10 (Euplotid Nuclear); 11
                                      (Bacterial); 12 (Alternative Yeast Nuclear);
                                      13 (Ascidian Mitochondrial); 14 (Flatworm
                                      Mitochondrial); 15 (Blepharisma
                                      Macronuclear); 16 (Chlorophycean
                                      Mitochondrial); 21 (Trematode
                                      Mitochondrial); 22 (Scenedesmus obliquus);
                                      23 (Thraustochytrium Mitochondrial))

      The basic USA syntax is one of:
        "file"
        "file:entry"
        "format::file"
        "format::file:entry"
        "database:entry"
        "database"
        "@file"
    """ 
    util.checkDir(dir)
    for file in os.listdir(dir):
        if '.ffn' in file:
            infasta = file
            outpep = file.split(".")[0] + ".faa"
            cmd = "transeq -sequence fasta::%s/%s -outseq fasta::%s/%s -table 11" % (dir, infasta, dir, outpep)
            util.runProcess(cmd)
    logger.info("Sequences translated.")

### ---------------------------------------------------------------------------
def chadoDump(dir):
    """
    Dump the polypeptide sequences of all organisms stored in geneDB/chado in FASTA format
    """
    util.createDir(dir)
    # Connect to geneDB as read only user using ropy.query
    query = ropy.query.QueryProcessor(connection=connectionFactory)
    query.setSQLFilePath(os.path.dirname(__file__) + "/sql/")
    
    # List of organisms
    query.addQueryFromFile("organism_query", "get_all_organisms_with_polyseq.sql")
    organism_rows = query.runQuery("organism_query")
    logger.info("Extracting %s organism sequences from geneDB. Please wait..." % len(organism_rows))
    
    # Add fasta query
    query.addQueryFromFile("fasta_query", "get_fasta_polyseq_for_organism.sql")
    
    for organism in organism_rows:
        organism_name = organism[1]
        organism_id = organism[0]
        if organism_name == "dummy":
            continue
        
        # Dump sequence of each organism into a fasta file
        logger.info("Extracting %s..." % organism_name)
        fasta_rows = query.runQuery("fasta_query", (organism_id, ))
        file_path = "%s/%s_%s.faa" % (dir, organism_id, organism_name)
        out = open(file_path, 'w')
        for row in fasta_rows:
            if not row[0] == None:
                out.write(row[0])
                out.write("\n")
        out.close()
        logger.info("    ...sequence extracted into %s." % file_path)


### ---------------------------------------------------------------------------
def runFasta(seq_dir, genomes_dir, fasta_dir):
    """
    Run FASTA on protein sequences between new genome against all in house genomes
    
    FASTA searches a protein or DNA sequence data bank
     version 35.04 Aug. 25, 2009
     W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448
    """
    util.createDir(fasta_dir)
    # List of in-house genomes
    util.checkDir(genomes_dir)
    genome_files = []
    logger.info("Create fasta results directory for each in-house reference genome")
    for genome_file in os.listdir(genomes_dir):
        if '.faa' in genome_file:
            genome_files.append(genome_file)
            # Create fasta results directory for each in-house genome
            util.createDir("%s/%s" % (fasta_dir, genome_file.split(".")[0]))
            logger.info(genome_file)

    util.checkDir(seq_dir)
    if IS_LSF:
        # Rename new genome sequences for job array to be mygenome_1.faa mygenome_2.faa ...
        seq_num = 0
        for seq_file in os.listdir(seq_dir):
            if not '.faa' in seq_file:
                continue
            seq_num += 1
            if 'mygenome_' in seq_file and '.faa' in seq_file:
                continue
            seq_newfilepath = "%s/mygenome_%s.faa" % (seq_dir, seq_num)
            seq_filepath = "%s/%s" % (seq_dir, seq_file)
            os.rename(seq_filepath, seq_newfilepath)
        # Submit bsub job array on mygenome_${LSB_JOBINDEX}.faa against one refgenome at a time
        bsub_dir = "bsub"
        util.checkDir(bsub_dir)
        for genome_file in genome_files:
            res_dir = "%s/%s" % (fasta_dir, genome_file.split(".")[0])
            cmd = "fasta35 -z 1 -Q -H -S -m 10 %s/mygenome_${LSB_JOBINDEX}.faa %s/%s > %s/mygenome_${LSB_JOBINDEX}.fa" % (seq_dir, genomes_dir, genome_file, res_dir)
            util.submitJobArray(jobname="genepy-fasta", jobnum=seq_num, jobdir=bsub_dir, cmd=cmd)
        util.submitJobDependency('genepy-fasta')
        logger.info("Fasta on LSF finished")
    else:
        # List of new genome sequences
        for seq_file in os.listdir(seq_dir):
            if not '.faa' in seq_file:
                continue
            res_file = seq_file.split(".")[0] + ".fa"
            for genome_file in genome_files:
                res_dir = "%s/%s" % (fasta_dir, genome_file.split(".")[0])
                cmd = "fasta35 -z 1 -Q -H -S -m 10 %s/%s %s/%s > %s/%s" % (seq_dir, seq_file, genomes_dir, genome_file, res_dir, res_file)
                util.runProcess(cmd)
            logger.info(seq_file)
        logger.info("Fasta finished")
        

### ---------------------------------------------------------------------------
def topFastaHits(res_dir, extractedseq_dir):
    """
    Extract top fasta alignment hits that cover at least 80% of the length of 
    both sequences with at least 30% identity.
    Creates an in-house fasta sequence file for each hit
    Returns a dictionnary of hits
    """
    # Identity cutoff for reciprocal searches
    ident_cutoff  = 0.3;
    # Length of hit cutoff for reciprocal searches
    len_cutoff = 0.8;
    # Extracted sequence directory
    util.createDir(extractedseq_dir)
    # TODO Create MSP crunch file
    # Top hits dictionnary
    fastahits_dict = {}
    # Loop over the fasta results
    util.checkDir(res_dir)
    for (path, dirs, files) in os.walk(res_dir):
        for file in files:
            if not '.fa' in file:
                continue
            res_file = path + "/" + file
            logger.info("Reading... " +  res_file)
            # Read the fasta alignment results with biopython AlignIO fasta-m10     
            alignments = AlignIO.parse(open(res_file), "fasta-m10", seq_count=2)
            for alignment in alignments:
                # Select the hit based on cutoffs
                if float(alignment._annotations["sw_ident"]) < ident_cutoff:
                    continue
                record_query = alignment[0]
                record_match = alignment[1]
                overlap = float(alignment._annotations["sw_overlap"])
                if overlap/float(record_query.annotations["original_length"]) < len_cutoff and overlap/float(record_match.annotations["original_length"]) < len_cutoff:
                    continue
                # Create SeqRecord of selected hit
                extractedseq_record = SeqRecord(seq=Seq(str(record_match.seq).replace('-', '')), id=record_match.id, description=res_file)
                extractedseq_file = "%s/%s.faa" % (extractedseq_dir, record_match.id)
                # Print match sequence of selected hit into fasta file
                output_handle = open(extractedseq_file, "w")
                SeqIO.write([extractedseq_record], output_handle, "fasta")
                output_handle.close()
                logger.info("    ...sequence extracted into %s" % extractedseq_file)
                record_query_region = "%s-%s" % (record_query._al_start, record_query._al_stop)
                record_match_region = "%s-%s" % (record_match._al_start, record_match._al_stop)
                # add hit into dictionnary
                key = "%s||%s" % (record_query.id, record_match.id)
                # value in MSP crunch format
                value = "%s %s %s %s %s %s" % (alignment._annotations["sw_score"], alignment._annotations["sw_ident"], record_query_region, record_query.id, record_match_region, record_match.id)
                fastahits_dict[key] = value
    logger.info("Extract fasta alignment hits finished")
    return fastahits_dict
        
### ---------------------------------------------------------------------------
def concatSeq(genome_file, dir):
    """
    Concatenate separated CDS sequence fasta files located in dir into one file
    """
    util.checkDir(dir)
    if os.path.exists(genome_file):
        os.remove(genome_file)
    cmd = "cat %s/*.faa > %s" % (dir, genome_file)
    util.runProcess(cmd)
    logger.info("concatSeq finished")
        
### ---------------------------------------------------------------------------
def runReciprocalFasta(seq_dir, genome_file, fasta_dir):
    """
    Run FASTA between extracted in-house protein sequences against new genome 
    
    FASTA searches a protein or DNA sequence data bank
     version 35.04 Aug. 25, 2009
     W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448
    """
    util.createDir(fasta_dir)
    # Check new genome
    util.checkFile(genome_file)
    # Check ref genome extracted sequences
    util.checkDir(seq_dir)
    res_dir = fasta_dir
    if IS_LSF:
        # Rename new genome sequences for job array to be refgenome_1.faa refgenome_2.faa ...
        seq_num = 0
        for seq_file in os.listdir(seq_dir):
            if not '.faa' in seq_file:
                continue
            seq_num += 1
            if 'refgenome_' in seq_file and '.faa' in seq_file:
                continue
            seq_newfilepath = "%s/refgenome_%s.faa" % (seq_dir, seq_num)
            seq_filepath = "%s/%s" % (seq_dir, seq_file)
            os.rename(seq_filepath, seq_newfilepath)
        # Submit bsub job array on refgenome_${LSB_JOBINDEX}.faa against mygenome
        bsub_dir = "bsub"
        util.checkDir(bsub_dir)
        cmd = "fasta35 -z 1 -Q -H -S -m 10 %s/refgenome_${LSB_JOBINDEX}.faa %s > %s/refgenome_${LSB_JOBINDEX}.fa" % (seq_dir, genome_file, res_dir)
        util.submitJobArray(jobname="genepy-recipfasta", jobnum=seq_num, jobdir=bsub_dir, cmd=cmd)
        util.submitJobDependency('genepy-recipfasta')
        logger.info("Reciprocal Fasta on LSF finished")
    else:
        # List of inhouse extracted genome sequences
        for seq_file in os.listdir(seq_dir):
            if not '.faa' in seq_file:
                continue
            res_file = seq_file.split(".")[0] + ".fa"
            cmd = "fasta35 -z 1 -Q -H -S -m 10 %s/%s %s > %s/%s" % (seq_dir, seq_file, genome_file, res_dir, res_file)
            util.runProcess(cmd)
            logger.info(seq_file)
        logger.info("Reciprocal Fasta finished")

### ---------------------------------------------------------------------------
def topReciprocalFastaHits(res_dir):
    """
    Extract top hits that cover at least 80% of the length of both sequences
    with at least 30% identity.
    Returns a dictionary of hits
    """
    # Identity cutoff for reciprocal searches
    ident_cutoff  = 0.3;
    # Length of hit cutoff for reciprocal searches
    len_cutoff = 0.8;
    # TODO Create MSP crunch file
    # Top hits dictionnary
    fastahits_dict = {}
    # Loop over the fasta results
    util.checkDir(res_dir)
    for (path, dirs, files) in os.walk(res_dir):
        for file in files:
            if not '.fa' in file:
                continue
            res_file = path + "/" + file
            logger.info("Reading... " +  res_file)
            # Read the fasta alignment results with biopython AlignIO fasta-m10     
            alignments = AlignIO.parse(open(res_file), "fasta-m10", seq_count=2)
            for alignment in alignments:
                # Select the hit based on cutoffs
                if float(alignment._annotations["sw_ident"]) < ident_cutoff:
                    continue
                record_query = alignment[0]
                record_match = alignment[1]
                overlap = float(alignment._annotations["sw_overlap"])
                if overlap/float(record_query.annotations["original_length"]) < len_cutoff and overlap/float(record_match.annotations["original_length"]) < len_cutoff:
                    continue
                
                record_query_region = "%s-%s" % (record_query._al_start, record_query._al_stop)
                record_match_region = "%s-%s" % (record_match._al_start, record_match._al_stop)
                # add hit into dictionnary
                key = "%s||%s" % (record_match.id, record_query.id) # inverted key to be comparable with fasta hits
                value = "%s %s %s %s %s %s" % (alignment._annotations["sw_score"], alignment._annotations["sw_ident"], record_query_region, record_query.id, record_match_region, record_match.id)
                fastahits_dict[key] = value
    logger.info("Extract reciprocal fasta alignment hits finished")
    return fastahits_dict

### ---------------------------------------------------------------------------
def printMSPCrunch(fasta_hits, reciprocal_hits):
    """
    Print an MSPCrunch format description of the reciprocal hit
    """
    for reciprocal_key in reciprocal_hits.keys():
        if fasta_hits.has_key(reciprocal_key):
            logger.info(fasta_hits[reciprocal_key]) 
    logger.info("MSP Crunch extracted")

### ---------------------------------------------------------------------------
def getHits(fasta_hits, reciprocal_hits):
    """
    Return two dictionaries of ortholog hits and similarity hits containing
    {'new_genome_CDS_name':inhouse_genome_feature_id}
    """
    ortholog_hits = {}
    for reciprocal_key in reciprocal_hits.keys():
        if fasta_hits.has_key(reciprocal_key):
            ortholog_hits[reciprocal_key.split("||")[0]] = reciprocal_key.split("||")[1]
            del fasta_hits[reciprocal_key]
    similarity_hits = {}
    for fasta_key in fasta_hits:
        new_genome_key = fasta_key.split("||")[0]
        if not ortholog_hits.has_key(new_genome_key):
            similarity_hits[new_genome_key] = fasta_key.split("||")[1]
    return {'ortholog':ortholog_hits, 'similarity':similarity_hits}
    logger.info("Hits processed")
        
    
### ---------------------------------------------------------------------------
def transferFeatures(hits):
    """
    In table: feature_cvterm
    RILEY              /class
    genedb_products    /product
    
    In table: featureprop
    EC_number          /EC_number
    colour             /colour
    gene               /gene
    """
    # Connect to geneDB as read only user using ropy.query
    query = ropy.query.QueryProcessor(connection=connectionFactory)
    query.setSQLFilePath(os.path.dirname(__file__) + "/sql/")
    
    for hit in hits:
        # Extract all cvterm related to a feature_id from feature_cvterm table
        query.addQueryFromFile("feature_cvterm_query", "get_cvterm_from_feature_cvterm.sql")
        feature_cvterm_rows = query.runQuery("feature_cvterm_query", (hits[hit],))
        logger.debug("--- %s" % hit)
        logger.debug('/ortholog="%s"' % hits[hit])
        for row in feature_cvterm_rows:
            cvterm_name = row[0]
            cv_name = row[1]
            if cv_name == "RILEY":
                logger.debug('/class="%s"' % (cvterm_name))
            elif cv_name == "genedb_products":
                logger.debug('/product="%s"' % (cvterm_name))
        # Extract all cvterm relected to a feature_id from featureprop
        query.addQueryFromFile("featureprop_query", "get_cvterm_from_featureprop.sql")
        featureprop_rows = query.runQuery("featureprop_query", (hits[hit],))
        for row in featureprop_rows:
            logger.debug('/%s="%s"' % (row[0], row[1]))
    logger.info("Features transfered")

### ---------------------------------------------------------------------------
def runHamapScan(seq_dir, hamap_dir):
    """
    HAMAP: High-quality Automated and Manual Annotation of microbial Proteomes
    ftp download site: ftp://ftp.expasy.org/databases/hamap/
     
    pfscan compares a protein or nucleic acid sequence against a profile 
    library. The result is an unsorted list of profile-sequence matches.
    download site: http://www.isrec.isb-sib.ch/ftp-server/pftools/pft2.3/
    """
    util.createDir(hamap_dir)
    util.checkDir(seq_dir)
    hamap_profile_file = "%s/hamap/hamap.prf" % os.path.dirname(__file__)
    if IS_LSF:
        # Rename new genome sequences for job array to be mygenome_1.faa mygenome_2.faa ...
        seq_num = 0
        for seq_file in os.listdir(seq_dir):
            if not '.faa' in seq_file:
                continue
            seq_num += 1
            if 'mygenome_' in seq_file and '.faa' in seq_file:
                continue
            seq_newfilepath = "%s/mygenome_%s.faa" % (seq_dir, seq_num)
            seq_filepath = "%s/%s" % (seq_dir, seq_file)
            os.rename(seq_filepath, seq_newfilepath)
        # Submit bsub job array on mygenome_${LSB_JOBINDEX}.faa against hamap profile
        bsub_dir = "bsub"
        util.checkDir(bsub_dir)
        cmd = "pfscan -klf %s/mygenome_${LSB_JOBINDEX}.faa %s > %s/mygenome_${LSB_JOBINDEX}.out" % (seq_dir, hamap_profile_file, hamap_dir)
        util.submitJobArray(jobname='genepy-hamap', jobnum=seq_num, jobdir=bsub_dir, cmd=cmd)
        util.submitJobDependency('genepy-hamap')
        logger.info("HAMAP scan on LSF finished")
    else:
        # List of new genome sequences
        for seq_file in os.listdir(seq_dir):
            if not '.faa' in seq_file:
                continue
            res_file = seq_file.split(".")[0] + ".out"
            cmd = "pfscan -klf %s/%s %s > %s/%s" % (seq_dir, seq_file, hamap_profile_file, hamap_dir, res_file)
            util.runProcess(cmd)
        logger.info("HAMAP scan finished")


### ---------------------------------------------------------------------------
### GENEPY main
### ---------------------------------------------------------------------------
def main():
    # Fasta file extension: 
    # .ffn for the untranslated nucleotide sequences for each CDS; .faa for protein coding sequences (CDS)
    # .fa for the fasta alignment results
    # .fna for whole genomic DNA sequences; .frn for nucleotide sequences of RNA related features
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-d", "--dna", metavar="FILE", help="input dna FILE in fasta format", action="store", type="string", dest="dna")
    parser.add_option("-t", "--tab", metavar="FILE", help="input tab FILE in embl format", action="store", type="string", dest="tab")
    parser.add_option("-e", "--embl", metavar="FILE", help="input embl FILE with CDS features in embl format", action="store", type="string", dest="embl")
    parser.add_option("--genedb", help="extract reference genome protein sequences from geneDB", action="store_true", dest="db")
    parser.add_option("--fasta", help="run fasta against each extracted in-house genomes", action="store_true", dest="fasta")
    parser.add_option("--hamap", help="run pfscan against HAMAP profiles", action="store_true", dest="hamap")
    parser.add_option("--clean", help="delete all results without deleting reference genomes", action="store_true", dest="clean")
    parser.add_option("--deepclean", help="delete all reference genomes and results", action="store_true", dest="deepclean")
    (options, args) = parser.parse_args()
    # Print help if no argument given
    if util.printHelp(options):
        parser.print_help()
        sys.exit()
    # Print command line
    cmdline = "$ python "
    for argv in sys.argv:
        cmdline += argv + " " 
    logger.debug(cmdline)
    
    # >>> ---------------------------------------------------------------------
    # >>> DATA PREPARATION
    # >>> ---------------------------------------------------------------------
    # List of needed software
    for softname in soft_lists:
        util.checkSoft(softname)
    # Prepare new genome data
    if options.dna and options.tab and not options.embl:
        util.checkFile(options.dna)
        mygenome_emblfile = fasta2embl(options.dna)
        mygenome_emblfile_withcds = concatFeatures(mygenome_emblfile, options.tab)
        splitSeq(mygenome_dir, mygenome_emblfile_withcds, "CDS")
        translateSeq(mygenome_dir)
    elif not options.dna and not options.tab and options.embl:
        mygenome_emblfile_withcds = options.embl
        splitSeq(mygenome_dir, mygenome_emblfile_withcds, "CDS")
        #splitSeqWithBiopython(mygenome_emblfile_withcds, "CDS") # does not work with testdata_01
        translateSeq(mygenome_dir)
    elif not options.deepclean:
        util.checkDir(mygenome_dir)
    # Extract in house genomes from chado db
    if options.db:
        chadoDump(refgenomes_dir)
    elif not options.deepclean:
        util.checkDir(refgenomes_dir)
    # bsub output directory
    if IS_LSF and not (options.clean or options.deepclean):
        util.createDir(bsub_dir)

    # >>> ---------------------------------------------------------------------
    # >>> ORTHOLOG SEARCH
    # >>> ---------------------------------------------------------------------
    # Run fasta & reciprocal fasta
    if options.fasta:
        runFasta(mygenome_dir, refgenomes_dir, fasta_dir)
        fasta_hits = topFastaHits(fasta_dir, refgenomes_extractedseq_dir)
        concatSeq(mygenome_fastafile_allcds, mygenome_dir)
        runReciprocalFasta(refgenomes_extractedseq_dir, mygenome_fastafile_allcds, reciprocalfasta_dir)
        reciprocalfasta_hits = topReciprocalFastaHits(reciprocalfasta_dir)
        printMSPCrunch(fasta_hits, reciprocalfasta_hits)
        hits = getHits(fasta_hits, reciprocalfasta_hits)
        logger.info("ORTHOLOGS")
        logger.info(hits['ortholog'])
        logger.info("SIMILARITY")
        logger.info(hits['similarity'])
        transferFeatures(hits['ortholog'])
    # Run hamap scan
    if options.hamap:
        runHamapScan(mygenome_dir, hamap_dir)

    # >>> ---------------------------------------------------------------------
    # >>> CLEANING OUTPUT DATA
    # >>> ---------------------------------------------------------------------
    # Clean results before a re-run
    if options.clean:
        # fasta results
        util.rmDir(fasta_dir)
        util.rmDir(reciprocalfasta_dir)
        util.rmDir(refgenomes_extractedseq_dir)
        util.rmFile(mygenome_fastafile_allcds)
        # hamap results
        util.rmDir(hamap_dir)
        # bsub outputs
        if IS_LSF:
            util.rmDir(bsub_dir)
    # Deep clean - remove all
    if options.deepclean:
        util.rmDir(refgenomes_dir)
        util.rmDir(mygenome_dir)
        util.rmDir(fasta_dir)
        util.rmDir(reciprocalfasta_dir)
        util.rmDir(refgenomes_extractedseq_dir)
        util.rmFile(mygenome_fastafile_allcds)
        util.rmDir(hamap_dir)

    # >>> ---------------------------------------------------------------------
    # >>> SIMILARITY SEARCH
    # >>> ---------------------------------------------------------------------
    #runBlastP() # BlastP simply compares a protein query to a protein database.
    #runMeropsScan() # Peptidases database
    
    # >>> ---------------------------------------------------------------------
    # >>> PROTEIN DOMAIN SEARCH
    # >>> ---------------------------------------------------------------------
    #runPrositeScan() # Database of protein domains, families and functional sites
    #runInterProScan() # Database of protein families, domains, regions, repeats and sites
    
    # >>> ---------------------------------------------------------------------
    # >>> STRUCTURE SEARCH
    # >>> ---------------------------------------------------------------------
    #runHelixTurnHelix()
    
    # >>> ---------------------------------------------------------------------
    # >>> HYDROPHOBIC FEATURES SEARCH
    # >>> ---------------------------------------------------------------------
    #runSignalP() # predicts the presence and location of signal peptide cleavage sites in amino acid sequences
    #runTMHMM() # Prediction of transmembrane helices in proteins
    
    # >>> ---------------------------------------------------------------------
    # >>> CLUSTERING
    # >>> ---------------------------------------------------------------------
    #runTribeMCL()
    
if __name__ == '__main__':
    main()
    
