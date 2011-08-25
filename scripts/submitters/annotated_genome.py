'''
Created on Mar 2, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.

For submitting EMBL annotated genomes to EBI

http://www.ebi.ac.uk/embl/Submission/genomes.html
http://www.ebi.ac.uk/embl/Documentation/User_manual/usrman.html
http://www.ebi.ac.uk/embl/Documentation/FT_definitions/feature_table.html
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
FTP_URL = 'ftp-private.ebi.ac.uk'
FTP_USERNAME = 'enaftp'
FTP_PASSWORD = 'L9Nwoqoz'
FTP_DIR = {'metahit':'mH4kl9FJh93HKqPoY63n/to_ena',
           'other':'hjGkHF83J33Mrt98TYaz/to_ena'
           }

CONTACTS = {}
CONTACTS['ap12'] = {'author':'Pajon A.'}
CONTACTS['maa'] = {'author':'Aslett M.'}

### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------
def doReplaceHeader(dep_file, locus_tag, ac):
    """
Replace

ID   XXX; XXX; linear; DNA; XXX; ; 3763317 BP.
XX
AC   XXX;

by

ID   XXX; XXX; linear; XXX; XXX; XXX; 3763317 BP.
XX
AC   ;
XX
AC * _AL1001

"""
    r_file = open(dep_file, 'r')
    lines = r_file.readlines()
    r_file.close()
    w_file = open(dep_file,'w')
    for line in lines:
        if line.startswith('ID'):
            line = line.replace(' ;', 'XXX;')
            line = line.replace('DNA', 'genomic DNA').replace('WGS', 'XXX')
            #line = line.replace('linear', 'XXX').replace('DNA', 'XXX')
            w_file.write(line)
        elif line.startswith('AC') and ac == 'XXX':
            line = line.replace('XXX', '')
            w_file.write(line)
            w_file.write('XX\n')
            w_file.write('AC * _%s001\n' % locus_tag)
        else:
            w_file.write(line)
    w_file.close()


### ---------------------------------------------------------------------------
def getMolType(embl_file):
    """
Extract mol_type from ID line

ID   XXX;XXX;XXX; genomic DNA;XXX; PRO; 1981535 BP.
"""
    f = open(embl_file, 'r')
    line = f.readline()
    assert line[:5].rstrip() == "ID" 
    fields = [data.strip() for data in line[5:].strip().split(";")] 
    assert len(fields) == 7 
    """ 
    The tokens represent: 
     0. Primary accession number 
     1. Sequence version number 
     2. Topology: 'circular' or 'linear' 
     3. Molecule type (e.g. 'genomic DNA') 
     4. Data class (e.g. 'STD') 
     5. Taxonomic division (e.g. 'PRO') 
     6. Sequence length (e.g. '4639675 BP.') 
    """ 
    f.close()
    return fields[3]


### ---------------------------------------------------------------------------
def doConvert(embl_file, dep_file, contact, project, genome_project_id, organism_name, strain, locus_tag, taxon_id, dna_source, authors, comment, ac, clean=False):
    record = SeqIO.read(open(embl_file), "embl")

    # ----------------------------------------
    # HEADER
    # ----------------------------------------
    # remove accession
    if 'accession' in record.annotations.keys():
        del record.annotations['accession']
    record.annotations['accession'] = [ac]
    # ID line
    record.id = "XXX"
    record.name = "XXX"
    record.annotations['data_file_division'] = 'PRO'
    record.annotations['data_file_class'] = 'WGS'
    # PR line
    record.dbxrefs = ["Project:%s" % genome_project_id]
    # OS line
    record.annotations["organism"] = "%s %s" % (organism_name, strain)
    # DE line
    if project == 'metahit':
        record.description = "%s %s draft genome." % (organism_name, strain)
    else:
        record.description = "%s %s genome." % (organism_name, strain)
    # RN & RL lines
    if dna_source == 'GHP':
        dna_source = 'Rowett Institute of Nutrition and Health, University of Aberdeen -- http://www.rowett.ac.uk/divisions/ghp/'
        authors = 'Pajon A., Turner K., Parkhill J., Duncan S., Flint H.'
    elif dna_source == 'INRA':
        dna_source = 'INRA Clermont-Ferrand-Theix -- http://www.clermont.inra.fr/'
        authors = 'Pajon A., Turner K., Parkhill J., Bernalier A.'
    elif dna_source == 'HCIR':
        dna_source = 'Helmholtz Centre for Infection Research -- http://www.helmholtz-hzi.de/'
        authors = 'Pajon A., Turner K., Parkhill J., Timmis K., Oxley A., Wurdemann D.'
    elif dna_source == 'DSMZ':
        dna_source = 'German Collection of Microorganisms and Cell Cultures -- http://www.dsmz.de/'
        authors = 'Pajon A., Turner K., Parkhill J.'
    elif dna_source == 'NCTC':
        dna_source = 'Health Protection Agency\'s National Collection of Type Cultures -- http://www.hpacultures.org.uk/'
        authors = 'Pajon A., Turner K., Parkhill J.'
    elif dna_source == 'DPM':
        dna_source = 'Departments of Periodontology and Microbiology, King\'s College London -- http://www.kcl.ac.uk/'
        authors = 'Pajon A., Turner K., Parkhill J., Wade W., Vartoukian S.'
    else:
        dna_source = dna_source
        authors = authors
    ref_journal = Reference()
    ref_journal.journal = 'Unpublished.'
    if project == 'metahit':
        ref_journal.consrtm = "metaHIT consortium -- http://www.metahit.eu/"
    ref_journal.title = 'The genome sequence of %s %s' % (organism_name, strain)
    ref_journal.authors = authors
    ref_dep = Reference()
    ref_dep.authors = CONTACTS[contact]['author']
    today = date.today()
    ref_dep.journal = "Submitted (%s) to the EMBL/GenBank/DDBJ databases. Sanger Institute, Wellcome Trust Genome Campus, Hinxton, Cambridge CB10 1SA, United Kingdom." % today.strftime("%d-%b-%Y")
    ref_dep.title = 'Direct submission'
    record.annotations['references'] = [ref_journal, ref_dep]
    # CC line
    record.annotations['comment'] = ['Data release policy http://www.sanger.ac.uk/legal/#t_2',
                                     'DNA source: %s' % dna_source,
                                     '%s' % comment]
    
    # ----------------------------------------
    # GAP FEATURE (only with clean option)
    # ----------------------------------------
    # Add FT gap 
    seq = record.seq
    in_N = False
    gap_features = []
    if clean:
        # TODO - Cope with a sequence which ends with N
        if seq[-1] != "N":
            print "WARNING: sequence ends with N"
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
                    gap_feature = SeqFeature(FeatureLocation(start_N,end_N), strand=1, type="gap")
                    gap_feature.qualifiers['estimated_length'] = [length]
                    gap_features.append(gap_feature)
                in_N = False
    
    # ----------------------------------------
    # OTHER FEATURE (only with clean option)
    # ----------------------------------------
    new_features = []
    first_source = True
    has_source = False
    removed_cds = 0
    for i in range(len(record.features)):
        feature = record.features[i]

        # Add strain into FT source 
        if feature.type == 'source' and first_source:
            has_source = True
            feature.location.end.position = len(record.seq)
            feature.qualifiers['organism'] = ["%s %s" % (organism_name, strain)]
            feature.qualifiers['strain'] = [strain]

        # Remove qualifier /note & /translation
        if clean:
            if 'note' in feature.qualifiers.keys():
                del feature.qualifiers['note']
            #if 'translation' in feature.qualifiers.keys():
            #    del feature.qualifiers['translation']

        # Rename locus_tag
        if clean:
            if 'locus_tag' in feature.qualifiers.keys():
                feature.qualifiers['locus_tag'] = [getLocusTag(feature.qualifiers['locus_tag'][0], locus_tag, feature.type)]

        # Check /EC_number="5.3.1.25" or /EC_number="1.1.2.-"
        if clean:
            if 'EC_number' in feature.qualifiers.keys():
                for i in range(len(feature.qualifiers['EC_number'])):
                    feature.qualifiers['EC_number'][i] = getEcNumber(feature.qualifiers['EC_number'][i])
            # Remove (EC 2.1.2.3) in /product and /function
            if 'product' in feature.qualifiers.keys():
                for i in range(len(feature.qualifiers['product'])):
                    (feature.qualifiers['product'][i], ec_list) = getValueWithoutEc(feature.qualifiers['product'][i], feature)
                    if ec_list:
                        for ec in ec_list:
                            if 'EC_number' not in feature.qualifiers.keys():
                                feature.qualifiers['EC_number'] = [ec]
                            else:
                                feature.qualifiers['EC_number'].append(ec)

        # Remove tRNA /product (not only when containing ???)
        if clean:
            if feature.type == 'tRNA' and 'product' in feature.qualifiers.keys():
                del feature.qualifiers['product']
                #for i in range(len(feature.qualifiers['product'])):
                #    if feature.qualifiers['product'][i].count('?') > 1:
                #        del feature.qualifiers['product'][i]
            if 'function' in feature.qualifiers.keys():
                for i in range(len(feature.qualifiers['function'])):
                    (feature.qualifiers['function'][i], ec_list) = getValueWithoutEc(feature.qualifiers['function'][i], feature)
                    if ec_list:
                        for ec in ec_list:
                            if 'EC_number' not in feature.qualifiers.keys():
                                feature.qualifiers['EC_number'] = [ec]
                            else:
                                feature.qualifiers['EC_number'].append(ec)

        # Remove FT gene & keep only one FT source per record & remove some CDS
        if clean:
            if not feature.type == 'gene':
                if feature.type == 'source':
                    if first_source:
                        new_features.append(feature)
                        first_source = False
                # Remove CDS that are not valid
                # CDS -- translation must start with M, nucleotide sequence without N's & no overlap with gap feature
                # CDS -- must not have internal stop codons
                # CDS -- must end with stop codons (TAG, TAA, or TGA) or add '<' or '>' e.g. complement(<1..174); 1399953..>1401221
                elif feature.type == 'CDS':
                    if not 'transl_table' in feature.qualifiers.keys():
                        feature.qualifiers['transl_table'] = 11
                    if feature.strand == 1:
                        stop_codon = record.seq[feature.location.nofuzzy_end-3:feature.location.nofuzzy_end]
                        if not str(stop_codon) in ['TAG', 'TAA', 'TGA']:
                            feature.location = FeatureLocation(ExactPosition(feature.location.nofuzzy_start), AfterPosition(feature.location.nofuzzy_end))
                    if feature.strand == -1:
                        stop_codon = record.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_start+3]
                        if not str(stop_codon) in ['CTA', 'TTA', 'TCA']:
                            feature.location = FeatureLocation(BeforePosition(feature.location.nofuzzy_start), ExactPosition(feature.location.nofuzzy_end))
                    translation = feature.extract(record.seq).translate(table=11)
                    if 'translation' in feature.qualifiers.keys():
                        if translation[-1] == '*':
                            if not len(translation) - 1 == len(feature.qualifiers['translation'][0]):
                                print 'WARNING: CDS %s translation length of different size' % feature.location
                                print translation
                                print feature.qualifiers['translation'][0]
                            else:
                                if not str(translation[:-1]) == str(feature.qualifiers['translation'][0]):
                                    print 'WARNING: CDS %s translation not identical' % feature.location
                                    print translation[:-1]
                                    print feature.qualifiers['translation'][0]
                        else:
                            if not len(translation) == len(feature.qualifiers['translation'][0]):
                                print 'WARNING: CDS %s translation length of different size' % feature.translation
                                print translation
                                print feature.qualifiers['translation'][0]
                            else:
                                if not str(translation) == str(feature.qualifiers['translation'][0]):
                                    print 'WARNING: CDS %s translation not identical' % feature.location
                                    print translation
                                    print feature.qualifiers['translation'][0]
                    #feature.qualifiers['translation'] = [translation]
                    if translation.startswith('M') and record.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_end].count('N') == 0:
                        if translation[:-1].count('*') >= 1:
                            print 'WARNING: CDS %s with internal stop codon' % feature.location
                            print translation
                        else:
                            new_features.append(feature)
                    else:
                        print 'WARNING: CDS %s does not start with M' % feature.location
                        print translation
                        removed_cds = removed_cds + 1
                else:
                    new_features.append(feature)

        if not clean:
             new_features.append(feature)

    # Add source feature
    if not has_source:
        feature = SeqFeature(FeatureLocation(0,len(record.seq)), type="source")
        feature.qualifiers['organism'] = ["%s %s" % (organism_name, strain)]
        feature.qualifiers['strain'] = [strain]
        feature.qualifiers['db_xref'] = ["taxon:%s" % taxon_id]
        feature.qualifiers['mol_type'] = getMolType(embl_file) 
        new_features.append(feature)

    if clean:
        print 'WARNING: %s CDSs have been removed' % removed_cds
    else:
        print "Only adding header, use '--clean' for cleaning features"
    
    new_features.extend(gap_features)
    new_features.sort(feature_compare)
    record.features = new_features
    
    # Write out new embl file
    SeqIO.write([record], open(dep_file, "w"), "embl")

### ---------------------------------------------------------------------------
def getValueWithoutEc(value, feature):
    # Remove (EC 2.1.2.3)
    pattern = '\(EC( |:)([0-9]|\-|\.)*?\)'
    ec_list = []
    if re.search(pattern, value):
        splitted_values = value.split('(EC')
        for x in splitted_values:
            if x.count(')') >= 1 and x.count('.') >= 1:
                ec_number = getEcNumber(x.split(')')[0]).strip().replace(':','')
                ec_number_found = False
                if 'EC_number' in feature.qualifiers.keys():
                    for i in range(len(feature.qualifiers['EC_number'])):
                        if ec_number == feature.qualifiers['EC_number'][i]:
                            ec_number_found = True
                if not ec_number_found:
                    ec_list.append(ec_number)
        return re.sub(pattern, '', value).strip(), ec_list
    else:
        return value, ec_list
 
### ---------------------------------------------------------------------------
def getEcNumber(value):
    # /EC_number="5.3.1.25" or /EC_number="1.1.2.-"
    if value.count(".") == 3:
        if value.endswith('.'):
            return "%s-" % value
        else:
            return value
    # /EC_number="5.-"
    elif value.count('.') == 1:
        if value.endswith('-'):
            return "%s.-.-" % value
        elif value.endswith('.'): 
            return "%s-.-.-" % value
        else:
            print "EC_number case not supported: %s" % value
            return value
    # /EC_number="2.7.-"
    elif value.count('.') == 2:
        if value.endswith('-'):
            return "%s.-" % value
        elif value.endswith('.'): 
            return "%s-.-" % value
        else:
            print "EC_number case not supported: %s" % value
            return value
    else:
        print "EC_number case not supported: %s" % value
        return value

### ---------------------------------------------------------------------------
def getLocusTag(value, locus_tag, feature_type):
    if feature_type == "tRNA":
        return '%s_T_%5s' % (locus_tag, value[len(value)-5:len(value)])
    elif feature_type == "rRNA":
        return '%s_R_%5s' % (locus_tag, value[len(value)-5:len(value)])
    else:
        return '%s_%5s' % (locus_tag, value[len(value)-5:len(value)])
        
### ---------------------------------------------------------------------------
def getLocation(feature):
    if feature.strand == -1:
        return "%s\t%s" % (feature.location.end.position, int(feature.location.start.position) + 1)
    else:
        return "%s\t%s" % (int(feature.location.start.position) + 1, feature.location.end.position)

### ---------------------------------------------------------------------------
def feature_compare(fx, fy):
    if int(fx.location.start.position) > int(fy.location.start.position):
        return 1
    elif int(fx.location.start.position) == int(fy.location.start.position):
        return 0
    else: 
        return -1

### ---------------------------------------------------------------------------
def doValidate(dep_file, val_file):
    # java -cp embl-client-DEV.jar uk.ac.ebi.client.EnaValidator -r <files>
    cmd = "java -classpath %s/validators/embl-client-DEV.jar uk.ac.ebi.client.EnaValidator -r %s &> %s" % (os.path.realpath(os.path.dirname(__file__)).replace('/submitters', ''), dep_file, val_file)
    
    process = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process_out, process_err = process.communicate()
    retval = process.poll()
    if not retval == 0:
        print "ERROR executing %s command: %s" % (cmd.split()[0], process_err)
        sys.exit(1)
    else:
        print process_err
        print process_out

### ---------------------------------------------------------------------------
def doSubmit(file, project):
    server = ftplib.FTP(FTP_URL, FTP_USERNAME, FTP_PASSWORD)
    server.cwd(FTP_DIR[project])
    
    f = open(file,'rb')
    print server.storlines('STOR %s' % file, f)
    
    f.close()
    server.quit()


### ---------------------------------------------------------------------------
def doRun(list, project, contact, clean=False, submit=False):
    print "Reading input file %s" % list
    processed_lines = 0
    for line in open(list, "r"):
        if line[0] == '!':
            continue
        if not line.count('||') == 9:
            continue
        # embl_file||genome_project_id||organim_name||strain||locus_tag||yaxon_id||dna_source||authors||comment||ac
        processed_lines += 1 
        line = line.strip()
        values = line.split('||')
        embl_file=values[0]
        genome_project_id = values[1]
        organism_name = values[2]
        strain = values[3]
        locus_tag = values[4]
        taxon_id = values[5]
        dna_source = values[6]
        authors = values[7]
        comment = values[8]
        ac = values[9]

        # Check if EMBL file exists
        if not os.path.exists(embl_file):
            print "ERROR: file %s does not exist" % embl_file
            sys.exit(1)
            
        # Create new EMBL file for deposition
        dep_file = "%s.4dep" % embl_file
        val_file = "%s.4dep.val" % embl_file
        
        if not submit:
            print "Converting file %s into %s before submission..." % (embl_file, dep_file)
            doConvert(embl_file, dep_file, contact, project, genome_project_id, organism_name, strain, locus_tag, taxon_id, dna_source, authors, comment, ac, clean)
            
            print "Replacing header in file %s before submission..." % (dep_file)
            doReplaceHeader(dep_file, locus_tag, ac)
            
            print "Validating file %s before submission..." % (dep_file)
            doValidate(dep_file, val_file)
            
            print "Only converting files, use '--submit' for submitting data"
        
        else:
            if not os.path.exists(dep_file):
                print "ERROR: file %s does not exist" % dep_file
                print "Please convert files before submitting them, remove '--submit' for converting data"
                sys.exit(1)
            print 'Submitting file %s to ftp://%s:%s@%s/%s' % (dep_file, FTP_USERNAME, FTP_PASSWORD, FTP_URL, FTP_DIR[project])
            doSubmit(dep_file, project)

    print "%s processed lines from input file" % processed_lines

### ---------------------------------------------------------------------------
def main():
    from optparse import OptionParser

    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", "--list", metavar="FILE", help="FILE containing the list of all organism names, and its associated information (embl_file||genome_project_id||organim_name||strain||locus_tag||taxon_id||dna_source||authors||comment||ac) you wish to submit", action="store", type="string", dest="list")
    parser.add_option("-c", "--contact", metavar="CONTACT", help="name of the contact %s submitting the data" % CONTACTS.keys(), action="store", choices=CONTACTS.keys(), dest="contact")
    parser.add_option("-p", "--project", metavar="PROJECT", help="name of the project %s defining to which ftp directory the data should go" % FTP_DIR.keys(), action="store", choices=FTP_DIR.keys(), dest="project")
    parser.add_option("--submit", help="To submit data, not only converting and validating EMBL files", action="store_true", dest="submit")
    parser.add_option("--clean", help="To clean features, not only adding header and validating EMBL files", action="store_true", dest="clean")

    (options, args) = parser.parse_args()
    
    # Print help if no argument given
    if not (options.list and options.contact and options.project):
        parser.print_help()
        sys.exit()
    
    # Check proxy settings
    if not os.environ.has_key('http_proxy'):
        print 'You may need to set your http_proxy env variable (e.g. setenv http_proxy "http://wwwcache.sanger.ac.uk:3128")'

    # Check input file list exists
    list = options.list
    if not os.path.exists(list):
        print "ERROR: file %s does not exist" % list
        sys.exit(1)    

    doRun(list, options.project, options.contact, options.clean, options.submit)
            
### ---------------------------------------------------------------------------
if __name__ == '__main__':
    main()
