'''
Created on Feb 11, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.

For submitting genome projects to NCBI 

'''

import sys, os
import urllib, urllib2

### ---------------------------------------------------------------------------
### constants
### ---------------------------------------------------------------------------
LOCUSTAG_URL = 'http://www.ncbi.nlm.nih.gov/genomes/checkltp.cgi'
PROJECT_URL = 'http://www.ncbi.nlm.nih.gov/genomes/mpfsubmission.cgi'

CONTACTS = {}
CONTACTS['ap12'] = {'firstname':'Anne', 'lastname':'Pajon', 'email':'ap12@sanger.ac.uk'}
CONTACTS['maa'] = {'firstname':'Martin', 'lastname':'Aslett', 'email':'maa@sanger.ac.uk'}

### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------
def isLocusTagAvailable(locus_tag):
    opener = urllib2.build_opener(urllib2.HTTPCookieProcessor())
    urllib2.install_opener(opener)
    args = {'ltp':locus_tag}
    request = urllib2.Request(LOCUSTAG_URL, urllib.urlencode(args))
    response = opener.open(request)
    html = response.read() 
    if 'is available at the moment' in html:
        return True
    else:
        return False
    
### ---------------------------------------------------------------------------
def doSubmit(contact, organism_name='', strain='', locus_tag='TEST', dna_source=''):
    opener = urllib2.build_opener(urllib2.HTTPCookieProcessor())
    urllib2.install_opener(opener)
    args = {'max_center_num':'1',
            'max_organ_num':'1',
            'max_contact_num':'1',
            'max_chromo_num':'0',
            'nextpage':'0',
            'mode':'subm',
            'type_prefix':'single_tag',
            'fname0':CONTACTS[contact]['firstname'],
            'lname0':CONTACTS[contact]['lastname'],
            'email0':CONTACTS[contact]['email'],
            'suburl':'http://',
            'suborgan0':'The Wellcome Trust Sanger Institute',
            'suburl0':'http://www.sanger.ac.uk',
            'sequrl':'http://',
            'seqcenter0':'The Wellcome Trust Sanger Institute',
            'sequrl0':'http://www.sanger.ac.uk',
            'organism':organism_name,
            'strain':strain,
            'tag_prefix':locus_tag,
            'dnasource':dna_source,
            }
    request = urllib2.Request(PROJECT_URL, urllib.urlencode(args))
    response = opener.open(request)
    print "Genome project for %s %s submitted. Please check your email!" % (organism_name, strain)
    
### ---------------------------------------------------------------------------
def doRun(list, contact, submit=False):
    print "Reading input file %s" % list
    processed_lines = 0
    for line in open(list, "r"):
        if line[0] == '!':
            continue
        if not line.count('||') == 3:
            continue
        # ! organism_name||strain||locus_tag||dna_source
        processed_lines += 1 
        line = line.strip()
        values = line.split('||')
        organism_name = values[0]
        strain = values[1]
        locus_tag = values[2]
        if values[3] == 'GHP':
            dna_source = 'Gut Health Programme, Rowett Institute of Nutrition and Health, University of Aberdeen. http://www.rowett.ac.uk/divisions/ghp/'
        elif values[3] == 'INRA':
            dna_source = 'INRA Clermont-Ferrand-Theix. http://www.clermont.inra.fr/'
        elif values[3] == 'DSMZ':
            dna_source = 'Deutsche Sammlung von Mikroorganismen und Zellkulturen. GmbH http://www.dsmz.de/'
        elif values[3] == 'NCTC':
            dna_source = "Health Protection Agency's National Collection of Type Cultures. http://www.hpacultures.org.uk/"
        else:
            dna_source = values[3]
        
        print "Processing %s..." % organism_name
        if isLocusTagAvailable(locus_tag):
            print "Locus tag %s is available" % locus_tag
            if submit:
                doSubmit(contact=contact, organism_name=organism_name, strain=strain, locus_tag=locus_tag, dna_source=dna_source)
        else:
            print "Locus tag %s is not available, please choose another one!" % locus_tag
            
        if not submit:
            print "Only checking for availabilities of locus_tag, use '--submit' for submitting data"
            
    print "%s processed lines from input file" % processed_lines

### ---------------------------------------------------------------------------
def main():
    from optparse import OptionParser

    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", "--list", metavar="FILE", help="FILE containing the list of all organism names, its compulsory associated information (organism_name||strain||locus_tag||dna_source)", action="store", type="string", dest="list")
    parser.add_option("-c", "--contact", metavar="CONTACT", help="name of the contact from %s" % CONTACTS.keys(), action="store", choices=CONTACTS.keys(), dest="contact")
    parser.add_option("--submit", help="To submit data, not only checking locus_tag", action="store_true", dest="submit")
    
    (options, args) = parser.parse_args()

    # Print help if no argument given
    if not (options.list and options.contact):
        parser.print_help()
        sys.exit()
    
    # Check proxy settings
    if not os.environ.has_key('http_proxy'):
        print 'You may need to set your http_proxy env variable (e.g. setenv http_proxy "http://wwwcache.sanger.ac.uk:3128")'
    
    # Check input file list exists
    list = options.list
    if not os.path.exists(list):
        print "file %s does not exist" % list
        sys.exit(1)    
    
    doRun(list, options.contact, options.submit)


if __name__ == '__main__':
    main()