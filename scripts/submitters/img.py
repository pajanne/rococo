'''
Created on Jan 20, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import sys
from optparse import OptionParser
from genepy import util
import urllib, urllib2
import http

### ---------------------------------------------------------------------------
### IMG constants
### ---------------------------------------------------------------------------
BASE_URL = 'http://merced.jgi-psf.org/cgi-bin/img_er_submit/main.cgi?'
USER = 'Anne.Pajon'
PASSWORD = 'nojaP.ennA'

### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------
def logger():
    print "Logging into IMG"
    opener = urllib2.build_opener(urllib2.HTTPCookieProcessor())
    urllib2.install_opener(opener)
    login_query_args = {'username':USER, 'pwd':PASSWORD, '_section_login':'Login' }
    login_encoded_args = urllib.urlencode(login_query_args)
    login_url = BASE_URL + login_encoded_args
    opener.open(login_url)
    return opener

### ---------------------------------------------------------------------------
def submitter(opener, filename, common_name, id):
    print "Submitting %s data to IMG" % common_name
    
    # page 1
    page1_query_args = {'section': 'ERSubmission', 'page':'showERPage'}
    page1_encoded_args = urllib.urlencode(page1_query_args)
    opener.open(BASE_URL, page1_encoded_args)
    
    # page 2
    page2_query_args = {'p_type': 'G', 'database': 'IMG ER', 'search_sp:gold_stamp_id':id, '_section_ERSubmission:selectSubProj':'Search Projects'}
    page2_encoded_args = urllib.urlencode(page2_query_args)
    page2 = opener.open(BASE_URL, page2_encoded_args)
       
    # page 3
    page3_query_args = {'project_oid': ''}
    # extract project_oid from previous page
    page3_query_args = http.extractValues(page2, page3_query_args)
    page3_query_args['database'] = 'IMG ER' 
    page3_query_args['_section_ERSubmission:submitProject'] = 'Select Project' 
    page3_encoded_args = urllib.urlencode(page3_query_args)
    opener.open(BASE_URL, page3_encoded_args)
    
    project_oid = page3_query_args['project_oid']
    
    # page 4 - MultipartForm
    form = http.MultiPartForm()
    form.add_field('database', 'IMG ER')
    form.add_field('project_oid', project_oid)
    form.add_file('genbank_file', filename, open(filename, "r"))
    form.add_field('img_ec_flag', 'Yes')
    form.add_field('gene_calling_flag', 'GeneMark')
    form.add_field('img_product_flag', 'Yes')
    form.add_field('is_img_public', 'No')
    form.add_field('seq_status', 'Draft')
    form.add_field('species_code', common_name)
    form.add_field('_section_ERSubmission:checkSubmission', 'Submit')
    request = urllib2.Request(BASE_URL)
    body = str(form)
    request.add_header('Content-type', form.get_content_type())
    request.add_header('Content-length', len(body))
    request.add_data(body)
    page4 = opener.open(request)
    
    for line in page4.readlines():
        if '(submission ID: ' in line:
            submission_id = line.split('(submission ID: ')[1].split(')')[0]
    print "Submitted to IMG with submission id %s" % submission_id
    return submission_id
    
### ---------------------------------------------------------------------------
def fetcher():
    pass

### ---------------------------------------------------------------------------
def doSubmit(common_name, filename, id):
    opener = logger()
    submitter(opener, filename, common_name, id)

### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-o", metavar="NAME", help="organism common name", action="store", type="string", dest="name")
    parser.add_option("-i", metavar="FILE", help="input organism sequence file in FASTA format", action="store", type="string", dest="input")
    parser.add_option("-p", metavar="ID", help="IMG project ID (GOLD Stamp ID)", action="store", type="string", dest="id")
    parser.add_option("-l", metavar="FILE", help="FILE containing the list of all organism common names, its associated sequence file and IMG project ID", action="store", type="string", dest="list")
    
    (options, args) = parser.parse_args()

    # Print help if no argument given
    if util.printHelp(options):
        parser.print_help()
        sys.exit()
        
    # Get and check input arguments
    if options.list:
        # Read organism common name and related fasta sequence file
        list_file = options.list
        util.checkFile(list_file)
        for line in open(list_file, "r"):
            if line[0] == '!':
                continue
            if line.count('||') < 1:
                continue
            # ! common_name||sequence_file
            line = line.strip()
            values = line.split('||')
            common_name = values[0]
            input_file = values[1]
            id = values[2]
            util.checkFile(input_file)
            doSubmit(common_name, input_file, id)
    else:
        common_name = options.name
        input_file = options.input
        id = options.id
        util.checkFile(input_file)
        doSubmit(common_name, input_file, id)


if __name__ == '__main__':
    main()