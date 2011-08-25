'''
Created on Jan 18, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import urllib, urllib2
import sys
from optparse import OptionParser
from genepy import util
import http

### ---------------------------------------------------------------------------
### RAST constants
### ---------------------------------------------------------------------------
BASE_URL = 'http://rast.nmpdr.org/rast.cgi?'
USER = 'ap12'
PASSWORD = 'hITt2Csz'

### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------
def checkValidInput(input_file, common_name):
    """
    Check if the input fasta sequence file is of correct format.
    RAST re-arrange the scaffolds if a splitted sequences is submitted
    Run EMBOSS union before if more than on '>' is found
    """
    util.checkFile(input_file)
    cmd = "grep '>' %s | wc -l" % input_file
    result = util.runProcess(cmd)
    if int(result) > 1:
        new_input_file = "%s.fna" % common_name
        util.checkSoft("union")
        util.checkSoft("descseq")
        cmd_union = "union -sequence %s -stdout Yes -auto Yes | descseq -filter Yes -name '%s' -auto Yes > %s" % (input_file, common_name, new_input_file)
        util.runProcess(cmd_union)
        return new_input_file
    else:
        return input_file

def logger():
    """
    urllib2 std lib doc: http://docs.python.org/library/urllib2.html
    urllib2 cookbook: http://personalpages.tds.net/~kent37/kk/00010.html
    """
    print "Logging into RAST"
    opener = urllib2.build_opener(urllib2.HTTPCookieProcessor())
    urllib2.install_opener(opener)
    login_query_args = { 'login':USER, 'password':PASSWORD, 'page':'Login', 'action':'perform_login' }
    login_encoded_args = urllib.urlencode(login_query_args)
    login_url = BASE_URL + login_encoded_args
    opener.open(login_url)
    return opener

### ---------------------------------------------------------------------------
def fetcher(opener, job_id):
    print "Fetching filename of job %s" % job_id
    response = opener.open("%spage=JobDetails&job=%s" % (BASE_URL, job_id))
    for line in response.readlines():
        if "\">EMBL</option>" in line:
            filename = line.split('"')[1]
    print filename

    print "Fetching data"
    fetch_query_args = { 'file':filename, 'page':'DownloadFile', 'job':job_id }
    fetch_encoded_args = urllib.urlencode(fetch_query_args)
    response = opener.open(BASE_URL + fetch_encoded_args)
    return response.read()

### ---------------------------------------------------------------------------
def submitter(opener, filename, organism_name):
    print "Submitting %s data to RAST" % organism_name
    
    # page1 - MultipartForm
    form = http.MultiPartForm()
    form.add_field('page', 'UploadGenome')
    form.add_file('upload', filename, open(filename, "r"))
    request = urllib2.Request(BASE_URL)
    body = str(form)
    request.add_header('Content-type', form.get_content_type())
    request.add_header('Content-length', len(body))
    request.add_data(body)
    page1 = opener.open(request)
    
    # page2 - simple form
    page2_query_args = {'upload_file':'', 'contig_count':'', 'upload_type':'', 'bp_count':'', 'gc_content':'', 'ambig_count':'', 'upload_check':''}
    # These have to be parsed out of the previous page    
    page2_query_args = http.extractValues (page1, page2_query_args)
    page2_query_args['organism_name'] = organism_name
    page2_query_args['genetic_code'] = '11'
    page2_query_args['page'] = 'UploadGenome'
    page2_query_args['domain'] = 'Bacteria'
    page2_query_args['taxonomy_id'] = ''
    page2_query_args['lineage'] = ''
    page2_query_args['laststep'] = 'Use this data and go to step 3'
    page2_encoded_args = urllib.urlencode(page2_query_args)
    page2 = opener.open(BASE_URL, page2_encoded_args)
    
    # page3 - simple form
    page3_query_args = {'organism_name':'', 'upload_file':'', 'contig_count':'', 'upload_type':'', 'gc_content':'', 'upload_check':'', 'laststep':'', 'bp_count':'', 'ambig_count':''}
    # These have to be parsed out of the previous page    
    page3_query_args = http.extractValues(page2, page3_query_args)
    page3_query_args['genetic_code'] = '11'
    page3_query_args['page'] = 'UploadGenome'
    page3_query_args['domain'] = 'Bacteria'
    page3_query_args['taxonomy_id'] = ''
    page3_query_args['lineage'] = ''
    page3_query_args['sequencing_method'] = 'Sanger'
    page3_query_args['coverage'] = 'unknown'
    page3_query_args['contigs'] = 'unknown'
    page3_query_args['average_read_length'] = ''
    page3_query_args['gene_caller'] = 'rast'
    page3_query_args['fix_errors'] = '1'
    page3_query_args['backfill_gaps'] = '1'
    page3_query_args['finish'] = 'Finish the upload'
    page3_encoded_args = urllib.urlencode(page3_query_args)
    page3 = opener.open(BASE_URL, page3_encoded_args)
    
    for line in page3.readlines():
        if "href='?page=JobDetails&job=" in line:
            id = line.split("href='?page=JobDetails&job=")[1].split("'")[0]
    print "Submitted to RAST with job id %s" % id
    return id

### ---------------------------------------------------------------------------
def doSubmit(common_name, filename):
    opener = logger()
    submitter(opener, filename, common_name)
    
def doFetch(common_name, job_id):
    opener = logger()
    embl = fetcher(opener, job_id)
    filename = "%s.rast.embl" % common_name
    f = open(filename, 'w')
    f.write(embl)
    

### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-o", metavar="NAME", help="organism common name", action="store", type="string", dest="name")
    parser.add_option("-i", metavar="FILE", help="input organism sequence file in FASTA format", action="store", type="string", dest="input")
    parser.add_option("-j", metavar="ID", help="input job ID to fetch results", action="store", type="string", dest="jobid")
    parser.add_option("-l", metavar="FILE", help="FILE containing the list of all organism common names and its associated sequence file", action="store", type="string", dest="list")
    parser.add_option("--fetch", help="To fetch results, job id must be provided", action="store_true", dest="fetch")
    
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
            if options.fetch:
                job_id = values[2]
                doFetch(common_name, job_id)
            else:
                input_file = checkValidInput(values[1], common_name)           
                doSubmit(common_name, input_file)
    else:
        common_name = options.name
        if options.fetch:
                job_id = options.jobid
                doFetch(common_name, job_id)
        else:
            input_file = checkValidInput(options.input, common_name)           
            doSubmit(common_name, input_file)

if __name__ == '__main__':
    main()