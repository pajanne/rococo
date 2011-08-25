'''
Created on Feb 11, 2010
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
### GenDB constants
### ---------------------------------------------------------------------------
BASE_URL = 'https://www.cebitec.uni-bielefeld.de/groups/brf/software/gendb-2.2/'
LOGIN_URL = 'https://www.CeBiTec.Uni-Bielefeld.DE/groups/brf/software/gendb-2.2/cgi-bin/login.cgi?'
USER = 'Anne_Pajon'
PASSWORD = 'sangerteam81'

### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------
def logger():
    print "Logging into GenDB"
    proxy_support = urllib2.ProxyHandler({})
    opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(), proxy_support)
    urllib2.install_opener(opener)
    # initial page redirect
    opener.open(BASE_URL)
    # cookie test redirect
    opener.open(BASE_URL + "?cookie_test=1")
    # login  page
    login_args = { 'login':USER, 'pass':PASSWORD }
    login_url = LOGIN_URL + urllib.urlencode(login_args)
    opener.open(login_url)
    # project page
    project_args = {'project':'GenDB_MetaHIT', 'width':'800', 'height':'530'}
    project_url = LOGIN_URL + urllib.urlencode(project_args)
    opener.open(project_url)
    print "Logged into %s as %s" % (BASE_URL, USER)
    return opener

### ---------------------------------------------------------------------------
def submitter(opener, filename, organism_name):
    print "Submitting %s data to GenDB" % organism_name
    

### ---------------------------------------------------------------------------
def doSubmit(common_name, filename):
    opener = logger()
    #submitter(opener, filename, common_name)

### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", metavar="FILE", help="FILE containing the list of all organism common names and its associated sequence file", action="store", type="string", dest="list")
    
    (options, args) = parser.parse_args()

    # Print help if no argument given
    if util.printHelp(options):
        parser.print_help()
        sys.exit()
        
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
        #util.checkFile(input_file)
        doSubmit(common_name, input_file)

if __name__ == '__main__':
    main()