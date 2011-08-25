'''
Created on Feb 3, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import sys
from optparse import OptionParser
from genepy import util
import urllib, urllib2

### ---------------------------------------------------------------------------
### constants
### ---------------------------------------------------------------------------
SUBMISSION_URL = 'https://www.ebi.ac.uk/embl/genomes/submission/app/submission/new-submission.jsf'
PROJECT_URL = 'https://www.ebi.ac.uk/embl/genomes/submission/app/project/project.jsf'
LOGIN_URL = 'https://www.ebi.ac.uk/embl/genomes/submission/login.jsf'
USER = 'ap12@sanger.ac.uk'
PASSWORD = 'sangerteam81'

### ---------------------------------------------------------------------------
### Main methods 
### ---------------------------------------------------------------------------
def doSubmit(organism_name='', strain='', locus_tag='TEST', seq_size='', seq_depth='', dna_source='', description='', submit=False):
    
    #urllib2.HTTPCookieProcessor()
    opener = urllib2.build_opener()
    urllib2.install_opener(opener)
    headers = {'Host':'www.ebi.ac.uk',
               'User-Agent':'Mozilla/5.0 (Macintosh; U; Intel Mac OS X 10.5; en-GB; rv:1.9.1) Gecko/20090624 (CK-WTSI) Firefox/3.5', 
               'Accept':'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
               'Accept-Language':'en-gb,en;q=0.5',
               'Accept-Encoding':'gzip,deflate',
               'Keep-Alive':'300',
               'Accept-Charset':'ISO-8859-1,utf-8;q=0.7,*;q=0.7',
               'Connection':'keep-alive',
               'Referer':'https://www.ebi.ac.uk/embl/genomes/submission/login.jsf'}
    # ----------
    # Login does not work
    #login_args = {'login-user':USER, 'login-password':PASSWORD, 
    #              'oracle.adf.faces.FORM': 'mainForm',
    #              'oracle.adf.faces.STATE_TOKEN':'1',
    #              'source':'submit-button',
    #              'event':''}
    #login_encoded_args = urllib.urlencode(login_args)
    # extract my session id by loading a page from the site
    #set_cookie = urllib2.urlopen(LOGIN_URL, login_encoded_args, headers).headers.getheader("Set-Cookie")
    #session_id = set_cookie[set_cookie.index("=")+1:set_cookie.index(";")]
    # ----------
    # Copy/paste your Cookie header from firebug
    headers['Cookie'] = 'JSESSIONID=51F7E7D5B966A0AAD627B1BE71C43FC6; __utmz=222765775.1258459105.1.1.utmccn=(organic)|utmcsr=google|utmctr=ebi+|utmcmd=organic; oracle.uix=0^^GMT-0:00; __utma=222765775.769608750.1258459105.1265048083.1265106200.6' 
    #launcher-active=false&submission-options=submission-project&oracle.adf.faces.FORM=mainForm&oracle.adf.faces.STATE_TOKEN=e&source=submit-button&event=
    project_args = {'launcher-active':'false',
                    'submission-options':'submission-project',
                    'oracle.adf.faces.FORM':'mainForm',
                    'oracle.adf.faces.STATE_TOKEN':'e',
                    'source':'submit-button',
                    'event':''}
    project_encoded_args = urllib.urlencode(project_args)
    request = urllib2.Request(SUBMISSION_URL, project_encoded_args, headers=headers)
    response = opener.open(request)
    next_url = response.geturl()
    response.close()
    
    #_id35=0&oracle.adf.faces.FORM=mainForm&oracle.adf.faces.STATE_TOKEN=f&source=submit-button&event=
    project_args = {'_id35':'0',
                    'oracle.adf.faces.FORM':'mainForm',
                    'oracle.adf.faces.STATE_TOKEN':'f',
                    'source':'submit-button',
                    'event':''}
    project_encoded_args = urllib.urlencode(project_args)
    request = urllib2.Request(next_url, project_encoded_args, headers=headers)
    response = opener.open(request)
    next_url = response.geturl()
    response.close()
    
    # Check locus_tag prefix
    #locus-form%3Alocustag-prefix=TEST&contact-form%3Acontact-table%3A_us=0&contact-form%3Acontact-table%3A_us=1&contact-form%3Acontact-table%3A1%3A_id55=Keith&contact-form%3Acontact-table%3A1%3A_id56=Turner&contact-form%3Acontact-table%3A1%3A_id57=akt%40sanger.ac.uk&contact-form%3Acontact-table%3ArangeStart=0&submitter-form%3Asubmitter-table%3A_us=0&submitter-form%3Asubmitter-table%3A0%3A_id61=The+Wellcome+Trust+Sanger+Institute&submitter-form%3Asubmitter-table%3A0%3A_id63=http%3A%2F%2Fwww.sanger.ac.uk&submitter-form%3Asubmitter-table%3ArangeStart=0&center-form%3Acenter-table%3A_us=0&center-form%3Acenter-table%3A0%3A_id67=The+Wellcome+Trust+Sanger+Institute&center-form%3Acenter-table%3A0%3A_id69=http%3A%2F%2Fwww.sanger.ac.uk&center-form%3Acenter-table%3ArangeStart=0&consortium-form%3Aconsortium-name=metaHIT&consortium-form%3Aconsortium-url=http%3A%2F%2Fwww.metahit.eu%2F&replicon-form%3Areplicon-table%3ArangeStart=0&project-form%3Aorganism-name=&project-form%3Astrain-name=&project-form%3A_id95=&project-form%3A_id96=1&project-form%3A_id102=&project-form%3A_id103=http%3A%2F%2Fwww.sanger.ac.uk%2Fpathogens%2Fmetahit%2F&project-form%3A_id104=ftp%3A%2F%2Fftp.sanger.ac.uk%2Fpub%2Fpathogens%2Fmetahit%2F&project-form%3Aseq-size=&project-form%3A_id107=&oracle.adf.faces.FORM=mainForm&oracle.adf.faces.STATE_TOKEN=1d&source=locus-form%3Acheck-prefix-button&event=&contact-form%3Acontact-table%3A_sm=&partialTargets=&partial=&submitter-form%3Asubmitter-table%3A_sm=&center-form%3Acenter-table%3A_sm=
    #locus-form%3Alocustag-prefix=CLS &contact-form%3Acontact-table%3A_us=0&contact-form%3Acontact-table%3A_s=1&contact-form%3Acontact-table%3A_us=1&contact-form%3Acontact-table%3A1%3A_id55=Keith&contact-form%3Acontact-table%3A1%3A_id56=Turner&contact-form%3Acontact-table%3A1%3A_id57=akt%40sanger.ac.uk&contact-form%3Acontact-table%3ArangeStart=0&submitter-form%3Asubmitter-table%3A_s=0&submitter-form%3Asubmitter-table%3A_us=0&submitter-form%3Asubmitter-table%3A0%3A_id61=The+Wellcome+Trust+Sanger+Institute&submitter-form%3Asubmitter-table%3A0%3A_id63=http%3A%2F%2Fwww.sanger.ac.uk&submitter-form%3Asubmitter-table%3ArangeStart=0&center-form%3Acenter-table%3A_s=0&center-form%3Acenter-table%3A_us=0&center-form%3Acenter-table%3A0%3A_id67=The+Wellcome+Trust+Sanger+Institute&center-form%3Acenter-table%3A0%3A_id69=http%3A%2F%2Fwww.sanger.ac.uk&center-form%3Acenter-table%3ArangeStart=0&consortium-form%3Aconsortium-name=metaHIT&consortium-form%3Aconsortium-url=http%3A%2F%2Fwww.metahit.eu%2F&replicon-form%3Areplicon-table%3ArangeStart=0&project-form%3Aorganism-name=&project-form%3Astrain-name=K10&project-form%3A_id95=Gut+Health+Programme%2C+Rowett+Institute+of+Nutrition+and+Health%2C+University+of+Aberdeen.+http%3A%2F%2Fwww.rowett.ac.uk%2Fdivisions%2Fghp%2F&project-form%3A_id96=1&project-form%3A_id102=20x&project-form%3A_id103=http%3A%2F%2Fwww.sanger.ac.uk%2Fpathogens%2Fmetahit%2F&project-form%3A_id104=ftp%3A%2F%2Fftp.sanger.ac.uk%2Fpub%2Fpathogens%2Fmetahit%2F&project-form%3Aseq-size=3.8&project-form%3A_id107=Shown+to+harbor+several+novel+antibiotic+resistance+genes.+%28Kazimierczak+et+al+2006%2C+Scott+et+al+2008%29.&oracle.adf.faces.FORM=mainForm&oracle.adf.faces.STATE_TOKEN=1b&source=submit-button&event=&contact-form%3Acontact-table%3A_sm=&partialTargets=&partial=&submitter-form%3Asubmitter-table%3A_sm=&center-form%3Acenter-table%3A_sm=
    
    project_args = {'locus-form:locustag-prefix':locus_tag,
                    'contact-form:contact-table:_us':'0',
                    'contact-form:contact-table:_us':'1',
                    'contact-form:contact-table:1:_id55':'Keith',
                    'contact-form:contact-table:1:_id56':'Turner',
                    'contact-form:contact-table:1:_id57':'akt@sanger.ac.uk',
                    'contact-form:contact-table:rangeStart':'0',
                    'submitter-form:submitter-table:_us':'0',
                    'submitter-form:submitter-table:0:_id61':'The Wellcome Trust Sanger Institute',
                    'submitter-form:submitter-table:0:_id63':'http://www.sanger.ac.uk',
                    'submitter-form:submitter-table:rangeStart':'0',
                    'center-form:center-table:_us':'0',
                    'center-form:center-table:0:_id67':'The Wellcome Trust Sanger Institute',
                    'center-form:center-table:0:_id69':'http://www.sanger.ac.uk',
                    'center-form:center-table:rangeStart':'0',
                    'consortium-form:consortium-name':'metaHIT',
                    'consortium-form:consortium-url':'http://www.metahit.eu/',
                    'replicon-form:replicon-table:rangeStart':'0',
                    'project-form:organism-name':organism_name,
                    'project-form:strain-name':strain,
                    'project-form:_id95':dna_source,
                    'project-form:_id96':'1',
                    'project-form:_id102':seq_depth,
                    'project-form:_id103':'http://www.sanger.ac.uk/pathogens/metahit/',
                    'project-form:_id104':'ftp://ftp.sanger.ac.uk/pub/pathogens/metahit/',
                    'project-form:seq-size':seq_size,
                    'project-form:_id107':description,
                    'oracle.adf.faces.FORM':'mainForm',
                    'oracle.adf.faces.STATE_TOKEN':'3f',
                    'source':'locus-form:check-prefix-button',
                    'event':'',
                    'contact-form:contact-table:_s':'1',
                    'partialTargets':'',
                    'partial':'',
                    'submitter-form:submitter-table:_s':'0',
                    'center-form:center-table:_s':'0'}
    project_encoded_args = urllib.urlencode(project_args)
    request = urllib2.Request(PROJECT_URL, project_encoded_args, headers=headers)
    response = opener.open(request)
    print "----------------------------------------"
    locus_page = response.read()
    response.close()
    if 'Provided locus tag prefix is available' in locus_page:
        print "Locus tag %s is available" % locus_tag
        if submit:
            print "Submitting Genome Project to EBI/EMBL Nucleotide Sequence Database"
            # Submit genome project
            #locus-form%3Alocustag-prefix=TEST&contact-form%3Acontact-table%3A_us=0&contact-form%3Acontact-table%3A_us=1&contact-form%3Acontact-table%3A1%3A_id55=Keith&contact-form%3Acontact-table%3A1%3A_id56=Turner&contact-form%3Acontact-table%3A1%3A_id57=akt%40sanger.ac.uk&contact-form%3Acontact-table%3ArangeStart=0&submitter-form%3Asubmitter-table%3A_us=0&submitter-form%3Asubmitter-table%3A0%3A_id61=The+Wellcome+Trust+Sanger+Institute&submitter-form%3Asubmitter-table%3A0%3A_id63=http%3A%2F%2Fwww.sanger.ac.uk&submitter-form%3Asubmitter-table%3ArangeStart=0&center-form%3Acenter-table%3A_us=0&center-form%3Acenter-table%3A0%3A_id67=The+Wellcome+Trust+Sanger+Institute&center-form%3Acenter-table%3A0%3A_id69=http%3A%2F%2Fwww.sanger.ac.uk&center-form%3Acenter-table%3ArangeStart=0&consortium-form%3Aconsortium-name=metaHIT&consortium-form%3Aconsortium-url=http%3A%2F%2Fwww.metahit.eu%2F&replicon-form%3Areplicon-table%3ArangeStart=0&project-form%3Aorganism-name=&project-form%3Astrain-name=&project-form%3A_id95=&project-form%3A_id96=1&project-form%3A_id102=&project-form%3A_id103=http%3A%2F%2Fwww.sanger.ac.uk%2Fpathogens%2Fmetahit%2F&project-form%3A_id104=ftp%3A%2F%2Fftp.sanger.ac.uk%2Fpub%2Fpathogens%2Fmetahit%2F&project-form%3Aseq-size=&project-form%3A_id107=&oracle.adf.faces.FORM=mainForm&oracle.adf.faces.STATE_TOKEN=1e&source=submit-button&event=&contact-form%3Acontact-table%3A_sm=&partialTargets=&partial=&submitter-form%3Asubmitter-table%3A_sm=&center-form%3Acenter-table%3A_sm=
            project_args['source'] = 'submit-button'
            project_encoded_args = urllib.urlencode(project_args)
            request = urllib2.Request(PROJECT_URL, project_encoded_args, headers=headers)
            response = opener.open(request)
            page = response.read()
            if '>Error<' in page:
                print "Errors when submitting genome project for %s %s" % (organism_name, strain)
                print "----------------------------------------"
                print page
                print "----------------------------------------"
            else:
                print "Submitting genome project for %s %s successful!" % (organism_name, strain)
            
            print "----------------------------------------"
            response.close()
        else:
            print "Only checking for availabilities of locus_tag"
    
    else:
        print "Locus tag %s is not available, please choose another one!" % locus_tag
        print "----------------------------------------"
        print locus_page
        print "----------------------------------------"
    

### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", metavar="FILE", help="FILE containing the list of all organism common names, its associated information", action="store", type="string", dest="list")
    parser.add_option("--submit", help="To submit data, not only checking locus_tag", action="store_true", dest="submit")
    
    (options, args) = parser.parse_args()

    # Print help if no argument given
    if util.printHelp(options):
        parser.print_help()
        sys.exit()
        
    # Get and check input arguments
    # Read organism common name and related fasta sequence file
    list_file = options.list
    util.checkFile(list_file)
    for line in open(list_file, "r"):
        if line[0] == '!':
            continue
        if line.count('||') < 1:
            continue
        # ! organism_name||strain||locus_tag||seq_size||seq_depth||dna_source||description
        line = line.strip()
        values = line.split('||')
        organism_name = values[0]
        strain = values[1]
        locus_tag = values[2]
        seq_size = values[3]
        seq_depth = values[4]
        if values[5] == 'GHP':
            dna_source = 'Gut Health Programme, Rowett Institute of Nutrition and Health, University of Aberdeen. http://www.rowett.ac.uk/divisions/ghp/'
        elif values[5] == 'INRA':
            dna_source = 'INRA Clermont-Ferrand-Theix. http://www.clermont.inra.fr/'
        elif values[5] == 'DSMZ':
            dna_source = 'Deutsche Sammlung von Mikroorganismen und Zellkulturen. GmbH http://www.dsmz.de/'
        elif values[5] == 'NCTC':
            dna_source = "Health Protection Agency's National Collection of Type Cultures. http://www.hpacultures.org.uk/"
        else:
            print "DNA source %s not found! Please provide details..." % values[5]
            continue
        
        #print dna_source
        description = values[6]
        doSubmit(organism_name=organism_name, strain=strain, locus_tag=locus_tag, 
                 seq_size=seq_size, seq_depth=seq_depth, dna_source=dna_source, 
                 description=description, submit=options.submit)


if __name__ == '__main__':
    main()