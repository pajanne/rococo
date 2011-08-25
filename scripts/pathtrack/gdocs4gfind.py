'''
Created on Oct 13, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

from optparse import OptionParser
import sys, os
import gdata.spreadsheet.service

### ---------------------------------------------------------------------------
values = {}

GFIND_TEMPLATE = """
root    => '%(root)s/%(common_name)s',
module  => 'PathTrack::GeneFinding',
prefix  => '_',
log	=> '%(root)s/log/%(common_name)s.log',

data => {
    fasta => '%(fasta_file)s',
    common_name => '%(common_name)s',
},

"""

GFUNC_TEMPLATE = """
root    => '%(root)s/%(common_name)s',
module  => 'PathTrack::GeneFunction',
prefix  => '_',
log	=> '%(root)s/log/%(common_name)s.log',

data => {
    embl => '%(root)s/%(common_name)s/GFIND/sequence.embl',
    common_name => '%(common_name)s',
},

"""

### ---------------------------------------------------------------------------
def getValuesFromGoogleDoc(doc='Annotation_Pipeline', email='pathpipe@gmail.com', password='pathpipeE204'):
    # connect
    gd_client = gdata.spreadsheet.service.SpreadsheetsService()
    gd_client.email = email
    gd_client.password = password
    gd_client.source = 'Reading %s from python' % doc
    gd_client.ProgrammaticLogin()

    # query the spreadsheet by name
    query  = gdata.spreadsheet.service.DocumentQuery()
    query['title'] = doc
    query['title-exact'] = 'true'
    feed = gd_client.GetSpreadsheetsFeed(query=query)
    
    # extract the unique spreadsheet and worksheet IDs
    spreadsheet_id = feed.entry[0].id.text.rsplit('/',1)[1]
    feed = gd_client.GetWorksheetsFeed(spreadsheet_id)
    worksheet_id = feed.entry[0].id.text.rsplit('/',1)[1]

    # get list of all rows
    rows = gd_client.GetListFeed(spreadsheet_id, worksheet_id).entry
    for row in rows:
        if row.custom['uniquecommonname'].text is None:
            continue
        # populate dict of values to process
        values[row.custom['uniquecommonname'].text] = row.custom['pathtofastafile'].text


### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-r", metavar="PATH", help="PATH to the root of the hierarchy (e.g. /lustre/scratch103/pathogen/gfind/)", action="store", type="string", dest="root")
    (options, args) = parser.parse_args()
    
    if not (options.root):
        parser.print_help()
        sys.exit()

    # populate dictionnary from Google Document
    getValuesFromGoogleDoc()

    # check root path
    if not os.path.exists(options.root):
        print "%s path do not exist! Please create root path first." % options.root
        sys.exit()
    
    # check log directory exists
    log_dir = "%s/log/" % (options.root)
    if not os.path.exists(log_dir):
        print "%s path do not exist! Please create first." % log_dir
        sys.exit()
    
    # check conf directory exists
    conf_dir = "%s/conf/" % (options.root)
    if not os.path.exists(conf_dir):
        print "%s path do not exist! Please create first." % conf_dir
        sys.exit()

    # open gfind_pipeline.conf
    conf_file = open('%s/conf/gfind_pipeline.conf' % options.root, 'w')
    for common_name in values.keys():
        print "%s\t%s" % (common_name, values[common_name])
        # open gfind conf file
        gfind_conf_filename = '%s/conf/%s_gfind.conf' % (options.root, common_name)
        gfind_conf_file = open(gfind_conf_filename, 'w')
        gfind_conf_file.write(GFIND_TEMPLATE % {'root':options.root,
                                                'common_name':common_name,
                                                'fasta_file':values[common_name]})
        gfind_conf_file.close()
        # open gfunc conf file
        gfunc_conf_filename = '%s/conf/%s_gfunc.conf' % (options.root, common_name)
        gfunc_conf_file = open(gfunc_conf_filename, 'w')
        gfunc_conf_file.write(GFUNC_TEMPLATE % {'root':options.root,
                                                'common_name':common_name})
        gfunc_conf_file.close()
        conf_file.write('GFIND\t%s\n' % gfind_conf_filename)
        conf_file.write('GFUNC\t%s\n' % gfunc_conf_filename)
    conf_file.close()
    
### ---------------------------------------------------------------------------
if __name__ == '__main__':
    main()

    
