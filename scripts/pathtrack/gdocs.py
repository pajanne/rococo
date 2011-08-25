'''
Created on May 26, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import gdata.spreadsheet.service

def getValues(doc='Metahit_454_Projects', email='pathpipe@gmail.com', password='pathpipeE204'):
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
    values = []
    for row in rows:
        if not row.custom['study'].text is None:
            if row.custom['study'].text[0] == '!':
                continue
        # create a list of project||genus||species-subspecies||strain||sample||library||run||paired||insert_size||path_to_sff||path_to_assembly
        value = "%s||%s||%s||%s||%s||%s||%s||%s||%s||%s||%s" % (row.custom['study'].text,
                                                        row.custom['genus'].text,
                                                        row.custom['species-subspecies'].text,
                                                        row.custom['strain'].text,
                                                        row.custom['sample'].text,
                                                        row.custom['library'].text,
                                                        row.custom['run'].text,
                                                        row.custom['paired'].text,
                                                        row.custom['insertsize'].text,
                                                        row.custom['pathtosff'].text,
                                                        row.custom['pathtoassembly'].text)
        values.append(value)
    return values
        
if __name__ == '__main__': 
    values = getValues()
    for value in values:
        print value
    
