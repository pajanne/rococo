#!/usr/bin/env python
# encoding: utf-8
'''
Created on Jan 5, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import os
from genepy import util
import sys
import loaders.db as db
import ropy.util
import ropy.query
from optparse import OptionParser
from ropy.log import LogSetup

### ---------------------------------------------------------------------------
### Logging setup
### ---------------------------------------------------------------------------
logsetup = LogSetup()
logsetup.logname = "genepy-multiloader"
logsetup.logpath = "%s/logs.txt" % os.path.realpath(os.path.dirname(__file__))
logsetup.setupLogging()

logger = logsetup.logger

### ---------------------------------------------------------------------------
### main
### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", "--list", metavar="FILE", help="FILE containing the list of all organisms", action="store", type="string", dest="list")
    parser.add_option("-D", action="store", dest="dbhost")

    (options, args) = parser.parse_args()
    
    # Print help if no argument given
    if util.printHelp(options):
        parser.print_help()
        sys.exit()
    
    # Print command line
    cmdline = "$ python "
    for argv in sys.argv:
        cmdline += argv + " " 
    logger.info(cmdline)
    
    # Print logger file info
    logger.info(logsetup.logpath)
    
    # Setup database connection
    host = ropy.util.getDArg("dbhost", raiseOnEmpty = True)
    database = ropy.util.getDArg("dbname", raiseOnEmpty = True)
    port = ropy.util.getDArg("dbport", raiseOnEmpty = True)
    user = ropy.util.getDArg("dbuser", raiseOnEmpty = True)
    password = ropy.util.getDArg("dbpassword")
    connectionFactory = ropy.query.ConnectionFactory(host, database, user, password, port)
    
    # Read organism list file and load it into the database
    data_path = options.list
    for line in open(data_path, "r"):
        if line[0] == '!':
            continue
        if line.count('||') < 1:
            continue
        # ! Genus||species||strain||taxonId
        line = line.strip()
        list = line.split('||')
        genus = list[0]
        species = list[1].replace('sp.', 'unknown')
        strain = list[2]
        taxonid = list[3]
        
        # Load organism
        chado_species = "%s (%s)" % (species, strain)
        common_name = getCommonName(genus, species, strain)
        abbreviation = common_name
        comment = None
        logger.info(common_name)
        logger.info(db.makeOrganism(connectionFactory, genus, chado_species, abbreviation, common_name, comment))
        
        # Load translation table
        logger.info(db.makeOrganismProp(connectionFactory, genus, chado_species, "genedb_misc", "translationTable", 11))
        
        # Load taxonomy id
        logger.info(db.makeOrganismProp(connectionFactory, genus, chado_species, "genedb_misc", "taxonId", taxonid))
        
        # Load HTML name fields for GeneDB web
        htmlFullName = getHtmlFullName(genus, species, strain)
        logger.info(db.makeOrganismProp(connectionFactory, genus, chado_species, "genedb_misc", "htmlFullName", htmlFullName))
        htmlShortName = getHtmlShortName(genus, species, strain)
        logger.info(db.makeOrganismProp(connectionFactory, genus, chado_species, "genedb_misc", "htmlShortName", htmlShortName))
    

### ---------------------------------------------------------------------------
def getHtmlFullName(genus, species, strain):
    if len(strain) > 1:
        return "<i>%s %s</i> %s" % (genus, species, strain) # <i>Genus species</i> strain
    else:
        return "<i>%s %s</i>" % (genus, species)

def getHtmlShortName(genus, species, strain):
    if len(strain) > 1:
        return "<i>%s. %s</i> %s" % (genus[0], species, strain) # <i>G. species</i> strain
    else:
        return "<i>%s. %s</i>" % (genus[0], species)

def getCommonName(genus, species, strain):
    # Gspecies_strain
    if len(strain) > 1:
        return "%s%s_%s" % (genus[0], species, getNameWithoutSpecialChar(strain))
    else:
        return "%s%s" % (genus[0], species)

def getNameWithoutSpecialChar(name):
    char_to_remove = ['/', '-', '.', ' ', '(', ')']
    for c in char_to_remove:
        name = name.replace(c, '')
    return name

### ---------------------------------------------------------------------------
if __name__ == '__main__':
    main()
