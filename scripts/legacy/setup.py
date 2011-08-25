'''
Created on Nov 23, 2009
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import sys, os
import ConfigParser
from ropy.query import ConnectionFactory
from ropy.log import LogSetup

### ---------------------------------------------------------------------------
### Load Configuration file
### ---------------------------------------------------------------------------
try:
    _config_file = "%s/conf/config.ini" % os.path.dirname(__file__)
    config = ConfigParser.ConfigParser()
    config.read(_config_file)
except Exception, e:
    print ("Could not read this file " + _config_file)
    sys.exit(1)

### ---------------------------------------------------------------------------
### Database Setup
### ---------------------------------------------------------------------------
host=config.get('Connection', 'host')
database=config.get('Connection', 'database')
user=config.get('Connection', 'user')
password=config.get('Connection', 'password')

connectionFactory = ConnectionFactory(host, database, user, password)

### ---------------------------------------------------------------------------
### Log Setup
### ---------------------------------------------------------------------------
try:
    logsetup = LogSetup()
    logsetup.logname = "genepy"
    logsetup.logpath = config.get('Logging', 'path')
    logsetup.setupLogging()
    logger = logsetup.logger
except Exception, e:
    print e
    print "Could not set up the logging, is the path set in the configuration file correct?"
    sys.exit(1)


### ---------------------------------------------------------------------------
def main():
    pass

if __name__ == '__main__':
    main()
