'''
Created on Mar 11, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import os
import logging.config
logging.config.fileConfig('%s/conf/logging.conf' % os.path.realpath(os.path.dirname(__file__)))
