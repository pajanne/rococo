'''
Created on Jan 20, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import mimetools, mimetypes
import itertools

### ---------------------------------------------------------------------------
### Utility methods
### ---------------------------------------------------------------------------
def extractValues(page, args):
    """
    Extract values associated to input dictionary keys from a web page 
    and if found store them back into args dictionary
    """
    for line in page.readlines():
        for key in args.keys():
            search_key = "name='%s'" % key
            if search_key in line:
                args[key] = line.split("value='")[1].split("'")[0]
    return args

### ---------------------------------------------------------------------------
### Class MultiPartForm
### ---------------------------------------------------------------------------
class MultiPartForm(object):
    """
    Accumulate the data to be used when posting a form.
    From http://www.doughellmann.com/PyMOTW/urllib2/
    """

    def __init__(self):
        self.form_fields = []
        self.files = []
        self.boundary = mimetools.choose_boundary()
        return
    
    def get_content_type(self):
        return 'multipart/form-data; boundary=%s' % self.boundary

    def add_field(self, name, value):
        """Add a simple field to the form data."""
        self.form_fields.append((name, value))
        return

    def add_file(self, fieldname, filename, fileHandle, mimetype=None):
        """Add a file to be uploaded."""
        body = fileHandle.read()
        if mimetype is None:
            mimetype = mimetypes.guess_type(filename)[0] or 'application/octet-stream'
        self.files.append((fieldname, filename, mimetype, body))
        return
    
    def __str__(self):
        """Return a string representing the form data, including attached files."""
        # Build a list of lists, each containing "lines" of the
        # request.  Each part is separated by a boundary string.
        # Once the list is built, return a string where each
        # line is separated by '\r\n'.  
        parts = []
        part_boundary = '--' + self.boundary
        
        # Add the form fields
        parts.extend(
            [ part_boundary,
              'Content-Disposition: form-data; name="%s"' % name,
              '',
              value,
            ]
            for name, value in self.form_fields
            )
        
        # Add the files to upload
        parts.extend(
            [ part_boundary,
              'Content-Disposition: file; name="%s"; filename="%s"' % \
                 (field_name, filename),
              'Content-Type: %s' % content_type,
              '',
              body,
            ]
            for field_name, filename, content_type, body in self.files
            )
        
        # Flatten the list and add closing boundary marker,
        # then return CR+LF separated data
        flattened = list(itertools.chain(*parts))
        flattened.append('--' + self.boundary + '--')
        flattened.append('')
        return '\r\n'.join(flattened)
