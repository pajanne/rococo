'''
Created on Jul 8, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

if __name__ == '__main__':
    import urllib2
    #passmgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
    #passmgr.add_password(None, 'http://proxy:8080', 'username', 'password')
    #authinfo = urllib2.ProxyBasicAuthHandler(passmgr)
    #proxy_support = urllib2.ProxyHandler({"http" : "http://proxy:8080"})
    
    #opener = urllib2.build_opener(proxy_support, authinfo)
    #urllib2.install_opener(opener)
    f = urllib2.urlopen('http://www.google.com')
    buf = f.read()
    print buf
    f.close()

