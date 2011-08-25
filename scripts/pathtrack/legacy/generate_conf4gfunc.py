'''
Created on Aug 19, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.


Takes a || delimited file of format:
! common_name||embl_file
and generate GeneFunction configuration files

Usage:
> python ~ap12/genlibpy/genepy/pathtrack/generate_conf4gfunc.py -l /lustre/scratch103/sanger/ap12/metahit_analysis/ORGANISMS.list -r /lustre/scratch103/sanger/ap12/metahit_analysis/

'''

from optparse import OptionParser
import sys, os

TEMPLATE = """
root    => '%(root)s/%(common_name)s',
module  => 'PathTrack::GeneFunction',
prefix  => '_',
log	=> '%(root)s/pipeline.log',

data => {
    embl => '%(embl_file)s',
    common_name => '%(common_name)s',
},

"""

### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", "--list", metavar="FILE", help="FILE containing the list of all organisms to process", action="store", type="string", dest="list")
    parser.add_option("-r", "--root", metavar="PATH", help="PATH to the root of the hierarchy", action="store", type="string", dest="root")
    
    (options, args) = parser.parse_args()

    if not (options.list and options.root):
        parser.print_help()
        sys.exit()

    # check root path
    if not os.path.exists(options.root):
        print "%s path do not exist" % options.root
        print "Create root path first."
        sys.exit()
    
    # check conf directory exists
    conf_dir = "%s/conf/" % options.root
    if not os.path.exists(conf_dir):
        print "%s path do not exist" % conf_dir
        print "Create conf directory first."
        sys.exit()
        
    # open gfunc_pipeline.conf
    top_conf_file = open('%s/conf/gfunc_pipeline.conf' % options.root, 'w')

    # check input assembly file and read it
    if not os.path.exists(options.list):
        print "%s file do not exist" % options.list
        sys.exit()
        
    lines = open(options.list, "r").readlines()
    for line in lines:
        if line[0] == '!':
            continue
        if not line.count('||') == 1:
            continue
        line = line.strip()
        values = line.split('||')
        common_name = values[0]
        embl_file = values[1]

        # create one gfunc conf file specific per organism
        gfunc_conf_file = '%s/conf/%s_gfunc.conf' % (options.root, common_name)
        gfunc_conf = open(gfunc_conf_file, 'w')
        gfunc_conf.write(TEMPLATE % {'root':options.root,
                                     'common_name':common_name,
                                     'embl_file':embl_file})
        gfunc_conf.close()
        print "GFUNC conf file %s has been generated." % gfunc_conf_file

        # update top conf file
        top_conf_file.write("GFUNC\t%s\n" % (gfunc_conf_file))

    top_conf_file.close()


### ---------------------------------------------------------------------------
if __name__ == '__main__':
    main()

