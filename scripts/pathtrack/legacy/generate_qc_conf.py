'''
Created on Jun 03, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.


Takes a tab delimited file of format:
genus||species_strain||assembly_id||contig||scaffold||run
and generate QC configuration files

Usage:
> python ~ap12/genlibpy/genepy/pathtrack/generate_qc_conf.py -a /lustre/scratch103/pathogen/pathpipe/tmp/metadata/metahit/assembly.index -r /lustre/scratch103/pathogen/pathpipe/ -c metahit

'''

from genepy import logsetup
from optparse import OptionParser
from genepy import util
import sys, os
import constants

### ---------------------------------------------------------------------------
### Logging configuration
### ---------------------------------------------------------------------------
import logging
log = logging.getLogger('genepy.pathtrack')

### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-a", "--assembly", metavar="FILE", help="FILE containing the list of all contigs and scaffolds to import", action="store", type="string", dest="assembly")
    parser.add_option("-r", "--root", metavar="PATH", help="PATH to the root of the hierarchy", action="store", type="string", dest="root")
    parser.add_option("-c", "--category", metavar="CATEGORY", help="name of the category from %s" % constants.CATEGORY, action="store", choices=constants.CATEGORY, dest="category")
    
    (options, args) = parser.parse_args()

    if not (options.assembly and options.root and options.category):
        parser.print_help()
        sys.exit()

    # check root path
    if not os.path.exists(options.root):
        log.error("%s path do not exist" % options.root)
        log.error("Create root path first, then run pipeline before importing assembly files.")
        sys.exit()
    
    # check log directory exists
    out_log = "%s/log/%s" % (options.root, options.category)
    util.checkDir(out_log)

    # open qc_pipeline.conf
    pipeline_qc = open('%s/conf/%s/qc_pipeline.conf' % (options.root, options.category), 'w')

    # check input assembly file and read it - one line per run (lane)
    util.checkFile(options.assembly)
    assembly_lines = open(options.assembly, "r").readlines()
    # compare project name - could have more than one run per project
    previous_project = ""
    for line in assembly_lines:
        if line[0] == '!':
            continue
        if not line.count('||') == 6:
            continue
        line = line.strip()
        values = line.split('||')
        project = values[0]
        genus = values[1]
        species = values[2]
        assembly_id = values[3]
        contig_file = values[4]
        scaffold_file = values[5]
        run = values[6]

        # check if new project
        if project != previous_project:
            # check if files are in place in the hierarchy
            species_path = "%s/%s/seq-pipelines/%s/%s" % (options.root, options.category, genus, species)
            assembly_path = "%s/ASSEMBLY" % species_path
            assembly_id_path = "%s/%s" % (assembly_path, assembly_id)
            scaffold_file_hierarchy = "%s/Scaffolds.fna" % assembly_id_path
            util.checkFile(scaffold_file_hierarchy)
            util.checkFile("%s.fai" % scaffold_file_hierarchy)
            util.checkFile("%s.refstats" % scaffold_file_hierarchy)
            util.checkFile("%s.bwt" % scaffold_file_hierarchy)

            # create one qc conf file specific per project
            qc_conf_filename = '%s/conf/%s/%s_qc.conf' % (options.root, options.category, project)
            qc_conf = open(qc_conf_filename, 'w')
            qc_conf.write(constants.QC_CONF_TEMPLATE % {'root':options.root,
                                                        'category':options.category,
                                                        'db':constants.DATABASE[options.category],
                                                        'db_host':os.getenv('VRTRACK_HOST'),
                                                        'db_port':os.getenv('VRTRACK_PORT'),
                                                        'db_rw_user':os.getenv('VRTRACK_RW_USER'),
                                                        'db_password':os.getenv('VRTRACK_PASSWORD'),
                                                        'project':project,
                                                        'ref':scaffold_file_hierarchy})
            qc_conf.close()

            log.info("QC conf file %s has been generated." % qc_conf_filename)

            # update qc_pipeline.conf
            pipeline_qc.write("__VRTrack_QC__\t%s\n" % (qc_conf_filename))

            # update previous project name
            previous_project = project

    pipeline_qc.close()
        


### ---------------------------------------------------------------------------
if __name__ == '__main__':
    main()

