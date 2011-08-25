'''
Created on Apr 21, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.


Takes a tab delimited file of format:
genus||species_strain||assembly_id||contig||scaffold||run
and import contig and scaffold files into the hierarchy
and generate required stats files needed for the QC pipeline

Usage:
> python ~ap12/genlibpy/genepy/pathtrack/import_assembly.py -a /lustre/scratch103/pathogen/pathpipe/tmp/metadata/metahit/assembly.index -r /lustre/scratch103/pathogen/pathpipe/ -c metahit

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
def has_same_md5(ori_file, copied_file):
    cmd = "md5sum %s"
    ori_md5 = util.runProcess(cmd % ori_file).split()[0]
    copied_md5 = util.runProcess(cmd % copied_file).split()[0]
    if ori_md5 == copied_md5:
        return True
    else:
        return False

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
    
    # check input assembly file and read it - one line per run (lane)
    util.checkFile(options.assembly)
    assembly_lines = open(options.assembly, "r").readlines()
    # compare project name - could have more than one run per project
    previous_project = ""
    is_new_project = True
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

        # get lane hierarchy path (run)
        run_path = util.runProcess("/software/bin/perl /nfs/users/nfs_a/ap12/genlibpy/genepy/pathtrack/get_lane_hierarchy_path.pl --lane=%s --db=%s" % (run, constants.DATABASE[options.category])).strip()
        fastq_file = "%s/%s/seq-pipelines/%s" % (options.root, options.category, run_path)

        # check if new project
        if project == previous_project:
            is_new_project = False
        else:
            previous_project = project
            is_new_project = True

        # check species path
        species_path = "%s/%s/seq-pipelines/%s/%s" % (options.root, options.category, genus, species)
        if not os.path.exists(species_path):
            log.error("%s path do not exist" % species_path)
            log.error("Run fastq import pipeline before importing assembly files.")
        else:
            # create assembly path
            assembly_path = "%s/ASSEMBLY" % species_path
            if not os.path.exists(assembly_path):
                os.makedirs(assembly_path)
                log.info("%s created" % assembly_path)
            else:
                log.info("%s path already exists" % assembly_path)

            # create assembly_id path (newbler_2009_06_29)
            assembly_id_path = "%s/%s" % (assembly_path, assembly_id)
            if not os.path.exists(assembly_id_path):
                os.makedirs(assembly_id_path)
                log.info("%s created" % assembly_id_path)
            else:
                log.info("%s path already exists" % assembly_id_path)

            # copy contigs file
            contig_file_hierarchy = "%s/LargeContigs.fna" % assembly_id_path
            util.checkFile(contig_file)
            cmd_cp = "cp %s %s" % (contig_file, contig_file_hierarchy)
            if not os.path.exists(contig_file_hierarchy):
                util.runProcess(cmd_cp)
                if not has_same_md5(contig_file, contig_file_hierarchy):
                    log.error("Copied file %s is not the same as original file %s" (contig_file, contig_file_hierarchy))
            else:
                log.info("%s file already exists" % contig_file_hierarchy)
                if not has_same_md5(contig_file, contig_file_hierarchy):
                    log.error("Copied file %s is not the same as original file %s" (contig_file, contig_file_hierarchy))

            # copy scaffolds file
            scaffold_file_hierarchy = "%s/Scaffolds.fna" % assembly_id_path
            util.checkFile(scaffold_file)
            cmd_cp = "cp %s %s" % (scaffold_file, scaffold_file_hierarchy)
            if not os.path.exists(scaffold_file_hierarchy):
                util.runProcess(cmd_cp)
                if not has_same_md5(scaffold_file, scaffold_file_hierarchy):
                    log.error("Copied file %s is not the same as original file %s" (scaffold_file, scaffold_file_hierarchy))
            else:
                log.info("%s file already exists" % scaffold_file_hierarchy)
                if not has_same_md5(scaffold_file, scaffold_file_hierarchy):
                    log.error("Copied file %s is not the same as original file %s" (scaffold_file, scaffold_file_hierarchy))

            # create fastqs path
            fastqs_path = "%s/fastqs" % assembly_id_path
            if not os.path.exists(fastqs_path):
                os.makedirs(fastqs_path)
                log.info("%s created" % fastqs_path)
            else:
                log.info("%s path already exists" % fastqs_path)
            
            # create simlinks to fastqs
            util.checkDir(fastq_file)
            fastq_name = run
            symlink = "%s/%s" % (fastqs_path, fastq_name)
            if not os.path.exists(symlink):
                os.symlink(fastq_file, symlink)
                log.info("%s symlink created" % symlink)
            else:
                log.info("%s symlink already exists" % symlink)
                
            # run samtools faidx, refstats, bwa to generate extra files required for the QC pipeline
            cmd = ""
            if not os.path.exists("%s.fai" % scaffold_file_hierarchy):
                cmd = "samtools faidx %s; " % scaffold_file_hierarchy
            else:
                log.info("%s.fai already exists" % scaffold_file_hierarchy)
            
            if not os.path.exists("%s.refstats" % scaffold_file_hierarchy):
                cmd = cmd + "ref-stats -r %s > %s.refstats; " % (scaffold_file_hierarchy, scaffold_file_hierarchy)
            else:
                log.info("%s.refstats already exists" % scaffold_file_hierarchy)
            
            if not os.path.exists("%s.bwt" % scaffold_file_hierarchy):
                cmd = cmd + "bwa index %s; " % scaffold_file_hierarchy
            else:
                log.info("%s.bwt already exists" % scaffold_file_hierarchy)

            # run stats.py on contigs and scaffolds
            cmd_stats = "python /nfs/users/nfs_a/ap12/genlibpy/genepy/pathtrack/stats.py -f %s; "
            if not os.path.exists("%s.stats" % contig_file_hierarchy):
                cmd = cmd + cmd_stats % contig_file_hierarchy
            else:
                log.info("%s.stats already exists" % contig_file_hierarchy)
            if not os.path.exists("%s.stats" % scaffold_file_hierarchy):
                cmd = cmd + cmd_stats % scaffold_file_hierarchy
            else:
                log.info("%s.stats already exists" % scaffold_file_hierarchy)

            # submit all jobs
            if is_new_project and not cmd == "":
                util.submitJob(jobname='stats_%s_%s' % (project, assembly_id), cmd=cmd, outdir=assembly_id_path)

### ---------------------------------------------------------------------------
if __name__ == '__main__':
    main()

