'''
Created on Mar 26, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.

Read a Google spreadsheet with:
project||genus||species-subspecies||strain||sample||library||run||paired||insert_size||path_to_sff||path_to_assembly
and convert sff file into fastq
if sample and library are undefined (='None'), sample = strain and library = run

Usage:
> python ~ap12/genlibpy/genepy/pathtrack/convert_sff2fastq.py -r /lustre/scratch103/pathogen/pathpipe/ -c metahit
'''
### ---------------------------------------------------------------------------
### Documentation header
# Script: convert_sff2fastq.py
# Read a Google spreadsheet and convert sff file into fastq

# Topic: Name


# Topic: Usage
# command line
# > python ~ap12/genlibpy/genepy/pathtrack/convert_sff2fastq.py -r /lustre/scratch103/pathogen/pathpipe/ -c metahit

# Topic: Description
# Read a Google spreadsheet and convert sff file into fastq
#
# Spreadsheet format
# > project||genus||species-subspecies||strain||sample||library||run||paired||insert_size||path_to_sff||path_to_assembly
#
# if sample and library are undefined (='None'), sample = strain and library = run

# Topic: Arguments

# Topic: See also

# Topic: Author
# Anne Pajon (ap12)

# Topic: Creation date
# Mar 26, 2010

# Topic: Copyright and License
# Copyright (C) 2009 Genome Research Limited. All Rights Reserved.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
### ---------------------------------------------------------------------------

from genepy import logsetup
from optparse import OptionParser
from genepy import util
import sys, os
import gdocs
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
    parser.add_option("-r", "--root", metavar="PATH", help="PATH to the root of the hierarchy", action="store", type="string", dest="root")
    parser.add_option("-c", "--category", metavar="CATEGORY", help="name of the category from %s" % constants.CATEGORY, action="store", choices=constants.CATEGORY, dest="category")
    
    (options, args) = parser.parse_args()

    if not (options.root and options.category):
        parser.print_help()
        sys.exit()
    
    # get the data from Goodgle spreadsheet
    lines = gdocs.getValues(doc='%s_454_Projects' % options.category.title())

    # check .tmp/fastq/ output path
    util.checkDir(options.root)
    outpath = "%s/.tmp/fastq/%s" % (options.root, options.category)
    util.checkDir(outpath)

    # process input data
    sample_count = 0
    for line in lines:
        line = line.strip()
        values = line.split('||')
        genus = values[1]
        species = values[2]
        strain = values[3]
        organism_name = "%s %s %s" % (genus, species, strain)
        sample_count = sample_count + 1
        sample = values[4]
        library = values[5]
        run = values[6]
        if values[7] == '1':
            paired = 'PAIRED'
        else:
            paired = 'SINGLE'
            log.error("Single read set for %s. Not implemented." % run)
            continue
        insert_size = values[8]
        sff_file = "/nfs/%s" % values[9]
        trim_status = "/nfs/%s/454TrimStatus.txt" % values[10]
        contigs_file = "/nfs/%s/454LargeContigs.fna" %values[10]
        scaffolds_file = "/nfs/%s/454Scaffolds.fna" %values[10]

        # check files
        util.checkFile(sff_file)
        util.checkFile(trim_status)
        util.checkFile(contigs_file)
        util.checkFile(scaffolds_file)

        # prepare commands to be executed
        dir_mh12scripts = "/nfs/users/nfs_m/mh12/git/python"
        ## convert sff into fastq
        outprefix = "%s/%s" % (outpath, run)
        cmd_sff2fastq = "%s/454TrimStatus2reads.py --pair_suffix=/1,/2 --sff %s %s %s" % (dir_mh12scripts, sff_file, trim_status, outprefix)
        fastq_pairs = "%s-pairs.fastq" % outprefix
        fastq_single = "%s-single.fastq" % outprefix
        ## split fastq pairs file
        fastq_1 = "%s_1.fastq" % outprefix
        fastq_2 = "%s_2.fastq" % outprefix
        cmd_splitfastq = "%s/fastn_unshuffle.py %s %s %s" % (dir_mh12scripts, fastq_pairs, fastq_1, fastq_2)
        ## rename fastq single file
        fastq_0 = "%s.fastq" % outprefix
        cmd_rename = "mv %s %s" % (fastq_single, fastq_0)
        ## tidy-up
        cmd_remove = "rm %s-info.txt; rm %s-pairs.fastq" % (outprefix, outprefix)
        ## gzip fastq files
        cmd_gzip = "gzip %s; gzip %s; gzip %s" % (fastq_1, fastq_2, fastq_0)
        ## all commands
        cmd = "%s; %s; %s; %s; %s" % (cmd_sff2fastq, cmd_splitfastq, cmd_rename, cmd_remove, cmd_gzip)

        log.info("> checking fastq files: %s" % run)
        # get lane hierarchy path (run)
        run_path = util.runProcess("/software/bin/perl /nfs/users/nfs_a/ap12/genlibpy/genepy/pathtrack/get_lane_hierarchy_path.pl --lane=%s --db=%s" % (run, constants.DATABASE[options.category])).strip()

        # check if fastq files have been loaded in db
        do_convert = False
        if run_path != "undefined":
            log.info("  loaded in db.")
            fastq_path = "%s/%s/seq-pipelines/%s" % (options.root, options.category, run_path)
            util.checkDir(fastq_path)
            # check if fastq files have been imported into the hierarchy
            if not (os.path.exists("%s/%s_1.fastq.gz" % (fastq_path, run)) and os.path.exists("%s/%s_2.fastq.gz" % (fastq_path, run)) and os.path.exists("%s/%s.fastq.gz" % (fastq_path, run))):
                log.info("  not imported into hierarchy.")
                # check if fastq files have been generated from sff into tmp dir
                if not (os.path.exists("%s.gz" % fastq_1) and os.path.exists("%s.gz" % fastq_2) and os.path.exists("%s.gz" % fastq_0)):
                    do_convert = True
                else:
                    log.info("  already generated from sff files.")
            else:
                log.info("  already imported into hierarchy.")
        else:
            log.info("  not loaded in db.")
            # check if fastq files have been generated from sff into tmp dir
            if not (os.path.exists("%s.gz" % fastq_1) and os.path.exists("%s.gz" % fastq_2) and os.path.exists("%s.gz" % fastq_0)):
                do_convert = True
            else:
                log.info("  already generated from sff files.")
        
        if constants.IS_LSF and do_convert:
            util.submitJob(jobname='sff2fastq_%s' % run, cmd=cmd, outdir=outpath)
                

### ---------------------------------------------------------------------------
if __name__ == '__main__':
    main()
