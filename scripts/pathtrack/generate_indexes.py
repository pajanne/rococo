'''
Created on Mar 26, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.

Read a Google spreadsheet with:
project||genus||species-subspecies||strain||sample||library||run||paired||insert_size||path_to_sff||path_to_assembly
and creates a sequence index and a sample info used as inputs to update_vrmeta.pl
and creates a assembly index used as input to import_assembly.py
if sample and library are undefined (='None'), sample = strain and library = run

see http://scratchy.internal.sanger.ac.uk/wiki/index.php/Creating_a_fake_sequence.index_for_making_a_VRTrack_meta_database

Usage:
> python ~ap12/genlibpy/genepy/pathtrack/generate_indexes.py -r /lustre/scratch103/pathogen/pathpipe/ -c metahit
'''

# Script: generate_indexes.py

# Topic: Name


# Topic: Usage
# > python ~ap12/genlibpy/genepy/pathtrack/generate_indexes.py -r /lustre/scratch103/pathogen/pathpipe/ -c metahit

# Topic: Description
# Read a Google spreadsheet with
# project||genus||species-subspecies||strain||sample||library||run||paired||insert_size||path_to_sff||path_to_assembly
# and creates a sequence index and a sample info used as inputs to update_vrmeta.pl
# and creates a assembly index used as input to import_assembly.py
# if sample and library are undefined (='None'), sample = strain and library = run

# Topic: Arguments

# Topic: See also
# http://scratchy.internal.sanger.ac.uk/wiki/index.php/Creating_a_fake_sequence.index_for_making_a_VRTrack_meta_database

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

from genepy import logsetup
from optparse import OptionParser
from genepy import util
import sys, os
from re import sub
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

    # check output path
    util.checkDir(options.root)
    out_metadata = "%s/.tmp/metadata/%s" % (options.root, options.category)
    util.checkDir(out_metadata)
    out_fastq = "%s/.tmp/fastq/%s" % (options.root, options.category)
    util.checkDir(out_fastq)

    # open output files
    sequence_index = open('%s/sequence.index' % out_metadata, 'w')
    samples_info = open('%s/samples.info' % out_metadata, 'w')
    assembly_index = open('%s/assembly.index' % out_metadata, 'w')

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
        if sample == 'None':
            sample = strain
        library = values[5]
        run = values[6]
        if library == 'None':
            library = run
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

        species_strain = sub('[^a-zA-Z0-9_]', '_', "%s_%s" % (species, strain)).replace('__', '_')
        strain_4_hierarchy = sub('[^a-zA-Z0-9_]', '_', strain).replace('__', '_')
        study = "%s_%s%s" % (values[0], genus[0], species_strain)

        instrument_platform = '454'
        empty = 'n/a'
        fastq_1_gz = "%s/%s_1.fastq.gz" % (out_fastq, run)
        fastq_2_gz = "%s/%s_2.fastq.gz" % (out_fastq, run)
        fastq_0_gz = "%s/%s.fastq.gz" % (out_fastq, run)

        # check that project name (study) is less than 40 char
        # mysql> desc project;
        # | name           | varchar(40)           | NO   | MUL |         |                | 
        # | hierarchy_name | varchar(40)           | NO   | MUL |         |                | 
        if len(study) > 40:
            log.warning("Project name %s has more than 40 char." % study)

        # checking files
        util.checkFile(sff_file)
        util.checkFile(trim_status)
        util.checkFile(contigs_file)
        util.checkFile(scaffolds_file)

        log.info("> checking fastq files: %s" % run)
        # get lane hierarchy path (run)
        run_path = util.runProcess("/software/bin/perl /nfs/users/nfs_a/ap12/genlibpy/genepy/pathtrack/get_lane_hierarchy_path.pl --lane=%s --db=%s" % (run, constants.DATABASE[options.category])).strip()

        # check if fastq files have been loaded in db
        do_generate_indexes  = False
        do_generate_assembly_indexes = False
        if run_path != "undefined":
            log.info("  loaded in db.")
            fastq_path = "%s/%s/seq-pipelines/%s" % (options.root, options.category, run_path)
            util.checkDir(fastq_path)
            # check if fastq files have been imported into the hierarchy
            if not (os.path.exists("%s/%s_1.fastq.gz" % (fastq_path, run)) and os.path.exists("%s/%s_2.fastq.gz" % (fastq_path, run)) and os.path.exists("%s/%s.fastq.gz" % (fastq_path, run))):
                log.info("  not imported into hierarchy.")
                # check if fastq files have been generated from sff into tmp dir
                if not (os.path.exists(fastq_1_gz) and os.path.exists(fastq_2_gz) and os.path.exists(fastq_0_gz)):
                    log.info("  not generated from sff files.")
                else:
                    log.info("  generate indexes.")
                    do_generate_indexes = True
            else:
                log.info("  already imported into hierarchy.")
                do_generate_assembly_indexes = True
        else:
            log.info("  not loaded in db.")
            # check if fastq files have been generated from sff into tmp dir
            if not (os.path.exists(fastq_1_gz) and os.path.exists(fastq_2_gz) and os.path.exists(fastq_0_gz)):
                log.info("  not generated from sff files.")
            else:
                log.info("  generate indexes.")
                do_generate_indexes = True
        
        # generate sequence and sample indexes
        if do_generate_indexes:
            # calculate md5
            md5_1 = util.runProcess("md5sum %s | cut -d ' ' -f 1" % fastq_1_gz).strip()
            md5_2 = util.runProcess("md5sum %s | cut -d ' ' -f 1" % fastq_2_gz).strip()
            md5_0 = util.runProcess("md5sum %s | cut -d ' ' -f 1" % fastq_0_gz).strip()

            # write to output files
            # sequence.index: fastq_file|md5|run_id|study_id|(study_name)|center_name|(submission_id)|(submission_date)|sample_id|sample_name|
            #                 (population)|(experiment_id)|instrument_platform|(instrument_model)|library_name|(run_name)|(run_block_name)|
            #                 insert_size|(library_layout)|paired_fastq|withdrawn|(withdrawn_date)|(comment)|read_count|base_count
            sequence_index.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                                 % (fastq_1_gz, md5_1, run, study, study, 'SC', empty, empty, sample, sample, strain, empty, instrument_platform,
                                    empty, library, empty, empty, insert_size, 'PAIRED', fastq_2_gz, '0', empty, empty, '0', '0'))
            sequence_index.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                                 % (fastq_2_gz, md5_2, run, study, study, 'SC', empty, empty, sample, sample, strain, empty, instrument_platform,
                                    empty, library, empty, empty, insert_size, 'PAIRED', fastq_1_gz, '0', empty, empty, '0', '0'))
            sequence_index.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                                 % (fastq_0_gz, md5_0, run, study, study, 'SC', empty, empty, sample, sample, strain, empty, instrument_platform,
                                    empty, library, empty, empty, insert_size, 'SINGLE', '', '0', empty, empty, '0', '0'))

            # samples.info: lookup_name|acc|individual_name|alias|population_name|species_name|taxon_id|sex
            samples_info.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (strain, strain, strain, strain, strain, organism_name, empty, empty))
        
        # generate assembly indexes
        if do_generate_indexes or do_generate_assembly_indexes:
            # assembly.index: genus||species_strain||assembly_id||contig||scaffold||fastq
            # /pyrodata01/assemblies/Ruminococcus/obeum/A2162/P_2009_07_12_22_16_23_runAssembly/454TrimStatus.txt
            assembly_id = "newbler_%s" % trim_status.split('/P_')[1][:10]
            assembly_index.write("%s||%s||%s||%s||%s||%s||%s\n" % (study, genus, species_strain, assembly_id, contigs_file, scaffolds_file, run))

    # close files
    sequence_index.close()
    samples_info.close()
    assembly_index.close()

### ---------------------------------------------------------------------------
if __name__ == '__main__':
    main()
