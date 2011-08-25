'''
Created on Mar 26, 2010
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.

Takes a delimited file of format:
!project||genus||species-subspecies||strain||run||paired||insert_size||sff||assembly_path
and creates a sequence index and a sample info used as inputs to update_vrmeta.pl
and creates a assembly index used as input to import_assembly.py
where sample = strain and library = run

see http://scratchy.internal.sanger.ac.uk/wiki/index.php/Creating_a_fake_sequence.index_for_making_a_VRTrack_meta_database

---------------------------------------------------------------------------
- POPULATING THE DATABASE
---------------------------------------------------------------------------
Documentation on wiki @ http://scratchy.internal.sanger.ac.uk/wiki/index.php/Updating_the_databases#454_projects

- to generate the index files
> python ~ap12/genlibpy/genepy/pathtrack/454projects2dataIndexes.py -l 454PROJECTS.list -o /lustre/scratch103/pathogen/pathpipe/tmp
> python ~ap12/genlibpy/genepy/pathtrack/454projects2dataIndexes.py -l 454PROJECTS.list -o /lustre/scratch103/pathogen/pathpipe/tmp --fastq
> python ~ap12/genlibpy/genepy/pathtrack/454projects2dataIndexes.py -l 454PROJECTS.list -o /lustre/scratch103/pathogen/pathpipe/tmp --md5

- to run the update
> update_vrmeta.pl --index sequence.index --samples samples.info --database pathogen_metahit_track

---------------------------------------------------------------------------
- PIPELINE TO CREATE DATA HIERARCHY
---------------------------------------------------------------------------
- conf files
> xemacs /lustre/scratch103/pathogen/pathpipe/conf/toplevel_454pipeline.conf
> xemacs /lustre/scratch103/pathogen/pathpipe/conf/metahit/import_metahit.conf

- to run the pipeline
> run-pipeline -c /lustre/scratch103/pathogen/pathpipe/conf/toplevel_454pipeline.conf -v -v -o

'''

from genepy import logsetup
from optparse import OptionParser
from genepy import util
import sys, os
import hashlib
from re import sub


IS_LSF = util.isLsf()

### ---------------------------------------------------------------------------
### Logging configuration
### ---------------------------------------------------------------------------
import logging
log = logging.getLogger('genepy.pathtrack')

### ---------------------------------------------------------------------------
def main():
    usage = "usage: %prog [Options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", metavar="FILE", help="FILE containing the list of all 454 runs and its associated information", action="store", type="string", dest="list")
    parser.add_option("-o", "--outpath", metavar="PATH", help="PATH where to generate indexes and temporary fastq files", action="store", type="string", dest="outpath")
    parser.add_option("--fastq", help="Do generate fastq files.", action="store_true", dest="fastq")
    parser.add_option("--md5", help="Do run md5sum on generated fastq files.", action="store_true", dest="md5")
    
    (options, args) = parser.parse_args()

    if not options.list and not options.outpath:
        parser.print_help()
        sys.exit()
    
    # input file
    input_file = options.list
    util.checkFile(input_file)
    input_lines = open(input_file, "r").readlines()

    # output path
    output_path = options.outpath
    util.checkDir(output_path)
    out_metadata = "%s/metadata" % output_path
    util.checkDir(out_metadata)
    out_fastq = "%s/fastq" % output_path
    util.checkDir(out_fastq)

    # checking file format first before processing it
    lines = []
    for line in input_lines:
        if line[0] == '!':
            continue
        elif not line.count('||') == 8:
            log.error("line is not well formated. Please check your input file.")
            log.error(line.count('||'))
            log.error(line)
            sys.exit()
        else:
            lines.append(line)
    log.debug(lines)

    # opening output files
    sequence_index_filename = '%s/sequence.index' % out_metadata
    sequence_index = open(sequence_index_filename, 'w')
    samples_info = open('%s/samples.info' % out_metadata, 'w')
    assembly_index = open('%s/assembly.index' % out_metadata, 'w')

    # processing input file
    sample_count = 0
    for line in lines:
        line = line.strip()
        values = line.split('||')
        log.info(line)
        genus = values[1]
        species = values[2]
        strain = values[3]
        organism_name = "%s %s %s" % (genus, species, strain)
        sample_count = sample_count + 1
        run = values[4]
        if values[5] == '1':
            paired = 'PAIRED'
        else:
            paired = 'SINGLE'
        insert_size = values[6]
        sff_file = "/nfs/%s" % values[7]
        trim_status = "/nfs/%s/454TrimStatus.txt" % values[8]
        contigs_file = "/nfs/%s/454LargeContigs.fna" %values[8]
        scaffolds_file = "/nfs/%s/454Scaffolds.fna" %values[8]

        species_strain = sub('[^a-zA-Z0-9_]', '_', "%s_%s" % (species, strain)).replace('__', '_')
        strain_4_hierarchy = sub('[^a-zA-Z0-9_]', '_', strain).replace('__', '_')
        study = "%s_%s%s" % (values[0], genus[0], species_strain)

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

        # convert sff into fastq
        outprefix = "%s/%s" % (out_fastq, run)
        cmd_sff2fastq = "/nfs/users/nfs_m/mh12/svn-repository/pathogen/user/mh12/python/454TrimStatus2reads.py --pair_suffix=/1,/2 --sff %s %s %s" % (sff_file, trim_status, outprefix)
        fastq_pairs = "%s-pairs.fastq" % outprefix
        fastq_single = "%s-single.fastq" % outprefix

        # split fastq pairs file
        fastq_1 = "%s_1.fastq" % outprefix
        fastq_2 = "%s_2.fastq" % outprefix
        cmd_splitfastq = "/nfs/users/nfs_m/mh12/svn-repository/pathogen/user/mh12/python/fastn_unshuffle.py %s %s %s" % (fastq_pairs, fastq_1, fastq_2)

        # rename fastq single file
        fastq_0 = "%s.fastq" % outprefix
        cmd_rename = "mv %s %s" % (fastq_single, fastq_0)

        # tidy-up
        cmd_remove = "rm %s-info.txt; rm %s-pairs.fastq" % (outprefix, outprefix)

        # gzip fastq files
        cmd_gzip = "gzip %s; gzip %s; gzip %s" % (fastq_1, fastq_2, fastq_0)

        # all commands
        cmd = "%s; %s; %s; %s; %s" % (cmd_sff2fastq, cmd_splitfastq, cmd_rename, cmd_remove, cmd_gzip)


        if IS_LSF:
            if not (os.path.exists("%s.gz" % fastq_1) and os.path.exists("%s.gz" % fastq_2) and os.path.exists("%s.gz" % fastq_0)):
                if options.fastq:
                    util.submitJob(jobname='sff2fastq_%s' % run, cmd=cmd, outdir=out_metadata)
                else:
                    log.info("fastq files do not exist, use '--fastq' to generate them.")
            else:
                log.info("fastq files already exist.")
        else:
            log.info("Need to be run on LSF.")

        instrument_platform = '454'
        empty = 'n/a'
        fastq_1_gz = "%s.gz" % fastq_1
        fastq_2_gz = "%s.gz" % fastq_2
        fastq_0_gz = "%s.gz" % fastq_0

        # write to output files
        # sequence.index: fastq_file|md5|run_id|study_id|(study_name)|center_name|(submission_id)|(submission_date)|sample_id|sample_name|
        #                 (population)|(experiment_id)|instrument_platform|(instrument_model)|library_name|(run_name)|(run_block_name)|
        #                 insert_size|(library_layout)|paired_fastq|withdrawn|(withdrawn_date)|(comment)|read_count|base_count
        sequence_index.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                             % (fastq_1_gz, 'md5', run, study, study, 'SC', empty, empty, strain, strain, strain, empty, instrument_platform,
                                empty, run, empty, empty, insert_size, 'PAIRED', fastq_2_gz, '0', empty, empty, '0', '0'))
        sequence_index.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                             % (fastq_2_gz, 'md5', run, study, study, 'SC', empty, empty, strain, strain, strain, empty, instrument_platform,
                                empty, run, empty, empty, insert_size, 'PAIRED', fastq_1_gz, '0', empty, empty, '0', '0'))
        sequence_index.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                             % (fastq_0_gz, 'md5', run, study, study, 'SC', empty, empty, strain, strain, strain, empty, instrument_platform,
                                empty, run, empty, empty, insert_size, 'SINGLE', '', '0', empty, empty, '0', '0'))

        # samples.info: lookup_name|acc|individual_name|alias|population_name|species_name|taxon_id|sex
        samples_info.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (strain, strain, strain, strain, strain, organism_name, empty, empty))

        # assembly.index: genus||species_strain||assembly_id||contig||scaffold||fastq
        # /pyrodata01/assemblies/Ruminococcus/obeum/A2162/P_2009_07_12_22_16_23_runAssembly/454TrimStatus.txt
        assembly_id = "newbler_%s" % trim_status.split('/P_')[1][:10]
        assembly_index.write("%s||%s||%s||%s||%s||%s||%s\n" % (study, genus, species_strain, assembly_id, contigs_file, scaffolds_file, run))

    # close files
    sequence_index.close()
    samples_info.close()
    assembly_index.close()

    if not options.fastq:
        log.info("Use '--fastq' for generating fastq files")

    if options.md5:
        # calculate md5 and modify sequence.index
        util.checkFile(sequence_index_filename)
        seq_lines = open(sequence_index_filename, "r").readlines()
        sequence_index = open(sequence_index_filename, 'w')
        for line in seq_lines:
            values = line.split('\t')
            fastq = values[0]
            run = values[2]
            if os.path.exists(fastq):
                md5 = util.runProcess("md5sum %s | cut -d ' ' -f 1" % fastq).strip()
                line = line.replace('md5', md5)
                sequence_index.write(line)
            else:
                log.info("fastq file %s does not exist, use '--fastq' for generating it." % fastq)

        # close file
        sequence_index.close()
    else:
        log.info("When all submitted jobs end, use '--md5' for updating sequence.index with md5sum.")

        
    
### ---------------------------------------------------------------------------
if __name__ == '__main__':
    main()
