#!/usr/local/bin/perl -w

# Copyright (C) 2005 Genome Research Limited. All Rights Reserved.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

=head1 Name

bacterial_pre-annotator.pl: does automated annotation of genes using HAMAP
                            and reciprocal fasta against in-house bacterial
							genomes

=head1 Synopsis

usage: bacterial_pre-annotator.pl [-q qualifier] [-r 'regex'] -e embl_file -o output_name
   or: bacterial_pre-annotator.pl [-q qualifier] [-r 'regex'] -t tab_file -d dna_file -o output_name
   or: bacterial_pre-annotator.pl -h (for help; explanation from perldoc $0)
  
  -q   the qualifier that contains the systematic id, eg. note or gene.
       default: locus_tag
  -r   a regular expression that can be used to detect the systematic id in
       the given qualifier. must be enclosed in single quotes. must have one
	   pair of brackets to capture the id. default: '(\\S+)'
	   
If you have separate files for each chromosome and can not join them, provide
a list of relevant file names by space separating them within double quotes.
When using -t and -d, there must be an equal number of files given to each, and
the first tab file listed in -t must have its corresponding dna file as the
first file listed in -d and so on.

=head1 Author

Sendu Bala (email: bix@sendu.me.uk)

=cut

################################################################################
# load modules
use strict;
use warnings;
use Getopt::Std;
use File::Basename;
use File::PSU::HashedDirs;
use LSF::PSU::JobArray;
use Bio::PSU::SeqFactory;
use Bio::PSU::SearchFactory;
use Bio::PSU::FastFasta;

# we use Data::Dumper instead of DB_File et al so we can store complex data structure,
# and instead of storable so it is compatible with different perl versions and
# architectures (hopefully future-proofed!).
use Data::Dumper;

# for log
use Log::Log4perl qw(get_logger);

# for config
use AppConfig qw(:argcount);

# catch ctrl-c so we can delete the lock file before exiting
$SIG{INT} = \&catch_ctrlc;

################################################################################
# show a usage message and then exit
sub usage {
    die <<EOF;

bacterial_pre-annotator.pl: does automated annotation of genes using HAMAP
                            and reciprocal fasta against in-house bacterial
			    genomes

usage: $0 [-q qualifier] [-r 'regex'] -e embl_file -o output_name
   or: $0 [-q qualifier] [-r 'regex'] -t tab_file -d dna_file -o output_name
   or: $0 -h (for help; explanation from perldoc $0)
  
  -q   the qualifier that contains the systematic id, eg. note or gene.
       default: locus_tag
  -r   a regular expression that can be used to detect the systematic id in
       the given qualifier. must be enclosed in single quotes. must have one
	   pair of brackets to capture the id. default: '(\\S+)'
	   
If you have separate files for each chromosome and can not join them, provide
a list of relevant file names by space separating them within double quotes.
When using -t and -d, there must be an equal number of files given to each, and
the first tab file listed in -t must have its corresponding dna file as the
first file listed in -d and so on.
  
EOF
}

################################################################################
# Configuration
my $config = AppConfig->new({
    CASE   => 1,
    GLOBAL => { 
        ARGCOUNT => ARGCOUNT_ONE,
    },
    PEDANTIC => 1,
});
# define all the variables we will use, with defaults were necessary
$config->define("genomes_dir", {
	        DEFAULT => '/nfs/pathdata/BacterialGenomes/',});
$config->define("tmp_storage_dir", {
	        DEFAULT => '/lustre/pathogen/scratch/programs/bacterial_pre_annotator/',});
$config->define("hamap_profile", {
	        DEFAULT => '/data/blastdb/expasy/hamap.prf',});
$config->define("hamap_annotation", {
	        DEFAULT => '/data/blastdb/expasy/alignment_id.dat',});
$config->define("log4perl_level", {
                DEFAULT => 'INFO',
		VALIDATE => 'DEBUG|INFO|WARN|ERROR|FATAL',});
# process file if it's available
my $config_file = "bacterial_pre-annotator.conf";
if (-f $config_file) {
    $config->file($config_file);
    print "Config file ", $config_file, " read.\n";
} else {
    print "Config file ", $config_file, " not found.\n";
    print "Default values used\n";
}

################################################################################
# Log
my $log_conf = q/ 
    log4perl.category = / . $config->log4perl_level() . q/, Logfile, Screen 
     
    log4perl.appender.Logfile = Log::Log4perl::Appender::File 
    log4perl.appender.Logfile.filename = sub { return get_log_fn(); }
    log4perl.appender.Logfile.mode   = write 
    log4perl.appender.Logfile.layout = Log::Log4perl::Layout::SimpleLayout
     
    log4perl.appender.Screen        = Log::Log4perl::Appender::Screen 
    log4perl.appender.Screen.layout = Log::Log4perl::Layout::SimpleLayout
/;

Log::Log4perl::init( \$log_conf );
my $logger = Log::Log4perl::get_logger();

$logger->info("Starting $0");

$logger->debug("genomes_dir = ", $config->genomes_dir());
$logger->debug("tmp_storage_dir = ", $config->tmp_storage_dir());
$logger->debug("hamap_profile = ", $config->hamap_profile());
$logger->debug("hamap_annotation = ", $config->hamap_annotation());
$logger->debug("log4perl_level = ", $config->log4perl_level());


################################################################################
# set up defaults and base variables
my $genomes_dir = $config->genomes_dir();
my $dbs_dir = $genomes_dir.'scriptdbs/';
my $lock_dir = $dbs_dir.'lock';
my $updates_file = $dbs_dir.'updates.db';
-e $updates_file or die "error: no genome database information was found; use bacterial_genome_database.pl first\n";
my $qualifier = 'locus_tag';
my $regex = '(\S+)';
my ($sec,$min,$hour,$month_day,$month,$year) = localtime(time);
$month++;
$year += 1900;
$month = sprintf("%02d", $month);
$month_day = sprintf("%02d", $month_day);
my $date = $year.$month.$month_day;
my $unique = $$;

# fasta variables from original reciprocal_fasta.pl script
my $fasta_args = '-z 1 -Q -H -S -m 10';
my $percent_id_cutoff = 30;
my $percent_hit_length_cutoff = 80;

# different cutoffs for the top fasta
my $top_evalue_cutoff = 1e-10;

# HAMAP variables
my $hamap_profile = $config->hamap_profile();
my $hamap_annotation = $config->hamap_annotation();
my $hamap_bin = '/nfs/users/nfs_a/ap12/bin/pfscan';
my $hamap_args = '-klf';
my $hamap_cutoff_level = 1;

# provide usage message in default cases
if (@ARGV > 0 && $ARGV[0] eq "-h") {
    usage;
}

# option handling, die with usage message on errors
my $input_options = "@ARGV";
my %options = ();
getopts ("q:r:e:t:d:o:b:", \%options);

my ($output_name, @embl_files, @tab_files, @dna_files);

if ((defined $options{t}) + (defined $options{d}) == 1) {
    warn "error: -t and -d must both be specified\n\n";
    usage;
}

if (defined $options{e} && defined $options{t}) {
    warn "error: -e and -t are mutually exclusive\n\n";
    usage;
}

if ((defined $options{o}) + (defined $options{e}) != 2 && (defined $options{o}) + (defined $options{t}) != 2) {
    warn "error: -o and -e or (-t and -d) are required\n\n";
    usage;
}

defined $options{o} and $output_name = $options{o};
my $error_file = $output_name.'.errors';
defined $options{q} and $qualifier = $options{q};
defined $options{r} and $regex = $options{r};
defined $options{e} and @embl_files = split(' ', $options{e});
defined $options{t} and @tab_files = split(' ', $options{t});
defined $options{d} and @dna_files = split(' ', $options{d});

# check the files exist
foreach my $file (@embl_files, @tab_files, @dna_files) {
    -s $file or die "error: input file '$file' is empty or can't be found\n";
}

# this is an internal-only use option that we use for the bsubed bits
defined $options{b} and $unique = $options{b};
my $tmp_storage_dir = $config->tmp_storage_dir();
my $unique_dir = $tmp_storage_dir.'bpa'.$unique;
my $bsub_fofn = $unique_dir.'/fofn';
my $single_file_base_dir = $unique_dir.'/sequences';
my $fasta_result_base_dir = $unique_dir.'/fasta_results';
my $genome_base_name = $unique_dir.'/genome';
my $error_base_name = $unique_dir.'/bsub_err';
my $output_base_name = $unique_dir.'/bsub_out';
my $max_sim_jobs = 16;

my $jobarray_dir = $unique_dir; # for hamap JobArray

(-d $unique_dir or mkdir($unique_dir)) or die "error: could not create storage dir '$unique_dir'\n";

$logger->info("Variables set up");


################################################################################
# the main code including reciprocal fastas is only carried out once (not each
# time this script is called in the bsub array for the hamap)
my @single_files;
my @genomes;
my %fastfastas;
my %reciprocal_hits;
my %best_hits;
my $clear_locks = 0;
unless (defined $options{b}) {
    # if someone else is using one of the bacterial_*.pl scripts, wait until they finish,
    # else make the lock file to stop others using this
    unless (check_locks()) {
	my $warn_count = 0;
	while (! mkdir($lock_dir, 0777)) {
	    $logger->error("ERROR:::: $!");
	    sleep(2);
	    $warn_count++;
	    if ($warn_count == 6) {
	    $logger->info("waiting for another user to finish their update...\n(you also need write access to /nfs/pathdata/BacterialGenomes)\n(as a last resort delete /nfs/pathdata/BacterialGenomes/scriptdbs/lock)");
	}
    }
}
$clear_locks = 1;
	
# now advertise to other bpas that we're still running and need the main lock
mkdir($lock_dir.$unique, 0777);

$logger->info("Lock checked and set up");

@genomes = get_inhouse_genomes($updates_file);
$logger->info("Array of in-house genomes created");


################################################################################
# get_inhouse_genomes
#
# get an array of ready-to-use in-house database names from the
# bacterial_genome_database.pl updates hash
#
################################################################################
sub get_inhouse_genomes {
    my $local_file = $_[1];
    open(DB, $local_file) or die "error: couldn't read from updates database '$local_file'\n";
    my $dumped;
    {
        local $/ = undef;
	$dumped = <DB>;
	close(DB);
    }
    my $VAR1;
    eval $dumped;
    die "database file '$updates_file' could not be understood\n" if $@;
		
    foreach my $key (keys %{$VAR1}) {
        push(@genomes, $key) unless $key eq 'bgd_warnings';
    }
    return @genomes;
}

	
################################################################################
# create symlinks to the genome files with simple incrementing-number based
# filenames to make the bsub array later on simple to set up.
# also make FastaFasta objects for each genome for fast extraction of specific
# proteins from them later.
my $genome_count = 0;
foreach my $genome (@genomes) {
    $genome_count++;
    $logger->info($genomes_dir . $genome);
    $logger->debug("ln -s $genomes_dir$genome.fa $genome_base_name$genome_count");
    system("ln -s $genomes_dir$genome.fa $genome_base_name$genome_count") == 0 or die "error: ln -s failed: $?\n";
    $fastfastas{$genome} = Bio::PSU::FastFasta->new( -infile => "$genomes_dir$genome.fa" );
}
$logger->info("Temporary symlinks to these genome files created");
$logger->info("FastFasta for each genome done");


################################################################################
# create a database of the proteins from the new given genome by calling bacterial_genome_database.pl
my $new_genome;
my $bgd_options = "-q $qualifier -r '($regex)' ";
if (@embl_files > 0) {
    $bgd_options .= "-e \"@embl_files\"";
    $new_genome .= basename($embl_files[0]);
} else {
    $bgd_options .= "-t \"@tab_files\" -d \"@dna_files\"";
    $new_genome .= basename($tab_files[0]);
}
$new_genome =~ s/\.(?:embl|tab|art)/.bpa_prot_db/;
$new_genome = $unique_dir.'/'.$new_genome;
unless (-e "$new_genome.fa") {
    $logger->debug("bacterial_genome_database.pl -o $new_genome $bgd_options");
    system("bacterial_genome_database.pl -o $new_genome $bgd_options") == 0 or die "error: bacterial_genome_database.pl failed: $?\n";
} else {
    warn "warning: '$new_genome.fa' already exists, not re-creating protein database from input embl/tab file(s)\n";
}
$new_genome .= '.fa';
$logger->debug("new genome: ", $new_genome);	

# we will store each single-file filename in a file of filenames for later use in hamap bsub array
open(FOFN, ">$bsub_fofn") or die "error: could not write to file '$bsub_fofn'\n";

$logger->info("Database of proteins from the new genome created");


################################################################################
# extract each sequence in the new genome into separate files
{
    #local $/ = "\n>";
    local $/ = ">";
    open(ALL, $new_genome) or die "error: could not open protein database '$new_genome'\n";
    while (<ALL>) {
	# print out a single sequence to file
	#s/^>|>$//;
	chomp;
	next unless length;
	my ($single_seq_file) = $_ =~ /^(\S+)/;

	# store each file within a hashed directory structure so we don't have too many files in a directory
	my $single_seq_file_dir = create_hashed_dirs($single_file_base_dir, $single_seq_file);
	$single_seq_file = $single_seq_file_dir.'/'.$single_seq_file;
	push(@single_files, $single_seq_file);
		
	open(SINGLE, ">$single_seq_file") or die "error: could not write to '$single_seq_file'\n";
	print SINGLE '>', $_;
	close(SINGLE);

	$logger->debug("$single_seq_file");
	$logger->debug("$_");
	
	# store its location in the fofn file
	print FOFN $single_seq_file, "\n";
    }
}
close(FOFN);
$logger->info("Extracted each sequences into separate files");


################################################################################
# run fasta of all the above sequences against each in-house database, simultaneously
$logger->debug("bsub -q basement -J bpa${unique}l1\"[1-$genome_count]\%$max_sim_jobs\" -o $output_base_name\%I -e $error_base_name\%I 'many_fasta.pl $fasta_args $bsub_fofn $genome_base_name\${LSB_JOBINDEX} $fasta_result_base_dir\${LSB_JOBINDEX}'");
system("bsub -q basement -J bpa${unique}l1\"[1-$genome_count]\%$max_sim_jobs\" -o $output_base_name\%I -e $error_base_name\%I 'many_fasta.pl $fasta_args $bsub_fofn $genome_base_name\${LSB_JOBINDEX} $fasta_result_base_dir\${LSB_JOBINDEX}'") == 0 or die "bsub failed: $?, $!\n";

# -K, which blocks until job completion, can't be used on a job array, so we use it on a trivial job dependent on the array completing:
$logger->debug("bsub -q normal -R \"select[mem > 1] rusage[mem=1]\" -w 'ended(bpa${unique}l1)' -K -J bpa${unique}l1_finish -o /dev/null -e /dev/null \"echo finished > /dev/null\"");
system("bsub -q normal -R \"select[mem > 1] rusage[mem=1]\" -w 'ended(bpa${unique}l1)' -K -J bpa${unique}l1_finish -o /dev/null -e /dev/null \"echo finished > /dev/null\"") == 0 or die "bsub -k failed: $?, $!\n";
	
# we're going to need a fofn file for each genome next, so open those now
my @fofn_filehandles;
for my $num (1..$genome_count) {
    my $fh;
    open($fh, ">$bsub_fofn$num") or die "error: could not open inhouse fofn file '$bsub_fofn$num'\n";
    push(@fofn_filehandles, $fh);
}
$logger->info("FASTA of all the above new sequences against each in-house genomes done");	


################################################################################
# get the top fasta hit (over certain thresholds) to each in-house genome;
# for any hits found extract the sequence from the in-house database and
# fasta it against the whole new database - if the top hit is the same as the
# original new_gene_id keep the result as a reciprocal hit.
#
# we also keep track of the best hit for each gene in case no reciprocal hit
# is found
my %top_hits;
my $inhouse_single_dir = "$unique_dir/inhouse";
foreach my $single_seq_file (@single_files) {
    $logger->debug($single_seq_file);
    # go through the fasta results and get the top hit details, extracting the
    # sequence from in-house genome for a fasta later
    my $single_result_name = basename($single_seq_file);
    my $id = $single_result_name;
    $single_result_name = get_hashed_dirs($single_result_name).'/'.$single_result_name;
    for my $num (1..$genome_count) {
	my $result_file = $fasta_result_base_dir.$num.'/'.$single_result_name;
	my $search = Bio::PSU::SearchFactory->make(-file => $result_file, -program => 'fasta');
			
	my @results = get_top_hit($search);
			
	# if we got a top hit extract the sequence and note the result
	if (@results == 10) {
	    my ($s_id, $score, $percent, $expect, $pc_query, $pc_subject, @others) = @results;
				
	    # store as the best hit if it is, and satisfies the top hit cutoffs
	    if ($expect <= $top_evalue_cutoff) {
		if (exists $best_hits{$id}) {
		    $score > ${$best_hits{$id}}[2] and $best_hits{$id} = [(($num - 1), $s_id, $score, $percent, $expect, @others)];
		} else {
		    $best_hits{$id} = [(($num - 1), $s_id, $score, $percent, $expect, @others)];
		}
	    }
	
  	    # Set up for reciprocal blast for those over the reciprocal cutoffs
	    if ($percent >= $percent_id_cutoff && $pc_query >= $percent_hit_length_cutoff && $pc_subject >= $percent_hit_length_cutoff) {
		$top_hits{$num}->{$id} = [($s_id, $score, $percent, $expect, @others)];
		my $seq_ref = $fastfastas{$genomes[$num - 1]}->get_seqs($s_id);
		my $sequence = ${$seq_ref}{$s_id};
					
		my $inhouse_seq_file_dir = create_hashed_dirs($inhouse_single_dir, $s_id);
		my $inhouse_seq_file = $inhouse_seq_file_dir.'/'.$s_id;
					
		open(SINGLE, ">$inhouse_seq_file") or die "error: could not write to '$inhouse_seq_file'\n";
		print SINGLE ">$s_id\n$sequence\n";
		close(SINGLE);
			
		# store this sequence location in the fofn file for this genome
		$logger->debug($inhouse_seq_file . " " . $s_id);
		print {$fofn_filehandles[$num - 1]} "$inhouse_seq_file\n";
	    }
	}
			
	# delete the result file
	unlink($result_file);
    }
}
# close the fofn files for each genome
@fofn_filehandles = ();

$logger->info("RECIPROCAL FASTA hits done");


################################################################################
# now do a fasta of the extracted inhouse sequences against the whole
# new genome db, keeping top hits that match back to $id as reciprocals
my @tops = keys %top_hits;
if (@tops > 0) {
    my $range = join(',', @tops);
		
    # the fasta - the sequences extracted from each genome vs the new genome
    $logger->debug("bsub -q basement -J bpa${unique}l2\"[$range]\%$max_sim_jobs\" -o $output_base_name\%I -e $error_base_name\%I 'many_fasta.pl $fasta_args $bsub_fofn\${LSB_JOBINDEX} $new_genome $fasta_result_base_dir\${LSB_JOBINDEX}'");
    system("bsub -q basement -J bpa${unique}l2\"[$range]\%$max_sim_jobs\" -o $output_base_name\%I -e $error_base_name\%I 'many_fasta.pl $fasta_args $bsub_fofn\${LSB_JOBINDEX} $new_genome $fasta_result_base_dir\${LSB_JOBINDEX}'") == 0 or die "bsub failed: $?, $!\n";
    $logger->debug("bsub -q normal -R \"select[mem > 1] rusage[mem=1]\" -w 'ended(bpa${unique}l2)' -K -J bpa${unique}l2_finish -o /dev/null -e /dev/null \"echo finished > /dev/null\"");
    system("bsub -q normal -R \"select[mem > 1] rusage[mem=1]\" -w 'ended(bpa${unique}l2)' -K -J bpa${unique}l2_finish -o /dev/null -e /dev/null \"echo finished > /dev/null\"") == 0 or die "bsub -k failed: $?, $!\n";
		
    # go through the fasta results and get the top hit details, see if we have a reciprocal and store
    for my $num (@tops) {
	while (my ($id, $result_ref) = each %{$top_hits{$num}}) {
	    my $path_to_orig_sid = get_hashed_dirs(${$result_ref}[0]).'/'.${$result_ref}[0];
	    my $result_file = $fasta_result_base_dir.$num.'/'.$path_to_orig_sid;
				
	    # check that the file exists (old bug now fixed?)
	    unless (-s $result_file) {
		warn "warning: fasta result file '$result_file' for ${$result_ref}[0] vs $id didn't exist!\n";
		next;
	    }
				
	    my $search = Bio::PSU::SearchFactory->make(-file => $result_file, -program => 'fasta');
       	    my ($s_id) = get_top_hit($search);
				
	    if ($s_id && $s_id eq $id) {
		push(@{$reciprocal_hits{$id}}, [(($num - 1), @{$result_ref})]);
	    }
	}
    }
}
}

# get the names of all the single sequence files
unless (@single_files > 0) {
    open(FOFN, $bsub_fofn) or die "error: could not open file '$bsub_fofn'\n";
    @single_files = <FOFN>;
    chomp(@single_files);
    close(FOFN);
}
$logger->info("FASTA of all extracted in-house sequences against new genome done");


################################################################################
# run HAMAP on all the single sequences in a bsub array

# without JobArray
$logger->debug("HAMAP without JobArray begin");
foreach my $hamap_input_file (@single_files) {
    $logger->debug($hamap_input_file);
    my $hamap_output_file = $hamap_input_file."_out";
    my $range = join(',', @single_files);
    $logger->debug("bsub -q basement -J bpa${unique}l1\"[$range]\%$max_sim_jobs\" -o $output_base_name.\%J.\%I -e $error_base_name.\%J.\%I '$hamap_bin $hamap_args $hamap_input_file $hamap_profile > $hamap_output_file'");
    # system("bsub -q basement -J bpa${unique}l1\"[1-$genome_count]\%$max_sim_jobs\" -o $output_base_name\%I -e $error_base_name\%I 'many_fasta.pl $fasta_args $bsub_fofn $genome_base_name\${LSB_JOBINDEX} $fasta_result_base_dir\${LSB_JOBINDEX}'") == 0 or die "bsub failed: $?, $!\n";

    # -K, which blocks until job completion, can't be used on a job array, so we use it on a trivial job dependent on the array completing:
    $logger->debug("bsub -q normal -R \"select[mem > 1] rusage[mem=1]\" -w 'ended(bpa${unique}l1)' -K -J bpa${unique}l1_finish -o /dev/null -e /dev/null \"echo finished > /dev/null\"");
    #system("bsub -q normal -R \"select[mem > 1] rusage[mem=1]\" -w 'ended(bpa${unique}l1)' -K -J bpa${unique}l1_finish -o /dev/null -e /dev/null \"echo finished > /dev/null\"") == 0 or die "bsub -k failed: $?, $!\n";

}
$logger->debug("HAMAP without JobArray end");

# with JobArray
$logger->debug("HAMAP with JobArray begin");
my $hamap_results_file = "${bsub_fofn}_hamap_results";
my $job_array = LSF::PSU::JobArray->new( bsub_options => "-q basement -R \"select[mem > 500] rusage[mem=500]\"",
                                         script_options => "-b $unique $input_options",
                                         jobs => scalar(@single_files),
                                         output_file => $hamap_results_file,
                                         error_file => $jobarray_dir . "/jobarray_err",
                                         work_dir => $jobarray_dir);

# submit the job array and get an output filename we can use
$job_array->bsub();
my $output_file = $job_array->output_file();

# since there are as many jobs as files, the array_fraction will be
# a single element and we actually only loop once here
foreach my $input_file ($job_array->array_fraction(\@single_files)) {
    $logger->debug("HAMAP in loop");
    $logger->debug("$hamap_bin $hamap_args $input_file $hamap_profile > $output_file");
    system("$hamap_bin $hamap_args $input_file $hamap_profile > $output_file") == 0 or die "hamap failed: $?, $!\n";
}

# stop array jobs now
$job_array->exit();

# cat the hamap output together
$job_array->cat();
$logger->debug("HAMAP with JobArray end");

# store the hamap results in a hash
my %hamap_hits;
if (-e $hamap_results_file) {
    open(HAMAP, $hamap_results_file) or die "error: could not open hamap results file '$hamap_results_file'\n";
	
    # make a hash of the hamap motif->info file
    open(ANNOT, $hamap_annotation) or die "error: could not open hamap annotation file '$hamap_annotation'\n";
    my %annotation;
    while (<ANNOT>) {
	my ($motif_id, $protein_name, $protein_desc) = $_ =~ /^\S+\s+(\S+)\s+\d+\s+\S+\s+(\S+)\s+(.+)/;
		
	# strip out the ec number, if present, and store it seperately
	$protein_desc =~ s/ \(EC (\S+)\)//;
		
	$annotation{$motif_id} = [$protein_name, $protein_desc, $1];
    }
    close(ANNOT);
	
    while (<HAMAP>) {
	my ($single_id, $start, $end, $motif_id, $score, $raw_score, $level) = $_ =~ /^>(.+)\/(\d+)-(\d+) motif=(.+)\|\S+ norm_score=(\d+\.?\d*) raw_score=(\d+) level=(\d)/;
		
	exists $annotation{$motif_id} or (warn "warning: hamap motif $motif_id not known, unable to add hamap annotation for $single_id!\n" and next);
		
	# store the highest scoring hamap hit for this protein
	if (exists $hamap_hits{$single_id}) {
	    if ($score > ${$hamap_hits{$single_id}}[2]) {
		$hamap_hits{$single_id} = [$start, $end, $score, $raw_score, $level, @{$annotation{$motif_id}}];
	    }
	} else {
	    $hamap_hits{$single_id} = [$start, $end, $score, $raw_score, $level, @{$annotation{$motif_id}}];
	}
    }
    close(HAMAP);
    unlink($hamap_results_file);
}

# delete all the single sequences
#$logger->debug("rm -fr $single_file_base_dir");
#system("rm -fr $single_file_base_dir");
#unlink($bsub_fofn);

$logger->info("HAMAP done");

################################################################################
# add orthologue links to the new embl/tab file for each gene_id with
# reciprocal_hits, transfering annotation for the highest scoring one
# open the file with bioperl, prefering hamap results
foreach my $feature_file (@embl_files, @tab_files) {
    my $in_tab = Bio::PSU::SeqFactory->make(-file => $feature_file, -format => 'embl');
    my $out_tab = Bio::PSU::SeqFactory->make(-file => ">$output_name", -format => 'embl');
	
    while (my $seq = $in_tab->next_seq()) {
	foreach my $feat ($seq->features) {
	    if ($feat->key eq 'CDS') {
		# don't do pseudo or partials
		($feat->qexists('pseudo') || $feat->qexists('partial')) and next;
				
		# get the systematic id
		my $sys_id;
		$feat->qexists($qualifier) or next;
		foreach my $possible_sys ($feat->qvalues($qualifier)) {
		    if ($possible_sys =~ /$regex/) {
			$sys_id = $1;
			last;
		    }
		}
		$sys_id or next;
				
		my @results;
		my $recip_mode = 0;
		if (exists $reciprocal_hits{$sys_id}) {
		    # get the reciprocal results sorted by score
		    @results = @{$reciprocal_hits{$sys_id}};
		    @results = sort { ${$b}[2] <=> ${$a}[2] } @results;
					
		    $recip_mode = 1;
		} elsif (exists $best_hits{$sys_id}) {
		    # get the best result
		    @results = ($best_hits{$sys_id});
		} elsif (! exists $hamap_hits{$sys_id}) {
		    next;
		}
				
		my $hamap_mode = 0;
		if (exists $hamap_hits{$sys_id}) {
		    # does it have an ec number?
		    defined ${$hamap_hits{$sys_id}}[7] ? ($hamap_mode = 2) : ($hamap_mode = 1);
		}
				
		# add ortholog links, similarity data to this feature for each result
		my $first = 1;
		foreach my $result_ref (@results) {
		    my ($genome_index, $s_id, $score, $pecent_id, $e_value, $subject_length, $overlap, $query_overlap, $subject_overlap) = @{$result_ref};
		    my $genome = $genomes[$genome_index];
					
		    # only an ortholog if it was a reciprocal best match
		    $feat->qadd('ortholog', "GeneDB_$genome:$s_id") if $recip_mode;
					
		    # always interested in similarity data
		    my $desc = $fastfastas{$genome}->desc($s_id);
		    my ($gene_name, $ec_number, $classification, $go, $colour, @product) = split(" ", $desc);
					
		    # debug sanity check
		    unless ($gene_name && $ec_number && $classification && $go && (defined $colour) && @product > 0) {
		        warn "something missing from the desc for $genome, $s_id : $desc\n";
			next;
		    }
					
		    my $product = "@product";
		    my $product_for_sim = $product;
		    $product_for_sim =~ s/;/\\;/g;
		    $feat->qadd('similarity', "fasta; $s_id; $genome; $product_for_sim; length $subject_length aa; id=$pecent_id\%; ; E()=$e_value; $overlap aa overlap; query $query_overlap aa; subject $subject_overlap aa");
					
		    # for the highest scorer also transfer annotation
		    if ($first) {
			# if it wasn't a reciprocal or hamap hit we only transfer colour and product
			if ($hamap_mode) {
			    $feat->qadd('primary_name', ${$hamap_hits{$sys_id}}[5]);
			}
			if ($recip_mode) {
			    unless ($hamap_mode) {
				unless($gene_name eq '-') {
				    $feat->qadd('primary_name', $gene_name);
				}
			    }
			    $hamap_mode == 2 and $ec_number = ${$hamap_hits{$sys_id}}[7];
			    unless($ec_number eq '-') {
				$feat->qadd('db_xref', 'EC:'.$ec_number);
				$feat->qadd('EC_number',  $ec_number);
			    }
			    unless($classification eq '-') {
				$feat->qadd('class', $classification);
			    }
			    unless($go eq '-') {
				$feat->qadd('GO', " ; GOid=$go; ; ; ; ; date=$date");
			    }
			}
			(defined $colour) or (warn "no colour!\n" and next);
			unless($colour eq '-') {
			    $feat->qadd('colour', $colour);
			}
			unless($product eq '-') {
			    unless ($recip_mode && $product !~ /^putative/) {
				# add 'putative' to the start of the product line for non-reciprocals
				$product = 'putative '.$product;
			    }
			    $hamap_mode and $product = ${$hamap_hits{$sys_id}}[6];
			    $feat->qadd('product', $product);
			}
						
			# say this is an automatic annotation and that it came
			# from this particular ortholog link
			if($hamap_mode) {
			    my $status = " ; evidence=IEA; db_xref=$genome:$s_id; date=$date; method=automatic:HAMAP";
			    $feat->qadd('status', $status);
			} else {
			    my $status = " ; evidence=IEA; db_xref=$genome:$s_id; date=$date; method=automatic:";
			    $recip_mode ? ($status .= "reciprocal fasta") : ($status .= "top fasta hit");
			    $feat->qadd('status', $status);
			}
						
			$first = 0;
		    }
		}
				
		# add a similarity and status line for the hamap result
		if ($hamap_mode) {
		    my $product_for_sim = ${$hamap_hits{$sys_id}}[6];
		    $product_for_sim =~ s/;/\\;/g;
		    $feat->qadd('similarity', "hamap; ; ; $product_for_sim; ; ; ; ; score=${$hamap_hits{$sys_id}}[2]; ; ; ; ${$hamap_hits{$sys_id}}[3]; ${$hamap_hits{$sys_id}}[0]-${$hamap_hits{$sys_id}}[1]");
					
		    # there might have been no fasta results at all so transfer gene_name and product
		    if (@results == 0) {
			$feat->qadd('primary_name', ${$hamap_hits{$sys_id}}[5]);
			$feat->qadd('product', ${$hamap_hits{$sys_id}}[6]);
		    }
			
		    # add the status line (hamap level translates to chado level by adding 1)
		    $feat->qadd('status', "level=".(${$hamap_hits{$sys_id}}[4] + 1)."; evidence=IEA; ; date=$date; method=automatic:hamap");
		}
	    }
	}
		
	# write out the seq with its modified features to output file
	$out_tab->write_seq($seq);
    }
}
$logger->info("Orthologue links done");


################################################################################
# copy out any errors
my $count = 0;
foreach my $genome (@genomes) {
    $count++;
    my $part_file = $error_base_name.$count;
    -s $part_file or next;
    open(ERRORS, ">>$error_file") or die "error: could not open error file '$error_file'\n";
    open(PART, "$part_file") or die "error: could not open part error file '$part_file'\n";
    while (<PART>) {
    	print ERRORS;
    }
    close(PART);
    close(ERRORS);
}
$logger->info("Errors copied to " . $error_base_name);


################################################################################
# delete all the working files
#system("rm -fr $unique_dir 2> /dev/null");
# it empties the dir but doesn't delete it?!... try again?
#sleep(5);
#system("rm -fr $unique_dir 2> /dev/null");

#print "...Deleted all working files...\n";
$logger->info("END");

exit;


################################################################################
# Get the top hit over threshold (based on code from reciprocal_fasta.pl)
# returns (s_id, score)
################################################################################
sub get_top_hit {
    my $search = shift;
	
    while (my $result = $search->next_result) {
        my @hits;
		
        while (my $hit = $result->next_hit) {
            push(@hits, $hit);
        }
		
        if (@hits > 0){
            # Sort the hits by score, percent
            @hits = sort { $a->sw_score <=> $b->sw_score || $a->percent <=> $b->percent } @hits;  
            my $top = pop(@hits);
			
            my $pc_query   = ($top->overlap / $top->q_len) * 100;
            my $pc_subject = ($top->overlap / $top->s_len) * 100;
 			
	    return ($top->s_id, $top->sw_score, $top->percent, $top->expect, $pc_query, $pc_subject, $top->s_len, $top->overlap, $top->q_begin.'-'.$top->q_end, $top->s_begin.'-'.$top->s_end);
        }
    }
    return ();
}


################################################################################
# Check to see if there are any bpa locks
################################################################################
sub check_locks {
    opendir(DIR, $dbs_dir) or die "error: could not open directory '$dbs_dir'\n";
    foreach my $file (readdir(DIR)) {
	-d $dbs_dir.$file or next;
	$file =~ /lock\d+/ and (return 1);
    }
    return 0;
}


################################################################################
#
################################################################################
sub catch_ctrlc {
	exit;
}


################################################################################
#
################################################################################
=head2 get_log_fn

Return Logfilename for this programm. It is the basename of this programm
with a ".log" ending.

=cut

sub get_log_fn {
    use File::Basename;
    return sprintf "%s.log", basename( $0, '.pl' );
}


END {
    if ($clear_locks) {
        rmdir($lock_dir.$unique);
	rmdir($lock_dir) unless check_locks();
    }
}

