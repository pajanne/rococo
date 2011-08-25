#!/software/bin//perl
#
# Created on Jun 02, 2010
# by
# author: Anne Pajon (ap12)
# Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
#

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

use VRTrack::VRTrack;

# Help / documentation
my $help = 0;

# Database connection details
my $host = $ENV{'VRTRACK_HOST'};
my $port = $ENV{'VRTRACK_PORT'};
my $user = $ENV{'VRTRACK_RW_USER'};
my $password = $ENV{'VRTRACK_PASSWORD'};
my $db = 'pathogen_metahit_track';

# Lane name
my $lane_name = '';

# Print help if no argument given
pod2usage(-verbose => 1) if (@ARGV < 1);

# Command line arguments
GetOptions ("lane=s" => \$lane_name,
	    "host=s" => \$host,
	    "port=i" => \$port,
            "user=s" => \$user,
	    "password=s" => \$password,
	    "db=s" => \$db,
	    "help"=> \$help) or pod2usage(-verbose => 1);

# Print help when --help
pod2usage(-verbose => 1) if $help;

# Print help if --lane not given
pod2usage(-verbose => 1) if ($lane_name eq '');

# Connect to database
my $track = VRTrack::VRTrack->new({ host => $host,
	                            port => $port,
	                            user => $user,
	                            password => $password,
	                            database => $db,
	                          });

# Return hierarchy path of a given lane if defined
if ($track->hierarchy_path_of_lane_name($lane_name)) {
  print $track->hierarchy_path_of_lane_name($lane_name), "\n";
} else {
  print "undefined\n";
}

__END__

=head1 NAME

 get_hierarchy_path_from_lane.pl - returns the hierarchy path based on the lane name

=head1 SYNOPSIS

 get_hierarchy_path_from_lane.pl [--lane='name'] [Options]

 Options:
 --lane='name' : name of the lane (required)

 --db='pathogen_metahit_track' : database name
 --host='web-mii-sha'p : database host
 --port=3303 : database port
 --user='pathpipe_ro' : database user
 --password='' : database password

 --help : help message

=cut
