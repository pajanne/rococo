use strict;
use warnings;

# Count the number of file in a directory
my $dir = 'test_root/Alistipes/shahii_WAL8301/GFIND/seq_chunks/';
my @files = <$dir/alistipes_shahii_wal8301_*.fasta>;
my $count = @files;
print $dir, "\n";
print $count, "\n";

my $chrom_regex = '(?-xism:^(?:\d+|X|Y)$)';
my $do_chrm = '^(?:\d+|X|Y)$';
my $chrom = 'scaffold00001';

my $chr_regex = '^(scaffold\d{5})';

# Chromosome distribution
#
print "chrom_regex : ", $chr_regex, "\n";
print "do_chrm     : ", $do_chrm, "\n";
print "chrom       : ", $chrom, "\n";
if ( $do_chrm && $chrom=~$chr_regex ) {
    print "In\n";
} else {
    print "Out\n";
}

my $chr_regex_empty = '';
if ($chr_regex_empty) {
   print "Exists\n";
} else {
   print "Does not exist\n";
}

