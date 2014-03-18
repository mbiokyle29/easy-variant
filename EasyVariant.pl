#!/usr/bin/perl
use lib "/home/kyle/lab/perlpipe/perl/moose/";
use Genome;
use Sam;
use Getopt::Long;
use File::Slurp;
use Data::Dumper;
use feature qw(say switch);

my ($start_pos, $end_pos, $sam_file, $reference);

GetOptions (
	"start=i" => \$start_pos,
	"end=i" => \$end_pos,
	"sam=s" => \$sam_file,
	"reference=s" => \$reference
);

unless($sam_file && $reference)
{
	say "--reference and --sam flags are required";
	die "EasyVariant.pl --reference ref --sam sam_file";
}

my $genome = Genome->new(
	fasta => $reference
);

unless($end_pos)
{
	$end_pos = $genome->length;
}

unless($start_pos)
{
	$start_pos = 1;
}

if($end_pos < $start_pos)
{
	say "Ending Position Cannot be less that starting position!";
	die "Please add valid values for --end and/or --start";
}

my @sam_lines = read_file($sam_file);
foreach my $sam_line (@sam_lines)
{
	next if ($sam_line =~ m/^@/);
	my $sam = Sam->new(raw_string => $sam_line);

	next if($end_pos < $sam->cigar->start_pos);
	next if($start_pos > $sam->cigar->end_pos);
}
