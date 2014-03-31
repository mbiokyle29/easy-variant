#!/usr/bin/perl
use warnings;
use strict;
use lib "/home/kyle/lab/easy-variant/";
use Genome;
use Sam;
use Getopt::Long;
use File::Slurp;
use Data::Printer;
use feature qw(say switch);

my $range_string = shift;
my $file = shift;

my @ranges = split(",", $range_string);
my %ranges_hash;

my $ranged_coverage = 0;

foreach my $range (@ranges)
{
	$range =~ m/(\d+)-(\d+)/;
	my $lower = $1;
	my $upper = $2;
	$ranges_hash{$lower}=$upper;
}

p(%ranges_hash);

open my $fh, "<", $file;
while(<$fh>)
{
	my $pos = $_;
	chomp($pos);
	my $in_range = 0;
	foreach my $lower (keys(%ranges_hash))
	{
		my $upper = $ranges_hash{$lower};
		if($lower <= $pos && $pos <= $upper )
		{
			$in_range = 1;
		}
	}
	if($in_range)
	{
		$ranged_coverage++;
	}	
}

print $ranged_coverage;
