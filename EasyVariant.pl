#!/usr/bin/perl
use warnings;
use strict;
use lib "/home/kyle/lab/easy-variant/";
use Genome;
use Sam;
use Getopt::Long;
use File::Slurp;
use Data::Dumper;
use feature qw(say switch);

# Predefine args
my ($start_pos, $end_pos, $sam_file, $reference, $min_depth, $indel_ratio, $output, $stats, $repeat_string, $ignore_ranges;

GetOptions (
	#Required arguments
	"sam=s" => \$sam_file,
	"reference=s" => \$reference,
	# Optional Arguemts
	"start=i" => \$start_pos,
	"end=i" => \$end_pos,
	"depth=i" => \$min_depth,
	"indel=f" => \$indel_ratio,
	"output=s" => \$output,
	"stats=i" => \$stats,
	"repeat-ranges" => \$repeat_string,
	"ignore-ranges" => \$ignore_ranges,
);

# Genome wide alignment results (This is where the magic happens!)
my %master_alignment;

# Check for required args
unless($sam_file && $reference)
{
	say "--reference and --sam flags are required";
	die "EasyVariant.pl --reference ref --sam sam_file";
}

my $genome = Genome->new( fasta => $reference );

# Set defualt values if non were defined from command line
unless($end_pos) { $end_pos = $genome->length; }
unless($start_pos) { $start_pos = 1; }
unless($min_depth) { $min_depth = 10; }
unless($indel_ratio) { $indel_ratio = .5; }
unless($stats) { $stats = 0; }

# Invalid range case
if($end_pos < $start_pos)
{
	say "Ending Position Cannot be less that starting position!";
	die "Please add valid values for --end and/or --start";
}

#################################################################
# SAM processing: cigar+alignment -> genome wide alignment matrix
#################################################################
my @sam_lines = read_file($sam_file);
foreach my $sam_line (@sam_lines)
{
	# Skip header lines and unaligned reads
	next if($sam_line =~ m/^@/);
	next if($sam_line =~ m/^(\S+\s+){2}\*/);
	
	# Create sam object
	chomp($sam_line);
	my $sam = Sam->new(raw_string => $sam_line);

	next if($end_pos < $sam->cigar->start_pos);
	next if($start_pos > $sam->cigar->end_pos);

	my $reference_pointer = $sam->cigar->start_pos;
	my $read_pointer = 1;
	my @cigar_stack = split(//,$sam->cigar->stack);

	# Flags
	my $first = 1;
	my $soft_clipping = 0;
	my $inserting = 0;
	my $start =  $sam->cigar->start_pos;
	my $end = $sam->cigar->end_pos;

	foreach my $cigar (@cigar_stack)
	{		
		given($cigar)
		{
			# To handle S, or soft clipping values that appear only at the start and end of the cigar string
			when("S")
			{
				$inserting = 0;

				# Handle a starting softclip, increment the insertion for
				# The postion one before the start of the cigar alignment
				if(!$soft_clipping && $first)
				{
					$soft_clipping = 1;
					$master_alignment{$sam->cigar->start_pos}{"I"} += 1;
					$first = 0;
				}

				# Handle the ending softclip, increment the insertion for
				# The postion at the end of the cigar alignment
				elsif(!$soft_clipping && !$first)
				{
					$soft_clipping = 1;
					$master_alignment{$sam->cigar->end_pos}{"I"} += 1;
				}
			}

			# When a D is popped we have a deletion
			# Increment the deletion counter at this position in the reference
			# Incremenet the reference pointer
			# Do no increment the read pointer becuase this was a deletion (ie not in the read)
			when("D")
			{
				# Accouting
				$soft_clipping = 0;
				$inserting = 0;
				if($first) { $first = 0;}
				
				$master_alignment{$reference_pointer}{"X"} += 1;
				$reference_pointer++;
			}

			when("M")
			{
				# Accouting
				$soft_clipping = 0;
				$inserting = 0;
				if($first) { $first = 0;}

				$master_alignment{$reference_pointer}{$sam->sequence->base_at($read_pointer)} += 1;
				$reference_pointer++;
				$read_pointer++;
			}

			when("I")
			{
				# Accouting
				$soft_clipping = 0;
				if($first) { $first = 0;}

				if(!$inserting)
				{
					if($reference_pointer == 1)
					{
						$master_alignment{$genome->length}{"I"} += 1;
					}
					else
					{
						$master_alignment{$reference_pointer-1}{"I"} += 1;
					}
					$inserting = 1;

				}
				$read_pointer++;
			}
		}
	}
}

###################
# SNP calling Block
###################

## Printing Args
my $genome_name = $genome->name;
my $genome_length = $genome->length;
unless ($output){ $output = $sam_file."output"; }
my $stat_file;

## Stats FILE
if($stats)
{
	my $stat_out = $output.".stats";
	open $stat_file, ">", $stat_out;
}

## INFO HEADER
open my $variants, ">", $output;
say $variants "EasyVariant caller results for $output file against $genome_name";
say $variants "With: \n\t depth cutoff = $min_depth \n\t call cutoff = $indel_ratio\n";

open my $all_base, ">", $output.".perbase";
say $all_base "EasyVariant caller base results for $output file against $genome_name";
say $all_base "With: \n\t depth cutoff = $min_depth \n\t call cutoff = $indel_ratio\n";

# Counter to track number of positions that had the right depth
my $passed = 0;
my $passed_positions;

my $ref_seq = $genome->seq;
foreach my $key (sort({ $a <=> $b} keys(%master_alignment)))
{
	my $wildtype = uc(substr($ref_seq, ($key-1), 1));
	my $denominator = 0;
	my $depth = 0;
	my $most = 0;
	my $most_base;

	foreach my $call (keys($master_alignment{$key}))
	{
		if($call ne "I")
		{
			$denominator += $master_alignment{$key}{$call};
		}

		$depth += $master_alignment{$key}{$call};
		if($master_alignment{$key}{$call} > $most)
		{
			$most = $master_alignment{$key}{$call};
			$most_base = $call;
		}
	}
	# Depth check
	next if($depth < $min_depth);
	$passed++;
	
	if($stats)
	{
		my $concat = $key."\n";
		$passed_positions.=$concat;
	}
	
	# Variant and ratio of variance check
	next if($wildtype eq uc($most_base));
	next if($most/$denominator < $indel_ratio);	

	say $variants "$key: $wildtype --> $most_base";

	say $all_base "At position: $key";
	my $ref = $master_alignment{$key};
	for my $base (keys(%$ref))
	{
		say $all_base "$base -> $$ref{$base}";
	}
	say $all_base "";
}

##########################################
# Calculate the Denominator for coverage #
##########################################
my $coverage_denominator = $genome_length;
if($ignore_ranges)
{
	# Convert entered string into hash(start)->end
	my $ignore = parse_ranges($ignore_ranges);
	my $ignore_count = subtract_ignore($ignore, $genome_length);
	$coverage_denominator -= $ignore_count;
	say $variants "Dropping $ignore_count from coverage, becuase --ignore $ignore_ranges";
}

# Depth coverage of all bases, minus any ignored region
my $depth_percent = ($passed/$coverage_denominator)*100;
say $variants "\n$passed bases had a depth of $min_depth or more, out of the total $coverage_denominator, ($depth_percent%)";





# Close all file handles
close $variants;
close $all_base;
if($stats && $stat_file)
{
	print $stat_file $passed_positions;
	close $stat_file;
}



###############
# SUBROUTINES #
###############
########################################
# Convert Range String into range hash #
########################################
sub parse_ranges
{
	# Formated range range string x-y,a-b,c-d
	my $range_string = shift;
	my @ranges = split(",", $range_string);
	my %ranges_hash;
	foreach my $range (@ranges)
	{
		$range =~ m/(\d+)-(\d+)/;
		my $lower = $1;
		my $upper = $2;
		$ranges_hash{$lower}=$upper;
	}
	return \%ranges_hash;
}

##################################################
# Calculate length of genome minus ignore ranges #
##################################################
sub subtract_ignore
{
	my ($ignore, $genome_length) = @_;
	my $ignore_count = 0;
	for my $start (%$ignore)
	{
		my $end = $$ignore{$start};
		unless($end < $start)
		{
			$ignore_count += (($end - $start)+1);
		}	
	}
	return $genome_length-$ignore_count;
}