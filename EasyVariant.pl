#!/usr/bin/perl
#################################################################################################################1
# EasyVariant.pl --sam samfile --reference fasta/fa 
# This 'script' is a miniature, basic variant caller, it processes a sam file containing matched alignments
# Each SAM read is formed into a Moose object, containing Sequence objects and Cigar objects
# The SAM is then manipulated, and the Cigar string is parsed to determine what exactly the alignment determined
# The number of A,T,C, and G is tracked for each base position, along with X for deletions and I for insertions
# An insertion at what would be the 'new' base 10 is counted at the old base 9 (the base preceding the insertion)
# The optional arguments allow for a decent amount of fine tuning, including:
# -choosing specific ranges to only include
# -ignoreing a set of ranges
# -designating a set of ranges as repeat regions (For coverage statistics)
# -getting not only the variant calls, but the ATCGXI tally for each call
#
# Notes|Warnings|Disclaimer
# 	Using the --start --end flags, with the other range functions may cause unintended behvaior
#	best to only use one (start/end) or (ignore and repeat)
# 
#	Be sure to set use lib to the correct dir for the Moose modules
#
#	Range strings must be entered as int-int,int-int !!!
#	
#	The Genome package is backed by mysql, and might take some configuration
#################################################################################################################
use warnings;
use strict;
use lib "/home/kyle/lab/easy-variant/lib/";
use Genome;
use Sam;
use Getopt::Long;
use File::Slurp;
use feature qw(say switch);

# Predefine args for GetOpts
my ($start_pos, $end_pos, $sam_file, $reference);
my ($min_depth, $indel_ratio, $output);
my ($repeat_ranges, $ignore_ranges, $all_bases);

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
	"repeat-ranges=s" => \$repeat_ranges,
	"ignore-ranges=s" => \$ignore_ranges,
	"all-bases=i" => \$all_bases,
);

# Check for required args
unless($sam_file && $reference)
{
	say "--reference and --sam flags are required";
	die "EasyVariant.pl --reference ref --sam sam_file";
}

# Genome wide alignment results (This is where the magic happens!)
my %master_alignment;

# Hashes for optional ignore and repeat regions key = start, val = end
# If the flages were given at command line, parse the strings into hashes
my ($ignore_hash, $repeat_hash);
if($ignore_ranges)
{
	$ignore_hash = parse_ranges($ignore_ranges);
}
if($repeat_ranges)
{
	$repeat_hash = parse_ranges($repeat_ranges);
}

# 'Build' the genome, hopefully it's already in mysql
my $genome = Genome->new( fasta => $reference );

# Set defualt values if non were defined from command line
unless($end_pos) { $end_pos = $genome->length; }
unless($start_pos) { $start_pos = 1; }
unless($min_depth) { $min_depth = 10; }
unless($indel_ratio) { $indel_ratio = .5; }
unless($all_bases) {$all_bases = 0; }

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

	# Check if sam falls completely within the chosen range
	next if($end_pos < $sam->cigar->start_pos);
	next if($start_pos > $sam->cigar->end_pos);

	#Check if sam falls completely within an ignore range
	if($ignore_hash)
	{
		next if ignore_sam($sam->cigar->start_pos, $sam->cigar->end_pos, $ignore_hash);
	}

	# Pointers to walk the Sam sequence and genome sequence
	my $reference_pointer = $sam->cigar->start_pos;
	my $read_pointer = 1;
	
	# The cigar stack, from 1M1D3M to MDMMM
	my @cigar_stack = split(//,$sam->cigar->stack);

	# Flags
	my $first = 1;
	my $soft_clipping = 0;
	my $inserting = 0;
	my $start =  $sam->cigar->start_pos;
	my $end = $sam->cigar->end_pos;

	# Iterate through each 'Letter' of the cigar stack
	foreach my $cigar (@cigar_stack)
	{		
		given($cigar)
		{
			# To handle S, or soft clipping values that appear only at the start and end of the cigar string
			# We do only one increment for a run of S's, Dont increment pointers b/c seq was clipped off
			when("S")
			{
				$inserting = 0;

				# Handle a starting softclip, increment the insertion for
				# The postion one before the start of the cigar alignment
				if(!$soft_clipping && $first)
				{
					# Set the nessecary flags
					$soft_clipping = 1;
					$first = 0;
					
					# If we are at base number 1, increment the I at last positions
					if($reference_pointer == 1)
					{
						$master_alignment{$genome->length}{"I"} += 1;
					}
					else
					{
						$master_alignment{$reference_pointer-1}{"I"} += 1;
					}
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

			# When an M is popped we have a match/mismatch
			# Increment the counter for what ever base is at the position of sam read pointer
			# Increment both counters
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
			
			# When an I is poped we have an insertion
			# Similar to the S, we just increment the I atprior base (check if its the first)
			# and increment the read pointer
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

#####################
# SNP calling Block #
#####################

## Printing Args
my $genome_name = $genome->name;
my $genome_length = $genome->length;
unless ($output){ $output = $sam_file."output"; }
my $stat_file;

## INFO HEADER
open my $variants, ">", $output;
say $variants "EasyVariant caller results for $sam_file file against $genome_name";
say $variants "With: \n\t depth cutoff = $min_depth \n\t call cutoff = $indel_ratio\n";

## If all_bases flag was set initalize the file handle and print header
my $all_base;
if($all_bases)
{
	open $all_base, ">", $output.".perbase";
	say $all_base "EasyVariant caller base results for $sam_file file against $genome_name";
	say $all_base "With: \n\t depth cutoff = $min_depth \n\t call cutoff = $indel_ratio\n";
}

# Counters to track number of positions that had the right depth
# Conditionally set the unique/repeat counters if enabled
my $passed = 0;
my ($unique_passed, $repeat_passed);
if($repeat_hash) 
{ 
	$repeat_passed = 0;
	$unique_passed = 0; 
}

# Get the whole reference sequence, total the tallys in the master alignment and compare
my $ref_seq = $genome->seq;
foreach my $key (sort({ $a <=> $b} keys(%master_alignment)))
{
	# First make sure base is not in the ignore range
	if($ignore_hash)
	{
		next if in_range($key, $ignore_hash);
	}

	my $wildtype = uc(substr($ref_seq, ($key-1), 1));
	
	# Set up things for the tallying
	my $most_base;
	my $denominator = 0; 
	my $depth = 0; 
	my $most = 0;
	foreach my $call (keys($master_alignment{$key}))
	{
		# The denominator is A+T+C+G+X
		if($call ne "I")
		{
			$denominator += $master_alignment{$key}{$call};
		}

		# The depth is just A+T+C+G
		if($call ne "X" && $call ne "I")
		{
			$depth += $master_alignment{$key}{$call};
		}
		if($master_alignment{$key}{$call} > $most)
		{
			$most = $master_alignment{$key}{$call};
			$most_base = $call;
		}
	}
	# Depth check
	next if($depth < $min_depth);
	$passed++;

	# Check repeat / Unique counters and increment
	if($repeat_hash)
	{
		if(in_range($key, $repeat_hash)) { $repeat_passed++; }
		else { $unique_passed++; }
	}

	# Variant and ratio of variance check
	next if($wildtype eq uc($most_base));
	next if($most/$denominator < $indel_ratio);	

	# If we have made it this far, we have a variant!
	say $variants "$key: $wildtype --> $most_base";
	
	# Conditionally print detailed base information
	if($all_bases)
	{
		say $all_base "At position: $key";
		my $ref = $master_alignment{$key};
		for my $base (keys(%$ref))
		{
			say $all_base "$base -> $$ref{$base}";
		}
		say $all_base "";
	}
}

###########################################
# Calculate the Denominators for coverage #
###########################################
say $variants "################# COVERAGE INFORMATION ##################\n";

# The total bases for coverage, can be less then genome length if there are ignore ranges
my $coverage_denominator = $genome_length;
if($ignore_hash)
{
	# Convert entered string into hash(start)->end
	my $ignore_count = count_ranges($ignore_hash);
	$coverage_denominator -= $ignore_count;
	$ignore_ranges =~ s/,/, /g;
	say $variants "\nDropping $ignore_count bases from total coverage, becuase --ignore $ignore_ranges";
}

# Depth coverage of all bases, minus any ignored region
my $depth_percent = ($passed/$coverage_denominator)*100; # TODO DONT KEEP TRACK
say $variants "\n$passed bases had a depth of $min_depth or more, out of the total $coverage_denominator, ($depth_percent%)";

# Now handle depth coverage for unique regions
if($repeat_hash)
{
	my $repeat_count = count_ranges($repeat_hash);
	$coverage_denominator -= $repeat_count;
	my $unique_percent = ($unique_passed/$coverage_denominator)*100;
	my $repeat_percent = ($repeat_passed/$repeat_count)*100;

	# Make it look nice
	$repeat_ranges =~ s/,/, /g;
	say $variants "Using the following repeat regions:  $repeat_ranges\n";
	say $variants "$unique_passed unique bases had a min depth of $min_depth, out of the unique total of $coverage_denominator, ($unique_percent%)";
	say $variants "$repeat_passed repeat bases had a min depth of $min_depth, out of the repeat total of $repeat_count, ($repeat_percent%)";
}

# Close all file handles
close $variants;
if($all_bases) { close $all_base; }

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

############################################
# Calculate the length of repeat range set #
############################################
sub count_ranges
{
	my $range_hash = shift;
	my $count = 0;
	foreach my $start (keys(%$range_hash))
	{
		my $end = $$range_hash{$start};
		# Sanity check
		unless($end < $start)
		{
			# Inclusive length of range
			$count+= (($end - $start)+1);
		}
	}
	return $count;
}

######################################################
# See if SAM falls entierly within the ignore ranges #
######################################################
sub ignore_sam
{
	my ($sam_start, $sam_end, $ignore_hash) = @_;
	my $ignore = 0;
	foreach my $start (keys(%$ignore_hash))
	{
		my $end = $$ignore_hash{$start};
		if($end >= $sam_end && $start <= $sam_start)
		{
			say "Ignoreing sam $sam_start -> $sam_end";
			$ignore = 1;
			last;
		}
	}
	return $ignore;
}

#################################################
# See if a given position is in a set of ranges #
#################################################
sub in_range
{
	my ($position, $ranges) = @_;
	my $in_ranges = 0;
	foreach my $start (keys(%$ranges))
	{
		my $end = $$ranges{$start};
		if($start <= $position && $position <= $end)
		{
			$in_ranges = 1;
			last;
		}
	}
	return $in_ranges;
}
