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

# Genome wide alignment results (This is where the magic happens!)
my %master_alignment;

# Hashes for optional ignore and repeat regions key = start val = end
my $ignore_hash;
my $repeat_hash;
if($ignore_ranges)
{
	$ignore_hash = parse_ranges($ignore_ranges);
}
if($repeat_ranges)
{
	$repeat_hash = parse_ranges($repeat_ranges);
}


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

# Counter to track number of positions that had the right depth
my $passed = 0;
my $unique_passed = 0;
my $repeat_passed = 0;

my $ref_seq = $genome->seq;
foreach my $key (sort({ $a <=> $b} keys(%master_alignment)))
{
	# First make sure base is not in the ignore range
	if($ignore_hash)
	{
		next if in_range($key, $ignore_hash);
	}

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
		if(in_range($key, $repeat_hash))
		{
			$repeat_passed++;
		}
		else
		{
			$unique_passed++;
		}
	}

	# Variant and ratio of variance check
	next if($wildtype eq uc($most_base));
	next if($most/$denominator < $indel_ratio);	


	say $variants "$key: $wildtype --> $most_base";
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
my $coverage_denominator = $genome_length;
say $variants "################# COVERAGE INFORMATION ##################\n";
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
close $all_base;

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