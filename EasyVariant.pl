#!/usr/bin/perl
use lib "/home/kyle/lab/easy-variant/";
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

my %master_alignment;

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

open my $output, ">", "/home/kyle/lab/results";
open my $inter, ">", "/home/kyle/lab/intermediate";

my @sam_lines = read_file($sam_file);
foreach my $sam_line (@sam_lines)
{
	next if ($sam_line =~ m/^@/);
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
	say $inter "$start  --  $end";
	my $counter = $sam->cigar->start_pos;

	say $inter "Before";
	while($counter <= $end)
	{
		if($master_alignment{$counter})
		{
			say $inter "$counter";
			say $inter Dumper %{ $master_alignment{$counter} };
		}
		$counter++;
	}

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

			# When a D is popped, 
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
				## TODO starting insertions can happen?!?!
				$soft_clipping = 0;
				if($first) { $first = 0;}

				if(!$inserting)
				{
					## TODO FIX FIX FIX 
					if($read_pointer = 1)
					{
						$master_alignment{$reference_pointer}{"I"} += 1;
					}
					$master_alignment{$reference_pointer-1}{"I"} += 1;
					$inserting = 1;
				}
				$read_pointer++;
			}
		}
	}
	say $inter "";
	say $inter $sam->cigar->raw_string;
	my $start =  $sam->cigar->start_pos;
	my $end = $sam->cigar->end_pos;
	say $inter "$start  --  $end";
	print $inter "sam seq: ";
	say $inter $sam->sequence->seq;
	print $inter "reference: ";
	my $counter = $sam->cigar->start_pos;
	while($counter <= $end)
	{
		my $base = $genome->base_at($counter);
		say $inter "$counter -> $base";
		$counter++;
	}

	$counter = $sam->cigar->start_pos;
	say $inter "After";
	while($counter <= $end)
	{
		say $inter "$counter";
		say $inter Dumper %{ $master_alignment{$counter} };
		$counter++;
	}

}

foreach my $key (sort( {$a <=> $b} keys(%master_alignment)))
{
	say $output "$key:";
	foreach my $base (keys(%{ $master_alignment{$key} }))
	{
		say $output "\t$base  => $master_alignment{$key}{$base}";
	}
}
close output;