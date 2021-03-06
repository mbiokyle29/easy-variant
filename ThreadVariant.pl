#!/usr/bin/perl
###################
## THREAD SET UP ##
###################
use warnings;
use strict;
use threads;
use threads::shared;
use Thread::Queue;
use lib "/home/kyle/lab/easy-variant/lib/";
use Sam;

my %master_alignment :shared;

# QUEUES #
my $work_queue = Thread::Queue->new();
my @results_queues :shared;
my @results_threads;

my $keep_working :shared = 1;
my $keep_results :shared = 1;

# WORKERS #
my @worker_threads;

for my $index (0..4)
{
    push(@results_queues, Thread::Queue->new());
    push(@worker_threads, threads->create('parse_sam', $index));
    push(@results_threads, threads->create('merge_results', $index));
}

###################
## MAIN SET UP   ##
###################
use Data::Printer;
use Genome;
use Getopt::Long;
use File::Slurp;
use feature qw(say switch);

# Predefine args
my ($start_pos, $end_pos, $sam_file, $reference, $output);
my ($indel_ratio, $min_depth);
GetOptions 
(
    #Required arguments
    "sam=s" => \$sam_file,
    "reference=s" => \$reference,
    
    # Optional Arguemts
    "start=i" => \$start_pos,
    "end=i" => \$end_pos,
    "output=s" => \$output,
    "depth=i" => \$min_depth,
    "indel=f" => \$indel_ratio,
);

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
    
    # Queue the SAM for processing
    $work_queue->enqueue($sam_line);
}
$keep_working = 0;
# Signal end of data, wait for threads to finish
$work_queue->end();
$_->join() for @worker_threads;
print "Work is done \n";
$keep_results = 0;
$_->join() for @results_threads;
p(%master_alignment);

###################
# SNP calling Block
###################

unless($indel_ratio) { $indel_ratio = .5; }
unless($min_depth) { $min_depth = 10; }
## Printing Args
my $genome_name = $genome->name;
my $genome_length = $genome->length;
unless ($output){ $output = $sam_file."output"; }

## INFO HEADER
open my $variants, ">", $output;
say $variants "EasyVariant caller results for $output file against $genome_name";
say $variants "With: \n\t depth cutoff = $min_depth \n\t call cutoff = $indel_ratio\n";

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

    foreach my $call (keys(% { $master_alignment{$key} }))
    {
        if($call ne "I")
        {
            $denominator += $master_alignment{$key}{$call};
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
        
    # Variant and ratio of variance check
    next if($wildtype eq uc($most_base));
    next if($most/$denominator < $indel_ratio); 

    say $variants "$key: $wildtype --> $most_base";
}

my $depth_percent = ($passed/$genome_length)*100;
say $variants "\n$passed bases had a depth of $min_depth or more, out of the total $genome_length, ($depth_percent%)";
close $variants;

sub merge_results
{
    my $index = shift;
    while($keep_results)
    {
    	next unless defined(my $local = $results_queues[$index]->dequeue_nb());
        foreach my $position (keys(%{ $local }))
        {
            foreach my $call (keys($$local{$position}))
            {
                lock(%master_alignment);
		unless(defined($master_alignment{$position}))
		{
		   print "Setting $position to new hash \n";
		   my %bases :shared =
		   (
		   	A => 0,
			T => 0,
			C => 0,
			G => 0,
			X => 0,
			I => 0,
		   );	
		   $master_alignment{$position} = \%bases;
		}
		$master_alignment{$position}{$call} += $$local{$position}{$call};
            }
        }
    }
}

sub parse_sam 
{
    my $index = shift;
    while($keep_working)
    {
        next unless defined(my $sam_line = $work_queue->dequeue_nb());
	my %local_alignment;
        my $sam = Sam->new(raw_string => $sam_line);
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
                        $local_alignment{$sam->cigar->start_pos}{"I"} += 1;
                        $first = 0;
                    }

                    # Handle the ending softclip, increment the insertion for
                    # The postion at the end of the cigar alignment
                    elsif(!$soft_clipping && !$first)
                    {
                        $soft_clipping = 1;
                        $local_alignment{$sam->cigar->end_pos}{"I"} += 1;
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

                    $local_alignment{$reference_pointer}{"X"} += 1;
                    $reference_pointer++;
                }

                when("M")
                {
                    # Accouting
                    $soft_clipping = 0;
                    $inserting = 0;
                    if($first) { $first = 0;}

                    $local_alignment{$reference_pointer}{$sam->sequence->base_at($read_pointer)} += 1;
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
                        $local_alignment{$reference_pointer-1}{"I"} += 1;
                        $inserting = 1;
                    }
                    $read_pointer++;
                }
            }
        }
        $results_queues[$index]->enqueue(\%local_alignment);
    }
    $results_queues[$index]->end();
}
