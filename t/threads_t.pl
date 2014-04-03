#!/usr/bin/perl 
use strict;
use warnings;
use threads;
use threads::shared;
use Thread::Queue;
use Data::Printer;
use lib '/home/kyle/Dev/easy-variant/';
use Sequence;
use Genome;

my $queue = Thread::Queue->new();    # A new empty queue
my %shared :shared;

my $max_theads = `nproc`;
$max_theads--;
my @threads;

for(1..$max_theads)
{
	my $thread = threads->create('worker');
	push(@threads, $thread);
}

my $counter = 0;
while($counter < 10)
{
	my $seq = Sequence->new(
		seq => "AAAAAAAAAAAA"
	);
	$queue->enqueue($seq);
	$counter++;
}

$queue->end();

$_->join for @threads;
p(%shared);

sub worker 
{
	my $thread_id = threads->tid();
	while(defined( my $sam = $queue->dequeue()))
	{
		my $base = $sam->base_at($thread_id);
		++$shared{$base};
	}
}