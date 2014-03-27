#!/usr/bin/perl
use Data::Printer;
use threads;
use threads::shared;
use Thread::Pool;
my @shared :shared;

my $t_pool = Thread::Pool->new(
	{
		do => 'align_sam',
		workers => 20,
	}
);

my $counter = 1;
while($counter < 10)
{
	$t_pool->job($counter, \@shared);
	$counter++;
}
$t_pool->shutdown;
p(@shared);

sub align_sam
{
	my $in = shift;
	my $ref = shift;
	$$ref[$in] = "Hello from $in";
	sleep 1;
}