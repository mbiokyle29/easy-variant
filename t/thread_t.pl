#!/usr/bin/perl
 
use strict;
use warnings;
 
use threads 1.39;
use threads::shared;
use Thread::Queue;
 
### Global Variables ###
 
# Maximum working threads
my $MAX_THREADS = 10;
 
# Flag to inform all threads that application is terminating
my $TERM :shared = 0;
 
# Threads add their ID to this queue when they are ready for work
# Also, when app terminates a -1 is added to this queue
my $IDLE_QUEUE = Thread::Queue->new();
 
 
### Signal Handling ###
 
# Gracefully terminate application on ^C or command line 'kill'
$SIG{'INT'} = $SIG{'TERM'} =
    sub {
        print(">>> Terminating <<<\n");
        $TERM = 1;
        # Add -1 to head of idle queue to signal termination
        $IDLE_QUEUE->insert(0, -1);
    };
 
 
### Main Processing Section ###
MAIN:
{
    ### INITIALIZE ###
 
    # Thread work queues referenced by thread ID
    my %work_queues;
 
    # Create the thread pool
    for (1..$MAX_THREADS) {
        # Create a work queue for a thread
        my $work_q = Thread::Queue->new();
 
        # Create the thread, and give it the work queue
        my $thr = threads->create('worker', $work_q);
 
        # Remember the thread's work queue
        $work_queues{$thr->tid()} = $work_q;
    }
 
 
    ### DO WORK ###
 
    # Manage the thread pool until signalled to terminate
    while (! $TERM) {
        # Wait for an available thread
        my $tid = $IDLE_QUEUE->dequeue();
 
        # Check for termination condition
        last if ($tid < 0);
 
        # Give the thread some work to do
        my $work = 5 + int(rand(10));
        $work_queues{$tid}->enqueue($work);
    }
 
 
    ### CLEANING UP ###
 
    # Signal all threads that there is no more work
    $work_queues{$_}->enqueue(-1) foreach keys(%work_queues);
 
    # Wait for all the threads to finish
    $_->join() foreach threads->list();
}
 
print("Done\n");
exit(0);
 
 
### Thread Entry Point Subroutines ###
 
# A worker thread
sub worker
{
    my ($work_q) = @_;
 
    # This thread's ID
    my $tid = threads->tid();
 
    # Work loop
    do {
        # Indicate that were are ready to do work
        printf("Idle     -> %2d\n", $tid);
        $IDLE_QUEUE->enqueue($tid);
 
        # Wait for work from the queue
        my $work = $work_q->dequeue();
 
        # If no more work, exit
        last if ($work < 0);
 
        # Do some work while monitoring $TERM
        printf("            %2d <- Working\n", $tid);
        while (($work > 0) && ! $TERM) {
            $work -= sleep($work);
        }
 
        # Loop back to idle state if not told to terminate
    } while (! $TERM);
 
    # All done
    printf("Finished -> %2d\n", $tid);
}
 
 