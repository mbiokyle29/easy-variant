package Sequence;
use Moose;
use namespace::autoclean;

has 'length' => (
	is  => 'rw',
	isa => 'Int',
	builder => '_build_length',
	lazy => 1
);

has 'seq' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'seq_arr' => (
	is => 'rw',
	isa => 'ArrayRef',
);

sub BUILD
{
	my $self = shift;
	my @seq = split(//, $self->seq);
	$self->seq_arr(\@seq);
}

sub base_at
{
	my $self = shift;
	my $pos = shift;
	$pos--;
	return substr($self->seq, $pos,1); 
}

sub range
{
	my $self = shift;
	my $start = shift;
	$start--;
	my $end = shift;
	return substr($self->seq, $start, $end);
}

sub _build_length
{
	return length(shift->seq);
}

__PACKAGE__->meta->make_immutable;
1;