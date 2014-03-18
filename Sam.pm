package Sam;
use Moose;
use Cigar;
use Sequence;
use namespace::autoclean;

has 'raw_string' => (
	isa => 'Str',
	required => 1,
	is => 'ro'
);

has 'sequence' => (
	isa => 'Sequence',
	is => 'rw'
);

has 'cigar' => (
	isa => 'Cigar',
	is => 'rw'
);

has 'start_pos' => (
	isa => 'Int',
	is => 'rw'
);

sub BUILD
{
	my $self = shift;
	my @sam_fields = split("\t", $self->raw_string);
	my $seq_string = $sam_fields[9];
	my $cigar_string = uc($sam_fields[5]);
	my $start_pos = $sam_fields[3];

	my $cigar = Cigar->new(
		raw_string => $cigar_string,
		start_pos => $start_pos
	);

	my $sequence = Sequence->new(
		seq => $seq_string
	);

	$self->sequence($sequence);
	$self->cigar($cigar);
	$self->start_pos($start_pos);
}

__PACKAGE__->meta->make_immutable;
1;