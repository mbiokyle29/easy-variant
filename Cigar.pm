package Cigar;
use Moose;
use namespace::autoclean;

has 'raw_string' => (
	isa => 'Str',
	required => 1,
	is => 'ro'
);

has 'start_pos' => (
	isa => 'Int',
	required => 1,
	is => 'ro'
);

has 'end_pos' => (
	isa => 'Int',
	is => 'rw',
	lazy => 1,
	builder => '_build_end'
);

has 'array' => (
	isa => 'ArrayRef',
	is => 'rw',
	auto_deref => 1,
);

has 'stack' => (
	isa => 'Str',
	is => 'rw',
);

has 'length' => (
	isa => 'Int',
	is => 'rw',
);

sub BUILD
{
	my $self = shift;	
	my $raw_string = $self->raw_string;
	
	my @array;
	while($raw_string =~ m/(\d+[SDMI])/g)
	{
		push(@array, $1);
	}
	
	$self->array(\@array);

	my $cigar_stack;
	foreach my $chunk ($self->array)
	{
		my $code = chop($chunk);
		my $push = $code x $chunk;
		$cigar_stack .= $push;
	}
	$self->stack($cigar_stack);

	my $length = 0;
	foreach my $chunk ($self->array)
	{
		if(chop($chunk) =~ m/[MD]/)
		{
			$chunk =~ m/(\d+)/;
			$length += $chunk;
		}
	}
	$self->length($length);
}

sub _build_end
{
	my $self = shift;
	if($self->length)
	{
		return ($self->start_pos+$self->length)-1;
	}
	return -1;
}

1;
