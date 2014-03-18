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
	builder => '_build_end_pos'
);

has 'array' => (
	isa => 'ArrayRef',
	is => 'rw',
	auto_deref => 1,
	builder => '_build_array'
);

has 'stack' => (
	isa => 'Str',
	builder => '_build_stack',
	is => 'rw',
	lazy => 1
);

has 'length' => (
	isa => 'Int',
	is => 'rw',
	lazy => 1,
	builder => '_build_length'
);

has 'start_pos' => (
	isa => 'Int',
	builder => '_build_start',
	is => 'rw',
	lazy    => 1
);

has 'end_pos' => (
	isa => 'Int',
	builder => '_build_end',
	is => 'rw',
	lazy    => 1
);

sub _build_array
{
	my $self = shift;
	my $raw_string = $self->raw_string;
	my @array;

	while($raw_string =~ m/(\d+[SDMI])/g)
	{
		push(@array, $1);
	}

	return \@array; 
}

sub _build_stack
{
	my $self = shift;	
	my $cigar_stack;
	foreach my $chunk ($self->array)
	{
		my $code = chop($chunk);
		my $push = $code x $chunk;
		$cigar_stack .= $push;
	}
	return $cigar_stack;
}

sub _build_length
{
	my $self = shift;
	my $length = 0;
	foreach my $chunk ($self->array)
	{
		if(chop($chunk) =~ m/[MD]/)
		{
			$chunk =~ m/(\d+)/;
			$length += $chunk;
		}
	}
	return $length;
}

sub _build_end
{
	my $self = shift;
	if($self->length)
	{
		return ($self->start_pos+$self->length);
	}
	return -1;
}

1;
