package Genome;
use Moose;
use namespace::autoclean;
use DBI;
use Data::Dumper;
use File::Slurp;

has 'length' => (
	is  => 'ro',
	isa => 'Int',
	builder => '_build_length',
	lazy => 1
);

has 'seq' => (
	is => 'rw',
	isa => 'Str',
	builder => '_build_seq',
	lazy => 1
);

has 'fasta' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'name' => (
	is => 'rw',
	isa => 'Str',
	builder => '_build_name',
);

sub _build_name
{
	my $self = shift;
	$self->fasta =~ m/(.*)\.fa.*$/;
	my $name = $1;
	$name =~ s/-/_/g;
	return $name;
}

sub _build_length
{
	my $self = shift;
	my $name = $self->name;
	my $countq = "SELECT position FROM $name ORDER BY position DESC LIMIT 1";
	my $dbh = DBI->connect('dbi:mysql:genomes','genomes','ebvHACK958$');
	my $sth = $dbh->prepare($countq);
	$sth->execute();
	
	my $return = $sth->fetchrow_arrayref();
	return $$return[0];
}

sub _build_seq
{
	my $self = shift;
	my $name = $self->name;
	my $ref_seq;
	my $query = "SELECT base FROM $name ORDER BY position ASC";
	my $dbh = DBI->connect('dbi:mysql:genomes','genomes','ebvHACK958$');
	my $sth = $dbh->prepare($query);
	$sth->execute;
	my $result = $sth->fetchall_arrayref;
	foreach my $base (@$result)
	{
		$ref_seq.=$$base[0];
	}
	return $ref_seq;
}

sub BUILD
{
	my $self = shift;
	my $name = $self->name;

	#Check if already exists save us some time
	my $dbh = DBI->connect('dbi:mysql:genomes','genomes','ebvHACK958$');
	my $table_check = "SHOW TABLES LIKE ?";
	my $check_sth = $dbh->prepare($table_check);
	$check_sth->execute($name);

	if($check_sth->fetchrow_arrayref)
	{
		# genome exists already
		return;
	}	

	# Get the reference
	my $fasta = $self->fasta;
	my @fasta_arr = read_file($fasta);

	my $ref_seq;       
	my $first = 1;
	foreach my $line (@fasta_arr)
	{
		chomp($line);
		if($first)
		{
			$first = 0;
			next;
		}

		elsif($line =~ m/^[^>\#\s]/)
		{
			$ref_seq.=$line;
		}
	}

	my $make_table = "CREATE TABLE $name (base VARCHAR(1), position BIGINT)";
	my $sth = $dbh->prepare($make_table);
	$sth->execute() or die "Error Creating table for genome $name $DBI::errstr\n";	

	my @bases = split(//,$ref_seq);
	my $counter = 1;
	my $insert = "INSERT INTO $name (base, position) VALUES(?, ?)";
	foreach my $base (@bases)
	{
		my $insert_h = $dbh->prepare($insert);
		$insert_h->execute($base, $counter);
		$counter++;
	}
}

sub base_at
{
	my $self = shift;
	my $pos = shift;
	$pos--;
	return substr($self->seq, $pos,1); 
}

__PACKAGE__->meta->make_immutable;
1;