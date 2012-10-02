package NGS::GenBank;

use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;

sub new {
	my ($class, %args) = @_;

	my $self = bless {}, ref($class) || $class;
	$self->_init(\%args);
	return $self;
}

sub _init {
	my ($self, $args) = @_;

	my $gb = Bio::SeqIO->new(-file => $args->{file},
							 -format => 'genbank');
	
	$self->{_gb} = $gb;
	$self->{file} = $args->{file};
	
	my $seq = $gb->next_seq;
	$self->{_seq_obj} = $seq;
	$self->{seq_id} = $seq->id;
}

sub format_seq_id {
	my ($self, %args) = @_;
	
	my $seq = $self->{_seq_obj};

	#
	# formats are a colon delimited string and should be the same name as function names for seq features.
	# They are processed in the order given in the option
	# e.g. format => 'display_id:seq_version'
	# $seq->display_id
	# $seq->seq_version
	#
	my @formats = split(':', $args{format});

	# create new id using the format names as the function names and joining them with a "."
	my $new_id = join ".", map { $seq->$_ } @formats;
	$seq->id( $new_id );
}

sub get_all_features {
	my ($self) = @_;

	my $seq = $self->{_seq_obj};
	
	# get all features
	my @features = $seq->get_SeqFeatures;
		
	# assign the gff3 ID attribute and assign the sequence identifier in case it was modified
	map { 
		if( $_->has_tag('locus_tag') && defined($_->each_tag_value('locus_tag')) ) { $_->primary_id( $_->each_tag_value('locus_tag') ); $_->seq_id( $seq->id ); } 
	} @features;

	return \@features;
}

sub get_features_from_type {
	my ($self, %args) = @_;

	my $seq = $self->{_seq_obj};
	my $type = $args{feat_type};
	
	# get all features for the given type
	my @features = grep { $_->primary_tag eq "$type" } $seq->get_SeqFeatures;
	
	# assign the gff3 ID attribute and assign the sequence identifier in case it was modified
	map { $_->primary_id( $_->each_tag_value('locus_tag') ); $_->seq_id( $seq->id ); } @features;

	return \@features;
}

sub get_genes {
	my ($self) = @_;

	my $seq = $self->{_seq_obj};
	
	# get all rRNA features
	my @features = grep { $_->primary_tag eq 'genes' } $seq->get_SeqFeatures;
	
	# assign the gff3 ID attribute and assign the sequence identifier in case it was modified
	map { $_->primary_id( $_->each_tag_value('locus_tag') ); $_->seq_id( $seq->id ); } @features;

	return \@features;
}

sub get_CDSs {
	my ($self) = @_;

	my $seq = $self->{_seq_obj};
	
	# get all rRNA features
	my @features = grep { $_->primary_tag eq 'CDS' } $seq->get_SeqFeatures;
	
	# assign the gff3 ID attribute and assign the sequence identifier in case it was modified
	map { $_->primary_id( $_->each_tag_value('locus_tag') ); $_->seq_id( $seq->id ); } @features;

	return \@features;
}

sub get_rRNAs {
	my ($self) = @_;

	my $seq = $self->{_seq_obj};
	
	# get all rRNA features
	my @features = grep { $_->primary_tag eq 'rRNA' } $seq->get_SeqFeatures;
	
	# assign the gff3 ID attribute and assign the sequence identifier in case it was modified
	map { $_->primary_id( $_->each_tag_value('locus_tag') ); $_->seq_id( $seq->id ); } @features;

	return \@features;
}

1;
