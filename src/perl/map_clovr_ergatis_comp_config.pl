#!/usr/bin/perl

=head1 NAME

map_clovr_ergatis_comp_config.pl - generates those variables that are different in clovr from that of ergatis into new pipeline_tmpl.config. 

=head1 SYNOPSIS

	./map_clovr_ergatis_comp_config.pl /full/path/to/clovr_component.config.file

	NOTE: if you are generating pipeline_tmpl.config file for the pipeline, recommended usage is:
	----	
	mv /path/to/pipeline_tmpl.config /path/to/pipeline_tmpl.config.bak &&
	touch /path/to/pipeline_directory/pipeline_tmpl.config &&
	for file in $(find /path/to/pipeline_directory/ -type f -name '*.config' -print) ;
	do
		perl map_clovr_ergatis_comp_config.pl $file 1>>/path/to/pipeline_directory/pipeline_tmpl.config
    							    2>>/path/to/some_error.log;
	done


=head1 AUTHOR
	
	Mahesh Vangala
	mvangala@som.umaryland.edu
	vangalamaheshh@gmail.com

=cut

use strict;
use warnings;

my $clovr_file = shift @ARGV;
my $ergatis_file = undef;
my $extension = undef;
my $component = undef;
my $dir = undef;

if( $clovr_file =~ /(.+)\/(.+)\.(.+)\.(.+)/) {
	$dir = $1;
	$extension = $3;
	$component = $2;
	$ergatis_file = '/opt/ergatis/docs/'.$component.'.'.$4;
}

unless( -e $clovr_file || -e $ergatis_file ) {
	die "Usage: program.pl clovr_component_file.default.config\n";
}

my $clovr = get_content( $clovr_file );
my $ergatis = get_content( $ergatis_file );
my $tmpl = get_content( "$dir/pipeline_tmpl.config.bak" );

check_for_var_diff( $clovr, $ergatis, "clovr", "ergatis" );
check_for_var_diff( $ergatis, $clovr, "ergatis", "clovr" );
print_ones_with_diff_vals( $clovr, $ergatis, $tmpl );

exit $?;


sub print_ones_with_diff_vals {
	my ($clovr, $ergatis, $tmpl) = @_;
	print STDOUT "[",$component," ", $extension,"]","\n";
	my $token = "[".$component." ".$extension."]";
	foreach my $key( keys %{$$tmpl{$token}} ) {
                print STDOUT $key,"=",$$tmpl{$token}{$key},"\n";
        }
	foreach my $key1( keys %$clovr ) {
		print STDERR $key1,"\n";
		foreach my $key2( keys %{$$clovr{$key1}} ) {
			unless( $$clovr{$key1}{$key2} eq $$ergatis{$key1}{$key2} || exists $$tmpl{$token}{$key2} ) {
				print STDERR "Clovr: ", $key2, "=", $$clovr{$key1}{$key2}, "\n"; 
				print STDOUT $key2,"=",$$clovr{$key1}{$key2},"\n";
			}
		}
	}
	print STDOUT "\n";
}

sub get_content {
	my ($file) = @_;
	open(FH,"<$file") or die "Error in opening the file, $file, $!\n";
	my $root = {};
	my $domain = undef;
	while(my $line = <FH>) {
		chomp $line;
		if( $line =~ /^(\[.+\])/ ) {
			$domain = $1;
		} elsif ( $line =~ /^(\$;.+\$;)\s*=\s*(.+)\s*$/ ) {
			my ($key,$value) = ($1,$2);
			if( $value =~ m/<<(.+)/ ) {
				my $token = $1;
				$value .= "\n";
				do {
					$line = <FH>;
					$value .= $line;
					chomp $line;
				} while( $line ne "$token" );
				chomp $value;
			}
			$$root{$domain}{$key} = $value;
		}		
	}
	close FH;
	return $root;
}

sub check_for_var_diff {
	my ($check, $check_against, $check_string, $check_against_string) = @_;
	print STDERR "Vars that are in $check_string but not in $check_against_string\n";
	print STDERR "----------------------------------------------------------------\n";
	foreach my $key1 ( keys %$check ) {
		print STDERR "$key1\n";
		foreach my $key2 ( keys %{$$check{key1}} ) {
			unless( exists $$check_against{ $key1 }{ $key2 } ) {
				print STDERR $key2,"=",$$check{$key1}{$key2},"\n";
			}
		}
		print STDERR "\n";
	}
}
