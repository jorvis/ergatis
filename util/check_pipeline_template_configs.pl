#!/usr/bin/env perl

=head1 NAME

check_pipeline_template_configs.pl - Will check the variables of a pipeline template's config files against a specific ergatis install.  Reports any differences in the variables between each.

=head1 SYNOPSIS

USAGE: create_evidence_file_mapping.pl
    --pipeline_layout=/path/to/pipeline.layout
    --ergatis_install=/path/to/ergatis_install/
  [ --help ]

=head1 OPTIONS

B<--pipeline_layout,-p>

    REQUIRED. Path to pipeline.layout file.

B<--ergatis_install,-e>

    REQUIRED Path to ergatis install directory to check pipeline template against

B<--help,-h>

    Print this message

=head1  DESCRIPTION

    Since the config files for pipeline templates are versioned separately, it is likely that
    the config files will be out of date or out of sync with a various version.  This script will 
    check for these inconsistencies.  
 
=head1  INPUT
    
    Just the paths to the pipeline layout file and the ergatis install to check it against is required.

=head1 OUTPUT

    Prints to standard out.  Will print missing variables from either the pipeline template config or
    the ergatis install config.  Will print a confirmation if the config file matches.


=head1  CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::OpenFile qw(open_file);
use File::Basename;
use Config::IniFiles;
use Data::Dumper;
use Pod::Usage;

my $pipeline_layout;
my $ergatis_install;
&get_options;

#Grab all the componet names
my @components;
my $in = open_file( $pipeline_layout, 'in' );
while(<$in>) {
    if(/\<name\>([^\<]+)\</ ) {
        next if( $1 eq 'start pipeline:' || $1 eq 'start' );
        push( @components, $1 );
    }
}
close($in);

#Locate the configs
my ($name, $path, $suffix) = fileparse( $pipeline_layout );

foreach my $component ( @components ) {
    my $comp_config = $path."/".$component.".config";
    die("Could not find config file for component name $component. Expected $comp_config")
        unless( -e $comp_config );

    my $component_name = $1 if( $component =~ /^([^\.]+)\./ );

    my $template_config = $ergatis_install."/docs/$component_name.config";
    die("Could not find template config file for componetn $component_name. Expected $template_config")
        unless( -e $template_config );

    &compare_config_vars( $comp_config, $template_config );
}

sub compare_config_vars {
    my ($comp, $temp) = @_;
    my $component_config = new Config::IniFiles( -file => $comp );
    my %cc_params;

    #Gather all the vars
    my @sections = $component_config->Sections();

    foreach my $section ( @sections ) {
        my @params = $component_config->Parameters( $section );
        my $tmp = {};
        map{ $tmp->{$_} = 1 } @params;
        $cc_params{$section} = $tmp;
    }

    undef $component_config;

    #Gather params from template config and check
    my $template_config = new Config::IniFiles( -file => $temp );
    my %tc_params;

    my @t_sections = $template_config->Sections();
    
    foreach my $t_section ( @t_sections ) {

        my @t_params = $template_config->Parameters($t_section);
        my $tc_param_hash = {};
        map { $tc_param_hash->{$_} = 1 } @t_params;
        $tc_params{ $t_section } = $tc_param_hash;

        if( exists( $cc_params{ $t_section } ) ) {
            my $config_params = $cc_params{$t_section};
            
            foreach my $t_param( @t_params ) {
                
                if( exists( $config_params->{$t_param} ) ) {
                    delete $config_params->{ $t_param };
                    delete $tc_param_hash->{ $t_param };
                }

            }
        }
    }

    &print_differences( $comp, \%cc_params, $temp, \%tc_params );

}

sub print_differences {
    my ($comp, $cc_params, $temp, $tc_params) = @_;

    my @c_diffs;
    map { push( @c_diffs, keys %{$cc_params->{$_}} ) } keys %$cc_params;

    my @t_diffs;
    map { push( @t_diffs, keys %{$tc_params->{$_}} ) } keys %$tc_params;

    if( @t_diffs == 0 and @c_diffs == 0 ) {
        #print "$comp OK\n";
    } else {
        print basename($comp)." not OK\n";
        if( @t_diffs > 0 ) {
            print "\tThere are ".scalar(@t_diffs)." variables in template config that are not found in current config\n\t\t";
            my $p_string = join( "\n\t\t", @t_diffs );
            print "$p_string\n";

        }
        if( @c_diffs > 0 ) {
            print "\tThere are ".scalar(@c_diffs)." variables in current config that are not found in template config\n\t\t";
            my $p_string = join( "\n\t\t", @c_diffs );
            print "$p_string\n";
        }
    }
}

sub get_options {
    my %options = ();
    my $results = GetOptions (\%options, 
                              'pipeline_layout|p=s',
                              'ergatis_install|e=s',
                              'help|h');

    if( $options{'help'} ) {
        &_pod;
    }

    if( $options{'pipeline_layout'} ) {
        $pipeline_layout = $options{'pipeline_layout'};
        die("File does not exist: $pipeline_layout")
            unless( -e $pipeline_layout && -f $pipeline_layout );
    } else {
        &_pod("Option pipeline_layout is required.", 1);
    }

    if( $options{'ergatis_install'} ) {
        $ergatis_install = $options{'ergatis_install'};
    } else {
        &_pod("Option ergatis_install is required.",1);
    }

    
}

sub _pod {
    my ($message, $exitval) = @_;
    $message = " " unless( $message );
    $exitval = 0 unless( $exitval );
    &pod2usage( { 
        -message => $message,
        -exitval => 0,
        -verbose => 5, } );
    
}
