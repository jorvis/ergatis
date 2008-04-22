#!/usr/local/bin/perl

=head1 NAME

reset_subgroups.pl - intelligently reset uncomplete jobs in incomplete ergatis groups

=head1 SYNOPSIS

    USAGE: reset_subgroups.pl

=head1 OPTIONS

    B<--iter_xml,-i> : the path to an iter xml file (i1.iter.xml.gz)

=head1  DESCRIPTION

When given an ergatis component iterator's workflow xml file, this script will examine the 
state of the subgroups.  Any subgroups not marked 'complete' will have all jobs reset 
as incomplete. this 'incomplete' state will cascade up through to the iterater state as well.

=head1  INPUT

A single iterator xml file path.

=head1  OUTPUT

temp files are generated, then copied over as needed.

=head1  CONTACT

    Jason Inman
    jinman@jcvi.org

=cut

use warnings;
use strict;
use XML::Twig;
use Getopt::Long;
use Pod::Usage;
use File::Copy;

umask 0002;

my $iter_xml;
my %opts;

GetOptions(\%opts,  'iter_xml|i=s',
                    'help|h',
          ) || die "Error getting options?  Looks like it";
pod2usage( {-exitval => 0, -verbose => 1}) if ($opts{'help'});

&check_params;

# setup the iter.xml parsing
my $ifh = &open_xml($iter_xml);

my $iter_tmp = "$iter_xml.tmp";
my $iter_tmp_fh = &open_out_xml($iter_tmp);

my $altered = 0;

my $iter_twig = new XML::Twig (twig_roots => {
                                'commandSet[@type="remote-serial"]' =>
                                    sub { &iter_twig_handler(@_,\$altered) },
                                             },
                               twig_print_outside_roots => $iter_tmp_fh,
                               pretty_print => 'indented',
                              );

# Can't use parsefile because we might be using a binmode of :gzip here:
$iter_twig->parse($ifh);

close $ifh;
close $iter_tmp_fh;

if ($altered) {

    move($iter_xml, "$iter_xml.original") if -e $iter_xml;
    copy($iter_tmp, $iter_xml) || die "Problems copying new iter.xml over: $!\n";

}

exit(0);

sub iter_twig_handler {
# top level of parsing... look at each incomplete group.  clear dcespec and send
# the g###.xml.gz to be parsed if not marked complete

    my ($tree, $elem, $altered) = @_;

    if ($elem->first_child('state')->text ne 'complete') {

        # mark our flag that we're going to change stuff
        $$altered++;
        &clear_dceSpec( $elem->first_child('dceSpec') );

        # and now send the group xml to be parsed
        my $group_xml = $elem->first_child('fileName')->text;

        &parse_group($group_xml);

        # oh yeah, and mark it as incomplete
        $elem->first_child('state')->set_text('incomplete');

    }
    
    # print this element now that we're done with it
    $elem->print();

}

sub clear_dceSpec {
# given a reference to a twig element (better be a dceSpec, cuz we don't check)
# get rid of everything except 'OS' and 'group' elements

    my $elem_ref = shift;

    foreach my $trash ('jobID','executionHost','log','runtime','evictable',
                       'duration','priority','reqStartTime','gridID') {
        $elem_ref->first_child($trash)->delete;
    }

}

sub parse_group {
# parse a given group xml.  In other words, get the g#.iter.xml.gz path and send
# it to be parsed.

    my $group_xml = shift;

    my $g_ifh = &open_xml($group_xml);

    my $group_tmp = "$group_xml.tmp";
    my $g_tmp_fh = &open_out_xml($group_tmp);

    my $altered = 0;

    my $g_twig = new XML::Twig (twig_roots => {
                                    'commandSet[@type="serial"]' =>
                                        sub { &g_twig_handler(@_,\$altered) },
                                                 },
                                   twig_print_outside_roots => $g_tmp_fh,
                                   pretty_print => 'indented',
                                  );

    # Can't use parsefile because we might be using a binmode of :gzip here:
    $g_twig->parse($g_ifh);

    close $g_ifh;
    close $g_tmp_fh;

    if ($altered) {
        move($group_xml, "$group_xml.original") if -e $group_xml;
        move($group_tmp, $group_xml) || die "Problems copying new group.xml over: $!\n";
    } 

}

sub g_twig_handler {
# second level of parsing... need to mark the g##.iter.xml.gz command incomplete
# and then send it off to be parsed
    my ($tree, $elem, $altered) = @_;

    $elem->first_child('state')->set_text('incomplete');
    my $g_iter_xml = $elem->first_child('commandSet')->first_child('fileName')->text;

    &parse_g_iter($g_iter_xml);

    # print this element now that we're done with it
    $elem->print();

}

sub parse_g_iter {
# parse the group iter xml, like g#.iter.xml.gz 
# in other words, send any subgroups to be parsed that aren't marked 'complete'
# mark them 'incomplete', too.

    my $g_iter_xml = shift;

    my $g_iter_fh = &open_xml($g_iter_xml);

    my $g_iter_tmp = "$g_iter_xml.tmp";
    my $g_iter_tmp_fh = &open_out_xml($g_iter_tmp);

    my $altered = 0;

    my $g_iter_twig = new XML::Twig (twig_roots => {
                                    'commandSet[@type="parallel"]' =>
                                        sub { &g_iter_twig_handler(@_,\$altered) },
                                                 },
                                   twig_print_outside_roots => $g_iter_tmp_fh,
                                   pretty_print => 'indented',
                                  );

    # Can't use parsefile because we might be using a binmode of :gzip here:
    $g_iter_twig->parse($g_iter_fh);

    close $g_iter_fh;
    close $g_iter_tmp_fh;

    if ($altered) {
        move($g_iter_xml, "$g_iter_xml.original") if -e $g_iter_xml;
        move($g_iter_tmp, $g_iter_xml) || die "Problems copying new g#.iter.xml over: $!\n";
    } 
  
}

sub g_iter_twig_handler {
# the third level of parsing... check each subgroup and for anything other than
# 'complete' state, mark 'incomplete' and send for parsing

    my ($tree, $elem, $altered) = @_;

    my $child = $elem->first_child('state');
    $child->set_text('incomplete');

    # loop through, sending any non 'complete' commandSet fileNames to be parsed
    while ($child = $child->next_sibling('commandSet') ) {
        
        if ($child->first_child('state')->text ne 'complete') {

            $child->first_child('state')->set_text('incomplete');

            my $filename = $child->first_child('fileName')->text;
            &parse_subgroup($filename);

        }

    }

    # and write out the whole element
    $elem->print();

}

sub parse_subgroup {
# For the given subgroup_xml, check the job states and 
# for any not marked complete, mark incomplete 
# again, done by twig handling :)

    my $sg_xml = shift;

    my $sg_fh = &open_xml($sg_xml);

    my $altered = 0;

    # set up a temp file
    my $sg_tmp = "$sg_xml.tmp";
    my $sg_tmp_fh = &open_out_xml($sg_tmp);

    my $sg_twig = new XML::Twig (twig_roots => {
                                    'commandSet' => sub{ &sg_twig_handler(@_,\$altered) }
                                               },
                                 twig_print_outside_roots => $sg_tmp_fh,
                                 pretty_print => 'indented',
                                );

    $sg_twig->parse($sg_fh);
    close $sg_fh;
    close $sg_tmp_fh;

    if ($altered) {
         move($sg_xml, "$sg_xml.original") if -e $sg_xml;
         move($sg_tmp, $sg_xml) || die "Problem moving subgroup.xml $!";
    }

}

sub sg_twig_handler {
# the lowest level of parsing...
# check an individual subgroup's commands.
# if anything is changed, increment the counter

    my ($tree, $elem, $altered) = @_;

    # somehow, this next line passes over subgroup xmls not even started yet.
    if (defined $elem->first_child('state')) { 

        # If the whole commandSet is not marked complete...
        if ($elem->first_child('state')->text ne 'complete') {

            # increment out counter so we know we have to do something...
            $$altered++;

            # loop through each command and set it's state to 'incomplete'
            # starting with the first one....
            my $child_elem = $elem->first_child('command');
            $child_elem->first_child('state')->set_text('incomplete'); 

            # and looping through the rest:        
            while ( $child_elem = $child_elem->next_sibling('command') ) {
                $child_elem->first_child('state')->set_text('incomplete');
            }
        }
    }
    # now print the whole commandSet block to the output
    $elem->print();
    
}

sub open_out_xml {
# really, open any file for closing.  This function returns a reference to a filehandle
# for writing in either :gzip binmode or normal depending on the extension

    my $file_target = shift;

    my $fh;
    if ($file_target =~ /\.gz.tmp$/) {
        open ( $fh, ">:gzip", $file_target) || die "Couldn't open to write to compressed $file_target: $!";    
    } else {
        open ( $fh, ">$file_target") || die "Couldn't open to write to  $file_target: $!";
    }

    return ($fh);

}

sub open_xml {
# really, open any file.  This function returns a reference to a filehandle
# for reading in either :gzip binmode or normal depending on the extension

    my $file_target = shift;

    my $fh;
    if ($file_target =~ /\.gz$/) {
        open ( $fh, "<:gzip", $file_target) || die "Couldn't open compressed $file_target: $!";    
    } else {
        open ( $fh, "<$file_target") || die "Couldn't open $file_target: $!";
    }

    return ($fh);

}

sub check_params {
# check to make sure options are ok.

    my $errors = '';

    if (defined $opts{'iter_xml'}) {
        $iter_xml = $opts{'iter_xml'};
    } else {
        my $errors .= "Please supply a --iter_xml\n";
    }

    die $errors if $errors;

    return(0);

}
