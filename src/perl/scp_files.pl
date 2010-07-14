#!/usr/bin/perl

=head1 NAME

scp_files.pl - Will scp some files.

=head1 SYNOPSIS

USAGE: scp_files.pl 
      --input_directory=/path/to/input_dir
      --output_directory=/path/to/output
    [ --file_extension_filters=btab,list
      --output_host=destination_server
      --tmp_dir=/tmp
      --log 
      --help
    ]

=head1 OPTIONS

B<--input_directory,-i>
    Directory where the files originate

B<--output_directory,-o>
    Directory where the files will end up

B<--file_extension_filters,-f>
    Script will filter files based on the presence of these extensions. Comma separated list.
    If not specified, will take all files

B<--output_host,-O>
    Server files will be copied to. If not specified, assumes local

B<--tmp,-t>
    Tmp directory. Files will be copied to this directory and a gzipped tarball will be created there.
    Default=/tmp

B<--log,-o> 
    Log file

B<--help,-h>
    This help message

=head1  CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut

############################## includes ##############################
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Find;
use Data::Dumper;
######################################################################

############################### global vars ##########################
my $destination_directory;
my $destination_host;
my $origin_directory;
my $tmp_dir = "/tmp";
my @file_extension_filters;
my $fe_regex;
my @files_to_copy;
######################################################################

######################### get options ################################
my %options = ();
my $results = GetOptions (\%options, 
                          'input_directory|i=s',
                          'output_directory|o=s',
                          'output_host|O=s',
                          'file_extension_filters|f=s',
                          'tmp_dir|t=s',
                          'help|h');

&check_options(\%options);
######################################################################

################################ main ################################

#find the files which need to be copied
find( \&find_files, $origin_directory );

#created the tmp directory
my $tmp_dirname = "$tmp_dir/${$}_".int(rand(2000));
eval {
    mkdir( $tmp_dirname );
    print "making temp directory: $tmp_dirname\n";
};
if( $@ ) {
    die("Could not create tmp directory $tmp_dirname ($@)");
}

#copy the files to the tmp directory
&copy_files_to_tmp_dir( $tmp_dirname, @files_to_copy );

#create the archive
my $tarbasename = "${$}_".int(rand(2000)).".tar.gz";
my $tar_file = "$tmp_dir/$tarbasename";
&create_tar( $tmp_dirname, $tar_file );

#create the output directory and set permissions
#so that no one can read/write to it
&create_directory( $destination_directory );

#copy the file over
my $cp_cmd = "cp $tar_file $destination_directory";
if( defined( $destination_host ) ) {
    $cp_cmd = "scp $tar_file $destination_host:$destination_directory";
}
run_cmd( $cp_cmd );

#untar the file
my $untar_cmd = "cd $destination_directory && tar -xzf $tarbasename";
if( defined( $destination_host ) ) {
    $untar_cmd = "ssh $destination_host '$untar_cmd'";
}
run_cmd( $untar_cmd );

#remove the remote copy of the tar file
my $rm_cmd = "rm $destination_directory/$tarbasename";
if( defined( $destination_host ) ) {
    $rm_cmd = "ssh $destination_host '$rm_cmd'";
}
run_cmd( "$rm_cmd" );

#open up permissions for the file
my $chmod_cmd = "chmod -R 774 $destination_directory";
if( defined( $destination_host ) ) {
    $chmod_cmd = "ssh $destination_host '$chmod_cmd'";
}
run_cmd( "$chmod_cmd" );

#Clean up the tmp directory
run_cmd( "rm -rf $tmp_dirname" );
run_cmd( "rm $tar_file" );
######################################################################

############################ subroutines #############################
sub create_directory {
    my ($destination_directory) = @_;
    
    #check to see if parent dirs exist
    my $parent_dir = $1 if( $destination_directory =~ m|(.*)/[^/]+$| );
    
    if( ! -e $parent_dir ) {
        my $cmd = "mkdir -p $parent_dir";
        if( defined( $destination_host ) ) {
            $cmd = "ssh $destination_host '$cmd'";
        }
        run_cmd( $cmd );
    }

    my $cmd = "mkdir -m 700 $destination_directory";
    if( defined( $destination_host ) ) {
        $cmd = "ssh $destination_host '$cmd'";
    }
    run_cmd( "$cmd" );

}
sub create_tar {
    my ($dir, $filename) = @_;
    my $cmd = "cd $dir && tar -czf $filename ./*";
    run_cmd( $cmd );
}

sub copy_files_to_tmp_dir {
    my ($tmp_dirname, @files_to_copy) = @_;

    map {
        my $filename = $_;
        $filename =~ s/$origin_directory\/?//;

        my $dest_subdir = $1 if( $filename =~ m|^(.*)/[^/]+$| );

        my $new_dir;
        if( defined( $dest_subdir ) ) {
            $new_dir = $tmp_dirname."/".$dest_subdir;
            
            if( !-d $new_dir ) {
                run_cmd( "mkdir -p $new_dir" );
            }

        }

        my $df = $tmp_dirname."/".$filename;
        my $cmd = "cp $_ $df";
        run_cmd($cmd);
    } @files_to_copy;
}

sub find_files {
    return if( $_ =~ /^\./ );
    if( !-d && /$fe_regex/ ) {
        push(@files_to_copy, $File::Find::name);
    }
}

sub run_cmd {
    my ($cmd) = @_;
    if( system( $cmd ) ) {
        die("Could not run cmd [$cmd]. Exited with error code $?");
    }
    return 1;
}

sub check_options {
    my ($opts) = @_;

    if( $opts->{'help'} ) {
        pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
        exit(0);
    }

    foreach my $req_opt( qw(input_directory output_directory) ) {
        die("Option $req_opt is required") unless( $opts->{$req_opt} );
    }

    $destination_directory = $opts->{'output_directory'};
    $origin_directory = $opts->{'input_directory'};
    
    $destination_host = $opts->{'output_host'} if( $opts->{'output_host'} );

    @file_extension_filters = split( /[,\s]+/, $opts->{'file_extension_filters'} )
        if( $opts->{'file_extension_filters'} );

    if( @file_extension_filters ) {
        map { $_ .= "\$"; } @file_extension_filters;
        my $t = join("|", @file_extension_filters);
        $fe_regex = qr/$t/;
    } else {
        $fe_regex = qr/.*/;
    }

    $tmp_dir = $opts->{'tmp_dir'} if( $opts->{'tmp_dir'} );

}

####################### EOF #########################################
