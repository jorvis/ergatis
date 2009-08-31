#!/usr/bin/perl

=head1 NAME

build_ergatis.pl - automatically install ergatis and its required packages

=head1 SYNOPSIS

USAGE: build_ergatis.pl 
            --install_base=/path/to/install/dir
          [ --htdocs_area=/path/to/www/ergatis
            --tmp_area=/path/to/tmp/dir
            --software_config=/path/to/config.file
            --id_generator={default, igs}
	    --lib=/path/to/perl/lib/dir
            --log=/path/to/file
            --help
          ]

=head1 OPTIONS

B<--install_base,-i>
    Location to install ergatis

B<--htdocs_area>
    optional.  Web files will be copied to this directory.  It should be accessible by your webserver.

B<--tmp_area,-t>
    optional.  Location to store temporary files.  Default is /tmp/build_ergatis.

B<--software_config,-s>
    optional.  If provided this script will keep any settings from this config file, overwriting the defaults from the SVN one.

B<--id_generator>
    optional.  The type of id generator perl module to use.  Will result in edits to Ergatis/IdGenerator/Config.pm & Coati/IdGenerator/Config.pm

B<--lib>
    optional. Perl library directory to be included in perl -I $opt_lib Makefile.pl

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1 CAVEATS

There's a good chance the script will not run automatically on your system the first time you try to use it.  This is because SVN requires
you to manually accept the server certificates from the SVN repositories.  If you run
 $ svn co 

=cut

use strict;
use warnings;
use Config::IniFiles;
use File::Path;

## this use lib needed for SVN::Agent
#use lib '/home/jorvis/lib/lib/perl5/site_perl/5.8.8';
# or run script w/ perl -I <above_dir>
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
#use SVN::Agent;

umask(0022);

# Valid ID generators the script knows about
# there are other modules, but they would need to be added to the script
my %ID_GENERATORS = (
		     igs => "IGSIdGenerator",
		     default => "DefaultIdGenerator"
		    );

my %opts = &parse_options();

# if an optional lib was provided, include it
my $opt_lib = '';
if ($opts{lib}) {
    $opt_lib = "-I $opts{lib}";
}

## open the log if requested
my $logfh;
if (defined $opts{log}) {
    open($logfh, ">$opts{log}") || die "can't create log file: $!";
}

##########################################
#my $install_base='/usr/local/projects/ergatis/package-nightly';
my $install_base=$opts{install_base};
# using http instead of https
my %SVN_CO_PATH = (
		   ergatis => 'http://ergatis.svn.sourceforge.net/svnroot/ergatis/release/trunk',
		   bsml => 'http://bsml.svn.sourceforge.net/svnroot/bsml/release',
		   coati_install => 'http://coati-api.svn.sourceforge.net/svnroot/coati-api/release/coati_install',
		   shared_prism => 'http://prism-api.svn.sourceforge.net/svnroot/prism-api/release/shared_prism',
		   chado_prism => 'http://prism-api.svn.sourceforge.net/svnroot/prism-api/release/chado_prism',
		   prok_prism => 'http://prism-api.svn.sourceforge.net/svnroot/prism-api/release/prok_prism',
		   euk_prism => 'http://prism-api.svn.sourceforge.net/svnroot/prism-api/release/euk_prism',
		   chado_schema => 'http://prism-api.svn.sourceforge.net/svnroot/prism-api/release/chado',
#		   igs_idgenerator => 'http://vader.igs.umaryland.edu/svnroot/ENGR/IGS-UIDGenerator/trunk/perl/lib/IGS',
		   );

## this directory will be created.  If it exists already, it and everything under it
#   will be removed
my $tmp_area = $opts{tmp_area};

## if this is set this script will keep any settings from this config file, overwriting
#   the defaults from the SVN one.
#my $software_config_kept = '/usr/local/projects/ergatis/package-latest/software.config';
my $software_config_kept = $opts{software_config};

##########################################

my $software_config = undef;
if ( $software_config_kept ) {
    $software_config = new Config::IniFiles( -file => $software_config_kept ) ||
                        die "failed to open old software config: $!";
}

clear_co_area($tmp_area);
clear_install_area($install_base);

install_ergatis($install_base);
(-e "$opts{install_base}/software.config") || die "Missing software config";

if ( $opts{htdocs_area} ) {
    install_ergatis_htdocs( $opts{htdocs_area} );
}

install_bsml($install_base);

install_coati($install_base);
install_shared_prism($install_base);
install_chado_prism($install_base);
install_prok_prism($install_base);
install_euk_prism($install_base);
install_chado_schema($install_base);
#HACK: unable to connect to vader svn from external sources
# once it's made public, we can fix this
#install_igs_idgenerator($install_base);
## for now, also recursively copy the IGS lib dir into the area
if ( $opts{id_generator} eq 'igs' ) {
    run_command("cp -r /usr/local/projects/ergatis/package-latest/lib/perl5/IGS $install_base/lib/perl5");
}
replace_software_config_values($software_config);

set_idgen_configuration($install_base);

_log("And we're done");

##########################################
## SUBROUTINES
##########################################

sub parse_options {
    my %options = ();
    my $results = GetOptions (\%options,
                              'install_base|i=s',
                              'tmp_area|t=s',
                              'htdocs_area=s',
			      'software_config|s=s',
			      'id_generator=s',
			      'lib=s',
                              'log|l=s',
                              'help|h') || pod2usage();

    ## display documentation
    if( $options{'help'} ){
        pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
    }

    ## make sure everything passed was peachy
    &check_parameters(\%options);

    return %options;
}

sub check_parameters {
    my $options = shift;

    ## make sure required arguments were passed
    my @required = qw( install_base );
    for my $option ( @required ) {
        unless  ( defined $options->{$option} ) {
            die "--$option is a required option";
        }
    }

    ##
    ## you can do other things here, such as checking that files exist, etc.
    ##
    if ($options->{software_config}) {
	unless ( -e $options->{software_config} ) {
	    die "Unable to read --software_config $options->{software_config}"
	}
    }

    ## handle some defaults
    $options->{tmp_area} = '/tmp/build_ergatis'  unless ($options->{tmp_area});

    $options->{id_generator} = 'default'  unless ($options->{id_generator});

    unless ( defined( $ID_GENERATORS{$options->{id_generator}} )) {
	die "Unknown --id_generator (".$options->{id_generator}.").  Must be one of ".join(', ', keys %ID_GENERATORS);
    }

}

sub _log {
    my $msg = shift;

    if ($logfh) {
        print {$logfh} "$msg\n";
    } else {
        print "LOG: $msg\n";
    }
}

sub clear_co_area {
    my $base = shift;
    
    _log("Removing and recreating tmp_area $tmp_area");
    rmtree($tmp_area);
    mkdir($tmp_area) || die "failed to create temp directory: $tmp_area";
}

sub clear_install_area {
    my $base = shift;

    _log("Clearing base $base");
    mkdir($base) unless (-d $base);
    for my $subdir ( qw( bin docs lib man samples ) ) {
        rmtree("$base/$subdir");
    }
    
    for my $file ( qw( software.config ) ) {
        if ( -e "$base/$file" ) {
            unlink( "$base/$file" );
        }
    }
};

sub install_ergatis {
    my $base = shift;    
    install_perl( $base, 'ergatis', 'trunk');
}

sub install_ergatis_htdocs {
    my $target = shift;
    my $source = "$tmp_area/ergatis/ergatis-trunk/htdocs";

    _log( "Moving htdocs to $source" );

    mkdir( $target ) unless (-e $target);

    run_command( "cp -r $source/* $target" );
}

sub install_bsml {
    my $base = shift;
    install_perl( $base, 'bsml', 'vNrNbN');
}

sub install_coati {
    my $base = shift;
    install_perl( $base, 'coati_install', 'vNrNbN');
}

sub install_shared_prism {
    my $base = shift;
    install_perl( $base, 'shared_prism', 'vNrNbN');
}

sub install_chado_prism {
    my $base = shift;
    install_perl( $base, 'chado_prism', 'vNrNbN');
}

sub install_prok_prism {
    my $base = shift;
    install_perl( $base, 'prok_prism', 'vNrNbN');
}

sub install_euk_prism {
    my $base = shift;
    install_perl( $base, 'euk_prism', 'vNrNbN');
}

sub install_chado_schema {
    my $base = shift;
 
    my $tmp_build_area = "$tmp_area/chado_schema";

    _log( "checking out Chado schema");    
    svn_checkout( $tmp_build_area, $SVN_CO_PATH{chado_schema} );
    
    _log( "installing Chado schema");
    chdir("$tmp_build_area/chado-vNrNbN/") || die "couldn't cd into $tmp_build_area/chado-vNrNbN";
    run_command( "perl install.pl INSTALL_BASE=$base" );
}


sub install_igs_idgenerator {
    my $base = shift;
    
    my $tmp_build_area = "$tmp_area/igs_idgenerator";
    
    _log( "checking out igs_idgenerator");
    svn_checkout( $tmp_build_area, $SVN_CO_PATH{igs_idgenerator} );

    my $cmd = "cp -r $tmp_area/igs_idgenerator/ENGR/IGS-UIDGenerator/trunk/perl/lib/IGS $install_base/lib/perl5";
    run_command( $cmd );

#    chdir("$tmp_build_area//") || die "couldn't cd into $tmp_build_area/chado-vNrNbN";
#    run_command( "perl install.pl INSTALL_BASE=$base" );    
}


sub replace_software_config_values {
    my $old_config = shift;
    
    my $new_config_path = "$install_base/software.config";
    
    my $new_config = new Config::IniFiles( -file => $new_config_path ) ||
                        die "failed to open new software config: $!";

    my $reset_permission_to = undef;

    if ( ! -w $new_config_path ) {
        $reset_permission_to = (stat $new_config_path)[2] & 07777;
        chmod(0644, $new_config_path) || die "couldn't chmod $new_config_path for writing";
    }  

    if (defined $old_config) {
	for my $section ( $new_config->Sections() ) {
	    if ( $old_config->SectionExists($section) ) {		
		for my $parameter ( $new_config->Parameters($section) ) {
		    if ( defined $old_config->val($section, $parameter) ) {
			$new_config->setval( $section, $parameter, $old_config->val($section, $parameter) );
		    }
		}
	    }
	}
    }
    
    $new_config->RewriteConfig;
    
    if ( defined $reset_permission_to ) {
        chmod($reset_permission_to, $new_config_path) || die "couldn't restore permissions on $new_config_path (tried to set to $reset_permission_to)";
    } 
};

sub install_perl {
    my ( $base, $package, $version ) = @_;
    
    my $tmp_build_area = "$tmp_area/$package";

    _log( "checking out $package-$version");
    svn_checkout( $tmp_build_area, $SVN_CO_PATH{$package} );

    _log( "installing $package-$version");
    chdir("$tmp_build_area/$package-$version/") || die "couldn't cd into $tmp_build_area/$package-$version";
    run_command( "perl $opt_lib Makefile.PL INSTALL_BASE=$base" );
    run_command( "make" );
    run_command( "make install" );
}


# ripped from SVN::Agent
sub svn_checkout {
    my ($path, $repository) = @_;
    mkdir($path) or die "Unable to create svn co " . $path.": $!";
#    my $svn_cmd = "svn co --trust-server-cert --non-interactive $repository $path";  # reqs SVN >= 1.6
    my $svn_cmd = "svn co --non-interactive $repository $path";
    run_command( $svn_cmd );
}

sub run_command {
    my $cmd = shift;
    
    _log("System Command: $cmd");

    system($cmd);
    
    if ( $? == -1 ) {
        die "failed to execute command ($cmd): $!\n";
    } elsif ( $? & 127 ) {
	my $out = sprintf "command ($cmd): child died with signal %d, %s coredump\n",
	    ($? & 127),  ($? & 128) ? 'with' : 'without';
	die($out);
    } elsif ( $? == 256 ){
	die "Error $?: you are probably missing perl module dependencies";
    } elsif ( $? != 0 ){
	die "Error $?: something probably went wrong";
    } else {
	_log( "Command ($cmd) exited with code $?" );
    }
}

sub set_idgen_configuration {
    my $base = shift;

    return if ( $opts{id_generator} eq 'default' );
    
    for my $conf_file ( "$base/lib/perl5/Ergatis/IdGenerator/Config.pm", 
                        "$base/lib/perl5/Coati/IdGenerator/Config.pm" ) {
    
        my $reset_permission_to = undef;
    
        if ( ! -w $conf_file ) {
            $reset_permission_to = (stat $conf_file)[2] & 07777;
            chmod(0644, $conf_file) || die "couldn't chmod $conf_file for writing";
        }
    
        my @lines = ();
    
        ## change $base/lib/perl5/Ergatis/IdGenerator/Config.pm
        open(my $cfh, $conf_file) || die "couldn't open $conf_file for reading";

        while ( my $line = <$cfh> ) {
            if ( $line !~ /^\#/ ) {
#                $line =~ s/\:\:DefaultIdGenerator/\:\:IGSIdGenerator/g;
                $line =~ s/\:\:DefaultIdGenerator/\:\:$opts{id_generator}/g;
            }
            
            push @lines, $line;
        }

        close $cfh;
        
        ## now stomp it with our current content
        open(my $ofh, ">$conf_file") || die "couldn't open $conf_file for writing\n";
            print $ofh join('', @lines);
        close $ofh;
        
        if ( defined $reset_permission_to ) {
            chmod($reset_permission_to, $conf_file) || die "couldn't restore permissions on $conf_file (tried to set to $reset_permission_to)";
        }
    }
}






