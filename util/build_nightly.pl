#!/usr/bin/perl

use strict;
use Config::IniFiles;
use File::Path;

## this use lib needed for SVN::Agent
use lib '/home/jorvis/lib/lib/perl5/site_perl/5.8.8';
use SVN::Agent;

umask(0022);

## 
#print STDERR "This script has been temporarily disabled to assist Emmanuel complete a critical pipeline.\n";
#print STDERR "To reinstate, please contact Joshua. (2012-07-30)\n\n";
#exit(1);

##########################################
my $install_base='/usr/local/projects/ergatis/package-nightly/';
my $ergatis_git_export_path='https://github.com/jorvis/ergatis.git';	# This uses Git currently.  Others use SVN.
my $bsml_svn_export_path='http://svn.code.sf.net/p/bsml/code/release';
my $coati_svn_export_path='http://svn.code.sf.net/p/coati-api/code/release/coati_install';
my $shared_prism_svn_export_path='http://svn.code.sf.net/p/prism-api/code/release/shared_prism';
my $chado_prism_svn_export_path='http://svn.code.sf.net/p/prism-api/code/release/chado_prism';
my $prok_prism_svn_export_path='http://svn.code.sf.net/p/prism-api/code/release/prok_prism';
my $euk_prism_svn_export_path='http://svn.code.sf.net/p/prism-api/code/release/euk_prism';
my $chado_schema_svn_export_path='http://svn.code.sf.net/p/prism-api/code/release/chado';

## this directory will be created.  If it exists already, it and everything under it
#   will be removed
my $tmp_area = '/tmp/build_nightly';

## if this is set this script will keep any settings from this config file, overwriting
#   the defaults from the SVN one.
my $software_config_kept = '/usr/local/projects/ergatis/software.config';


##########################################

my $software_config;
if ( $software_config_kept ) {
    $software_config = new Config::IniFiles( -file => $software_config_kept ) ||
                        die "failed to open old software config: $!";
}

clear_export_area($tmp_area);
clear_install_area($install_base);

install_ergatis($install_base);
install_bsml($install_base);

install_coati($install_base);
install_shared_prism($install_base);
install_chado_prism($install_base);
install_prok_prism($install_base);
install_euk_prism($install_base);
install_chado_schema($install_base);

replace_software_config_values($software_config);
set_idgen_configuration($install_base);

## for now, also recursively copy the IGS lib dir into the area
#`cp -r /usr/local/projects/ergatis/package-latest/lib/perl5/IGS $install_base/lib/perl5`;

sub clear_export_area {
    my $base = shift;
    
    rmtree($tmp_area);
    mkdir($tmp_area) || die "failed to create temp directory: $tmp_area";
}

sub clear_install_area {
    my $base = shift;

    for my $subdir ( qw( bin docs lib man samples ) ) {
        rmtree("$base/$subdir");
    }
    
    for my $file ( qw( software.config ) ) {
        if ( -e "$base/$file" ) {
            unlink( "$base/$file" );
        }
    }
};


sub install_bsml {
    my $base = shift;
    
    my $tmp_build_area = "$tmp_area/bsml";

    print STDERR "exporting out BSML\n";    
    #my $svn = SVN::Agent->new( {path => "$tmp_build_area"} );
    #$svn->checkout( $bsml_svn_export_path );
    my $cmd = "echo p | svn export $bsml_svn_export_path $tmp_build_area";
    system($cmd);
    print STDERR "installing BSML\n";
    chdir("$tmp_build_area/bsml-vNrNbN/") || die "couldn't cd into $tmp_build_area/bsml-vNrNbN";
    run_command( "perl -I /home/jorvis/lib/perl5/5.8.8/ Makefile.PL INSTALL_BASE=$base SCHEMA_DOCS_DIR=$base/docs" );
    run_command( "make" );
    run_command( "make install" );

}

sub install_chado_schema {
    my $base = shift;
 
    my $tmp_build_area = "$tmp_area/chado_schema";

    print STDERR "exporting out Chado schema\n";    
    #my $svn = SVN::Agent->new( {path => "$tmp_build_area"} );
    #$svn->checkout( $chado_schema_svn_export_path );
    
    my $cmd = "echo p | svn export $chado_schema_svn_export_path $tmp_build_area";
    system($cmd);
    
    print STDERR "installing Chado schema\n";
    chdir("$tmp_build_area/chado-vNrNbN/") || die "couldn't cd into $tmp_build_area/chado-vNrNbN";
    run_command( "perl install.pl INSTALL_BASE=$base" );
}

sub install_chado_prism {
    my $base = shift;
 
    my $tmp_build_area = "$tmp_area/chado_prism";

    print STDERR "exporting out Prism (chado)\n";    
    #my $svn = SVN::Agent->new( {path => "$tmp_build_area"} );
    #$svn->checkout( $chado_prism_svn_export_path );
    
    my $cmd = "echo p | svn export $chado_prism_svn_export_path $tmp_build_area";
    system($cmd);

    print STDERR "installing Prism (chado)\n";
    chdir("$tmp_build_area/chado_prism-vNrNbN/") || die "couldn't cd into $tmp_build_area/chado_prism-vNrNbN";
    run_command( "perl -I /home/jorvis/lib/perl5/5.8.8/ Makefile.PL INSTALL_BASE=$base" );
    run_command( "make" );
    run_command( "make install" );
}

sub install_coati {
    my $base = shift;
    
    my $tmp_build_area = "$tmp_area/coati";

    print STDERR "exporting out Coati\n"; print "$coati_svn_export_path\n"; 
    #my $svn = SVN::Agent->new( {path => "$tmp_build_area"} );
    #$svn->checkout( $coati_svn_export_path );
    
    my $cmd = "echo p | svn export $coati_svn_export_path $tmp_build_area";
    system($cmd);

    print STDERR "installing Coati\n";
    chdir("$tmp_build_area/coati_install-vNrNbN/") || die "couldn't cd into $tmp_build_area/coati_install-vNrNbN";
    run_command( "perl -I /home/jorvis/lib/perl5/5.8.8/ Makefile.PL INSTALL_BASE=$base" );
    run_command( "make" );
    run_command( "make install" );

}

sub install_ergatis {
    my $base = shift;
    
    my $tmp_build_area = "$tmp_area/ergatis";

    print STDERR "exporting out Ergatis\n";    
    #my $svn = SVN::Agent->new( {path => "$tmp_build_area"} );
    #$svn->checkout( $ergatis_svn_export_path );
    
    my $cmd = "echo p | git clone $ergatis_git_export_path $tmp_build_area/ergatis-trunk";
    system($cmd);

    print STDERR "installing Ergatis\n";
    chdir("$tmp_build_area/ergatis-trunk") || die "couldn't cd into $tmp_build_area/ergatis-trunk";
    
    # If we are pulling Ergatis from github we are going to need to move some directories around in order
    # to replace what svn:externals did for us
    # TODO: Figure out a better way to achieve this without the special case just for Ergatis 
    run_command("mv install/* .; mv src/R .; mv src/python .; mv templates/pipelines global_pipeline_templates; mv src/perl bin; mv src/shell .; mv src/c .; rm -rf src; mv c src; cp -R lib/Ergatis htdocs/cgi/");
    
    run_command( "perl -I /home/jorvis/lib/perl5/5.8.8/ Makefile.PL INSTALL_BASE=$base" );
    run_command( "make" );
    run_command( "make install" );
}

sub install_euk_prism {
    my $base = shift;
 
    my $tmp_build_area = "$tmp_area/euk_prism";

    print STDERR "exporting out Prism (euk)\n";    
    #my $svn = SVN::Agent->new( {path => "$tmp_build_area"} );
    #$svn->checkout( $euk_prism_svn_export_path );
    
    my $cmd = "echo p | svn export $euk_prism_svn_export_path $tmp_build_area";
    system($cmd);

    print STDERR "installing Prism (euk)\n";
    chdir("$tmp_build_area/euk_prism-vNrNbN/") || die "couldn't cd into $tmp_build_area/euk_prism-vNrNbN";
    run_command( "perl -I /home/jorvis/lib/perl5/5.8.8/ Makefile.PL INSTALL_BASE=$base" );
    run_command( "make" );
    run_command( "make install" );
}

sub install_prok_prism {
    my $base = shift;
 
    my $tmp_build_area = "$tmp_area/prok_prism";

    print STDERR "exporting out Prism (prok)\n";    
    #my $svn = SVN::Agent->new( {path => "$tmp_build_area"} );
    #$svn->checkout( $prok_prism_svn_export_path );
    
    my $cmd = "echo p | svn export $prok_prism_svn_export_path $tmp_build_area";
    system($cmd);

    print STDERR "installing Prism (prok)\n";
    chdir("$tmp_build_area/prok_prism-vNrNbN/") || die "couldn't cd into $tmp_build_area/prok_prism-vNrNbN";
    run_command( "perl -I /home/jorvis/lib/perl5/5.8.8/ Makefile.PL INSTALL_BASE=$base" );
    run_command( "make" );
    run_command( "make install" );
}

sub install_shared_prism {
    my $base = shift;
 
    my $tmp_build_area = "$tmp_area/shared_prism";

    print STDERR "exporting out Prism (shared)\n"; print "$shared_prism_svn_export_path\n";    
    #my $svn = SVN::Agent->new( {path => "$tmp_build_area"} );
    #$svn->checkout( $shared_prism_svn_export_path );
    
    my $cmd = "echo p | svn export $shared_prism_svn_export_path $tmp_build_area";
    system($cmd);

    print STDERR "installing Prism (shared)\n";
    chdir("$tmp_build_area/shared_prism-vNrNbN/") || die "couldn't cd into $tmp_build_area/shared_prism-vNrNbN";
    run_command( "perl -I /home/jorvis/lib/perl5/5.8.8/ Makefile.PL INSTALL_BASE=$base" );
    run_command( "make" );
    run_command( "make install" );
}

sub replace_software_config_values {
    my $old_config = shift;
    
    my $new_config_path = "$install_base/software.config";
    
    my $new_config = new Config::IniFiles( -file => $new_config_path ); 
#                        die "failed to open old software config: $!";

    my $reset_permission_to = undef;

    if ( ! -w $new_config_path ) {
        $reset_permission_to = (stat $new_config_path)[2] & 07777;
        chmod(0644, $new_config_path) || die "couldn't chmod $new_config_path for writing";
    }  
                              
    for my $section ( $new_config->Sections() ) {
        if ( $old_config->SectionExists($section) ) {
            
            for my $parameter ( $new_config->Parameters($section) ) {
                if ( defined $old_config->val($section, $parameter) ) {
                    $new_config->setval( $section, $parameter, $old_config->val($section, $parameter) );
                }
            }
        }
    }
    
    $new_config->RewriteConfig;
    
    if ( defined $reset_permission_to ) {
        chmod($reset_permission_to, $new_config_path) || die "couldn't restore permissions on $new_config_path (tried to set to $reset_permission_to)";
    } 
};

sub run_command {
    my $cmd = shift;
    
    system($cmd);
    
    if ( $? == -1 ) {
        die "failed to execute command ($cmd): $!\n";
    } elsif ( $? & 127 ) {
         my $out = sprintf "command ($cmd): child died with signal %d, %s coredump\n",
                    ($? & 127),  ($? & 128) ? 'with' : 'without';
         die($out);
    }
}

sub set_idgen_configuration {
    my $base = shift;
    
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
                $line =~ s/\:\:DefaultIdGenerator/\:\:IGSIdGenerator/g;
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






