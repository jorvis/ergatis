#!/usr/bin/env perl

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use File::Basename;
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options,
              'username=s',
              'password=s',
              'db_name=s',
              'db_url=s',
              'pg_data_dir=s',
              'install_postgres',
              'start_postgres',
              'install_sybil',
              'create_site_config',
              'create_db',
              'cache_data',
              'install_mongo',
              'start_mongo',
              'config_url=s',
              'root_dir=s',
              'config_template=s',
              'sitename=s',
              'server=s',
              'conf_path=s',
              'cache_dir=s',
              'schema=s',
              'image_dir=s',
              'image_url=s',
              'java_path=s',
              'clustalw_path=s',
              'output_dir=s',
              'create_archive',
              'load_archive=s',
              'archive_is_tag',
              'help|h') || pod2usage();


my $PG_VERSION= '8.4';
my $ROOT = $options{root_dir} ? $options{root_dir} : '/mnt';
my $PG_DATA_DIR = $options{pg_data_dir} ? $options{pg_data_dir} : "$ROOT/pg_data";
my $SYBIL_TARBALL = 'http://cb2.igs.umaryland.edu/sybil_v1r6b1.tgz';
my $PATH_TO_SYBIL = '/var/www/sybil/current/';
my $CONF_PATH = $options{conf_path} ? $options{conf_path} : "$PATH_TO_SYBIL/cgi/conf/";
my $CONF_TEMPLATE = $options{config_template} ? $options{config_template} : "$CONF_PATH/master_template.conf";
my $JAVA_PATH = $options{java_path} ? $options{java_path} : `which java`;
my $CLUSTALW_PATH = $options{clustalw_path} ? $options{clustalw_path} : `which clustalw`;
my $SERVER = $options{server} ? $options{server} : 'localhost';
my $SITENAME = $options{sitename};

# Install Postgres
if(! -e "$PG_DATA_DIR/PG_VERSION" && $options{install_postgres}) {
    &run_cmd("apt-get update");
    &run_cmd("dpkg --configure -a"); # Need this in case something goes wrong above?
    &install_postgres();
    # Create user/database
#    my $cmd = "sudo -u postgres psql -c \"CREATE USER sybilcreate WITH PASSWORD 'createsybils' SUPERUSER\";sudo -u postgres createdb --owner=sybilcreate sybilcreate";
#    &run_cmd($cmd);
}
else {
    print "Skipped postgres install\n";
}

# Start postgres
if($options{start_postgres}) {
&start_postgres();
}

# Install Mongo
if($options{install_mongo}) {
    &install_mongo();
}
# Start Mongo
if($options{start_mongo}) {
    &start_mongo();
}

# Create Postgres DB
if($options{username} && $options{password} && $options{db_name} && $options{create_db}) {
    &create_db();
    # Create user/database
#    my $cmd = "export PGPASSWORD=createsybils;psql -U sybilcreate -h localhost -c \"CREATE USER $options{username} WITH PASSWORD '$options{password}'\"";
#    &run_cmd($cmd);
#    $cmd = "export PGPASSWORD=createsybils;psql -h localhost -U sybilcreate -c \"CREATE DATABASE $options{db_name} OWNER $options{username}\"";
#    &run_cmd($cmd);
}
else {
    print STDERR "Not creating a database since username password and db_name were not specified\n";
}

# Populate Postgres DB from remote source
if($options{db_url}) {
    # Download database
    my $cmd = "wget -N -P $ROOT/ $options{db_url}";
    &run_cmd($cmd);

    $options{db_url} =~ /\/([^\/]+)$/;
    my $filename = $1;
    
    # Load database
    &load_db("$ROOT/$filename");
}
else {
    print STDERR "Not loading a database because not db_url was specified\n";
}

# Install Sybil
if($options{install_sybil}) {
    &install_sybil();
}

# Configure Sybil site
if($options{config_url}) {
    my $cmd = "wget -N -P $PATH_TO_SYBIL/cgi/conf/ $options{config_url}";
    &run_cmd($cmd);

    $options{config_url} =~ /([^\/]+)\.conf$/;
    my $SITENAME = $1;
}
elsif($options{db_name} && $options{username} && $options{password} && $SITENAME && $options{create_site_config}) {
    &create_config();
}

if($options{cache_data}) {
    # Pre-cache some data
    print "Pre-caching data\n";
#    &run_cmd("cd $PATH_TO_SYBIL/cgi/shared/;sudo -u www-data ./search_proteins.cgi site=$SITENAME user=$options{username} password=$options{password}");
    &run_cmd("cd $PATH_TO_SYBIL/cgi/shared/;./search_proteins.cgi site=$SITENAME user=$options{username} password=$options{password}");

    &run_cmd("mkdir /var/www/sybil/$SITENAME");
    open OUT, ">/var/www/sybil/$SITENAME/index.html";

    print OUT "<META HTTP-EQUIV=\"Refresh\"\n\tCONTENT=\"0; URL=/sybil/current/cgi/shared/index.cgi?site=$SITENAME\">";
    close OUT;

    &run_cmd("perl -pi -e \"s[#\\s+'db_cache' => 'file',][\\t'db_cache' => 'file',]\" $PATH_TO_SYBIL/cgi/conf/$SITENAME.conf");
    print "Complete\n";
}

if($options{create_archive}) {
    &create_archive();
}
if($options{load_archive}) {
    &load_archive();
}

#####################
# Subroutines
#####################

# Generic run command sub.
sub run_cmd {

    my $cmd = shift;

    print `$cmd`;
    if($?) {
        print STDERR "$cmd\n\n$?";
    }
   # print "$cmd\n";
}

sub cache_data {

    # Pre-cache some data
    print "Pre-caching data\n";

    # First turn off caching as we're going to populate the mongo db
    &run_cmd("perl -pi -e \"s[^\\s+'db_cache' => 'file',][#\\t'db_cache' => 'file',]\" $PATH_TO_SYBIL/cgi/conf/$SITENAME.conf");
#    &run_cmd("cd $PATH_TO_SYBIL/cgi/shared/;sudo -u www-data ./search_proteins.cgi site=$SITENAME user=$options{username} password=$options{password}");

    # This puts a bunch of stuff in mong and also writes the protein name lookups
    &run_cmd("cd $PATH_TO_SYBIL/cgi/shared/;./search_proteins.cgi site=$SITENAME user=$options{username} password=$options{password}");

    # Add a shortcut url
    &run_cmd("mkdir /var/www/sybil/$SITENAME");
    open OUT, ">/var/www/sybil/$SITENAME/index.html";

    print OUT "<META HTTP-EQUIV=\"Refresh\"\n\tCONTENT=\"0; URL=/sybil/current/cgi/shared/index.cgi?site=$SITENAME\">";
    close OUT;

    # Turn on file caching
    &run_cmd("perl -pi -e \"s[#\\s+'db_cache' => 'file',][\\t'db_cache' => 'file',]\" $PATH_TO_SYBIL/cgi/conf/$SITENAME.conf");
    print "Complete\n";

}
# Create a postgres db
sub create_db {
    # Create user/database
    my $cmd = "export PGPASSWORD=createsybils;psql -U sybilcreate -h localhost -c \"CREATE USER $options{username} WITH PASSWORD '$options{password}'\"";
    &run_cmd($cmd);
    # Drop the database just in case it exists
    &run_cmd("export PGPASSWORD=createsybils;psql -U sybilcreate -h localhost -c \"DROP DATABASE \\\"$options{db_name}\\\"\"");
    $cmd = "export PGPASSWORD=createsybils;psql -h localhost -U sybilcreate -c \"CREATE DATABASE \\\"$options{db_name}\\\" OWNER $options{username}\"";
    &run_cmd($cmd);
}

# Install postgres
sub install_postgres {
    # Install postgres
    my $cmd = "apt-get -y install postgresql-$PG_VERSION";
    &run_cmd($cmd);

    # Need to prevent postgres from start at boot.
    &run_cmd("find /etc/rc* -name '*postgresql*' -exec rm {} \;");
}

sub start_postgres {

    if(! -d $PG_DATA_DIR || ! -e "$PG_DATA_DIR/PG_VERSION") {
    my $cmd = "mkdir $PG_DATA_DIR;chown postgres $PG_DATA_DIR;chmod 700 $PG_DATA_DIR";
    &run_cmd($cmd);

    $cmd = "/etc/init.d/postgresql-$PG_VERSION stop";
    &run_cmd($cmd);

    $cmd = "pg_dropcluster --stop $PG_VERSION main";
    &run_cmd($cmd);

    $cmd = "pg_createcluster -d $PG_DATA_DIR --start $PG_VERSION main";
    &run_cmd($cmd);
    `sleep 10`;
    my $cmd = "sudo -u postgres psql -c \"CREATE USER sybilcreate WITH PASSWORD 'createsybils' SUPERUSER\";sudo -u postgres createdb --owner=sybilcreate sybilcreate";
    &run_cmd($cmd);
    }
    else {
        &run_cmd("/etc/init.d/postgresql-$PG_VERSION start");
    }
}

sub install_mongo {

    &run_cmd("xstow /usr/local/stow/mongodb-linux-x86_64-1.8.1");
    &run_cmd("mkdir -p $ROOT/sybilmongo");

}

# Start a mongo server on a non-default port
sub start_mongo {

    # Start a mongo database
    &run_cmd("mkdir -p $ROOT/sybilmongo");
    &run_cmd("rm $ROOT/sybilmongo/mongod.lock");
    &run_cmd("mongod --port=10000 --logpath=$ROOT/sybilmongo/log --dbpath=$ROOT/sybilmongo/ --fork");
    sleep 10;
    `mongo -eval "db.currentOp()" --port 10000`;
    if( $? ) {
        &run_cmd("mongod --port=10000 --logpath=$ROOT/sybilmongo/log --dbpath=$ROOT/sybilmongo/ --fork");
    }
}

# Install sybil
sub install_sybil {

    # Install Sybil
    my $cmd = "wget -N -P $ROOT/sybil_install/ $SYBIL_TARBALL";
    &run_cmd($cmd);

    $SYBIL_TARBALL =~ /\/([^\/]+)$/;
    my $fn =$1;
    my $base = $fn;
    $base =~ s/\..*$//;
    $cmd = "cd $ROOT/sybil_install;tar xzvf $fn";
    &run_cmd($cmd);

    # Make the Sybil web directory
    $cmd = "mkdir -p $PATH_TO_SYBIL";
    &run_cmd($cmd);

    # Copy over the cgi/htdocs directories
    $cmd = "cp -r $ROOT/sybil_install/$base/cgi-bin $PATH_TO_SYBIL/cgi";
    &run_cmd($cmd);

    $cmd = "cp -r $ROOT/sybil_install/$base/htdocs $PATH_TO_SYBIL/htdocs";
    &run_cmd($cmd);

    # Move the default Sybil.conf and urls_shared.xml into place
    $cmd = "cp $PATH_TO_SYBIL/cgi/conf/demo.Sybil.conf $PATH_TO_SYBIL/cgi/conf/Sybil.conf";
    &run_cmd($cmd);

    $cmd = "cp $PATH_TO_SYBIL/cgi/shared/urls/clovr.urls_shared.xml $PATH_TO_SYBIL/cgi/shared/urls/urls_shared.xml";
    &run_cmd($cmd);

    # Open permissions on the htdocs directory for purposes of caching
    &run_cmd("chmod 777 /var/www/sybil/current/htdocs");
    &run_cmd("chmod 777 /var/www/sybil");

    # Open permissions on conf file so we can write new confs.
    &run_cmd("chmod 777 /var/www/sybil/current/cgi/conf");

    # Create a tmp area that is on the data volume but that is web accessible
    &run_cmd("mkdir -p $ROOT/sybiltmp/web;chmod 777 $ROOT/sybiltmp;chmod 777 $ROOT/sybiltmp/web/");
    &run_cmd("ln -s $ROOT/sybiltmp/web $PATH_TO_SYBIL/htdocs/tmp");

    # Need to change any cgi scripts that are using /usr/bin/env perl to /usr/bin/env perl
    &run_cmd("perl -pi -e 's[#!/usr/bin/env perl][#!/usr/bin/env perl]' $PATH_TO_SYBIL/cgi/shared/*.cgi");

    # Need to install an older version of Bio::Graphics
    &run_cmd("cpan -i LDS/Bio-Graphics-1.96.tar.gz");

    # Pull down an older version of extjs.
    &run_cmd("wget -N -P /var/www/ http://cb2.igs.umaryland.edu/ext-3.2.1.tgz");
    &run_cmd("cd /var/www/;tar xzvf ext-3.2.1.tgz;ln -s ext-3.2.1 ext");
    
    #HACK to fix the mongo config
    &run_cmd("perl -pi -e 's[tettelin-lx.igs.umaryland.edu][localhost:10000]' $PATH_TO_SYBIL/cgi/Sybil/ChadoMongoSybilDB.pm");

    # Install MongoDB perl module
    &run_cmd("export PERL_MM_USE_DEFAULT=1;cpan -if MongoDB");

}

# Create a config file from a template
sub create_config {

    # Step 1: Create the config file
    my $output_conf = "$CONF_PATH/$SITENAME.conf";
    &run_cmd("cp $CONF_TEMPLATE $output_conf");

    # Step 2: Replace the usual parameters

    ## site
    print `perl -pi -e 's[SITE_NAME][$SITENAME]' $output_conf`;

    ## db
    print `perl -pi -e 's[DB_NAME][$options{db_name}]' $output_conf`;

    ## username/password
    my $username = $options{username} ? "\\'user\\' => \\'$options{username}\\'," : "# No username"; 
    print `perl -pi -e "s[\'USERNAME_PARAM\'][$username]" $output_conf`;

    my $password = $options{password} ? "\\'password\\' => \\'$options{password}\\'," : "# No password"; 
    print `perl -pi -e "s[\'PASSWORD_PARAM\'][$password]" $output_conf`;

    ## cache directory
    # Create a tmp area that is on the data volume but that is web accessible
    &run_cmd("mkdir -p $ROOT/sybiltmp/web;chmod 777 $ROOT/sybiltmp;chmod 777 $ROOT/sybiltmp/web/");
    &run_cmd("ln -s $ROOT/sybiltmp/web $PATH_TO_SYBIL/htdocs/tmp");
    if(!$options{local_conf} && ! -e "$options{cache_dir}/$options{site}") {
        `mkdir $options{cache_dir}/$options{site}`;
    }
    print `perl -pi -e 's[CACHE_DIR][$options{cache_dir}]' $output_conf`;

    ## server
    print `perl -pi -e 's[SERVER][$SERVER]' $output_conf`;

    ## schema
    print `perl -pi -e 's[SCHEMA][$options{schema}]' $output_conf`;

    ## generate organism colors
    # Non-conf path
    $CONF_PATH =~ /^(.*\/)conf/;
    print "perl -I $1 -I $PATH_TO_SYBIL/cgi/shared/ $PATH_TO_SYBIL/cgi/shared/makeOrganismColors.cgi site=$SITENAME";
    my $colors = `perl -I $1 -I $PATH_TO_SYBIL/cgi/shared/ $PATH_TO_SYBIL/cgi/shared/makeOrganismColors.cgi site=$SITENAME`;

    $colors =~ s/\n/\\n/g;
    $colors =~ s/"/\"/g;
    $colors =~ s/'/#/g;
    print "$colors\n";
    print `perl -pi -e 's[.ORGANISM_COLORS.][$colors]' $output_conf`;

    ## image directory/url
    if(!$options{local_conf} && ! -e "$options{image_dir}/$SITENAME") {
        `mkdir -m 777 $options{image_dir}/$SITENAME`;
    }

    print `perl -pi -e 's[IMAGE_DIR][$options{image_dir}/$SITENAME]' $output_conf`;
    print `perl -pi -e 's[IMAGE_URL][$options{image_url}/$SITENAME]' $output_conf`;

    ## java/clustalw path
    chomp $JAVA_PATH;
    chomp $CLUSTALW_PATH;
    my $java_path = $JAVA_PATH ? "\\'path_to_java\\' => \\'$JAVA_PATH\\'," : '# default java path';
    my $clustalw_path = $CLUSTALW_PATH ? "\\'clustalw_path\\' => \\'$CLUSTALW_PATH\\'," : '# default clustalw path';
    print `perl -pi -e "s[\'JAVA_PATH\'][$java_path]" $output_conf`;
    print `perl -pi -e "s[\'CLUSTALW_PATH\'][$clustalw_path]" $output_conf`;

}

# Create a sybil archive that can be restored.
sub create_archive {

    # create a dump directory
    &run_cmd("mkdir $options{output_dir}/$SITENAME\_sybil");

    # dump the database
    &run_cmd("export PGPASSWORD=$options{password};pg_dump -h localhost -U $options{username} $options{db_name} > $options{output_dir}/$SITENAME\_sybil/$options{db_name}.dump");

    # compress the dump
    &run_cmd("gzip -f $options{output_dir}/$SITENAME\_sybil/$options{db_name}.dump");

    # copy the config file
    &run_cmd("cp $CONF_PATH/$SITENAME.conf $options{output_dir}/$SITENAME\_sybil/");

    # tar it up
    &run_cmd("tar -c -C $options{output_dir} -zvf $options{output_dir}/$SITENAME\_sybil.tgz $SITENAME\_sybil");
}

# Load a gzipped postgres dump into a database
sub load_db {

    my $db_file = shift;
    print STDERR "$db_file was the database archive\n";
    # Load database
    my $cmd = "export PGPASSWORD=$options{password};zcat $db_file | psql -h localhost -d $options{db_name} -U $options{username}";
    &run_cmd($cmd);


}

# Tag a sybil archive and deploy it.
sub load_archive {


    die "Need to specify an output directory\n" if(!$options{output_dir});
    # Unpack the archive
    my $archive = $options{load_archive};


    # Create a tmp area that is on the data volume but that is web accessible
    &run_cmd("mkdir -p $ROOT/sybiltmp/web;chmod 777 $ROOT/sybiltmp;chmod 777 $ROOT/sybiltmp/web/");
    &run_cmd("ln -s $ROOT/sybiltmp/web $PATH_TO_SYBIL/htdocs/tmp");

    if($options{archive_is_tag}) {
        my @lines = `vp-describe-dataset --tag-name=$options{load_archive}`;
        foreach my $line (@lines) {
        if($line =~ /FILE\s+(.*)/) {
		$archive = $1;
		chomp $archive;
	}
        elsif($line =~ /METADATA\s+compressed_file\s+(.*)/) {
            my $tarball = $1;
            my ($name,$path,$suff) = fileparse($tarball,('.tar.gz','.tgz'));
            &run_cmd("tar xCvf $options{output_dir} $tarball");
            $archive = `find $options{output_dir}/$name/ -name '*.*'`;
            chomp $archive;
            print STDERR "Found $archive\n";
        }
  
	}

    }
    # Looks like we might have a list instead
    elsif( $archive !~ /.tgz|.tar.gz/  ) {
        $archive = `head -1 $options{load_archive}`;
        chomp $archive; 
    }

    if(-e $archive) {
        &run_cmd("tar xCzvf $options{output_dir} $archive");
    } 
    else { die "Couldn't find archive $options{load_archive} $archive\n";}


    # find the conf file
    my $config = `find $options{output_dir}/ -name '*.conf'`;
    chomp $config;

    # Find the tarball
    my $db_archive = `find $options{output_dir}/ -name '*dump.gz'`;
    chomp $db_archive;

    # Get the sitename
    $config =~ /([^\/]+)\.conf$/;
    $SITENAME = $1;
    `echo $SITENAME > $options{output_dir}/sitename`; 
    # Copy the config file into place
    &run_cmd("cp $config $CONF_PATH/");
    
    # Get the db_name, username and password. Very HACKy
    open IN, "<$config" or die "Couldn't open $config\n";
    while(<IN>) {
        if(/'db' => '(\S+)',/) {
            $options{db_name} = $1;
        }
        elsif(/'user' => '(\S+)',/) {
            $options{username} = $1;
        }
        elsif(/'password' => '(\S+)',/) {
            $options{password} = $1;
        }
    }
    close IN;

    # Create the database
    &create_db();

    # Load the database
    &load_db($db_archive);
    
    # Cache data
    &cache_data();

    print STDERR "Deployed $SITENAME\n";
}
