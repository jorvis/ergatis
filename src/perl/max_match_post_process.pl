#!/usr/local/bin/perl


use lib("shared");
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use English;
use File::Path;

my %options = ();
my $results = GetOptions (\%options, 'db_dir=s', 'pep_dir=s', 'mini_db|i=s', 'help|h', 'project|p=s');
###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $mini_db   = $options{'mini_db'};
my $db_dir    = $options{'db_dir'};
$db_dir =~ s/\/+$//;       #remove terminating '/'s
my $pep_dir   = $options{'pep_dir'};
$pep_dir =~ s/\/+$//;       #remove terminating '/'s
my $project   = $options{'project'};

if(!$project or !$mini_db or !$db_dir or !$pep_dir or exists($options{'help'})) {
    &print_usage();
}

###-------------------------------------------------------###

if(! -d $db_dir ) {
    mkdir $db_dir;
}
chmod 0777, $db_dir;


if(! -d $pep_dir ) {
    mkdir $pep_dir;
}else {
    rmtree([<$pep_dir/*>], 0, 1); 
}
chmod 0777, $pep_dir;



unless(open (MINIDB, "< $mini_db")) {
    print STDERR "Unable to Open \"$mini_db\" due to $!\n";
    exit 5;
}else {
    my $pep_db = "$db_dir/$project.pep";
    open(DB, ">$pep_db") || die "Cant open $pep_db due to $!";
    my $line = <MINIDB>;
    my ($sequence, $header);
    while(defined($line)) {
        $sequence = undef;
        if($line =~ /^>/) {
	    $header = $line;
            chomp($header);
	    $header =~ s/>//;
            my ($gene, $asmbl_id) = split(/\s+/, $header);
	    while(defined($line=<MINIDB>) and $line !~ /^>/ ) {
		next if($line =~/^\s+$/);                   #skip blank lines
                #chomp($line);
                $sequence .= $line;
            }
            print DB ">$asmbl_id:$gene\n$sequence";
            
	    my $each_pep_dir = "$pep_dir/$asmbl_id";
	    if(! -d $each_pep_dir ) {
		mkdir $each_pep_dir;
		chmod 0777, $each_pep_dir;
	    }	

	    my $each_pep ="$each_pep_dir/$gene.fsa"; 
	    open (EACH_PEP, "> $each_pep") or die "Unable to write to $each_pep du to $!";
	    print EACH_PEP ">$gene\n$sequence";
            close EACH_PEP;
        } else {
	    $line = <MINIDB>;
        }
    }
    close DB;
    qx(setdb $pep_db);
    chmod 0777, <$pep_db*>;

}














