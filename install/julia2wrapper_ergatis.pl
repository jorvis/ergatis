use File::Basename;
use File::Find;

sub get_julia_bins{
    my($instdir,$workflowdocsdir,$schemadocsdir);
    my(@binfiles);
    my($wrapper_str);
    my $julia_path = `which julia`;  ## can be overwritten below
    chomp $julia_path; 

    foreach my $arg (@ARGV){
        if($arg =~ /INSTALL_BASE/){
            ($instdir) = ($arg =~ /INSTALL_BASE=(.*)/);
        }
        if($arg =~ /WORKFLOW_DOCS_DIR/){
            ($workflowdocsdir) = ($arg =~ /WORKFLOW_DOCS_DIR=(.*)/);
        }
        if($arg =~ /SCHEMA_DOCS_DIR/){
            ($schemadocsdir) = ($arg =~ /SCHEMA_DOCS_DIR=(.*)/);
        }
        if($arg =~ /JULIA_PATH/){
            ($julia_path) = ($arg =~ /JULIA_PATH=(.*)/);
        }
    }

    open FILE, 'MANIFEST' or die "MANIFEST is missing!\n";
    open SYMS, "+>README.symlinks" or die "Can't save symlinks for silly sadmins";
    print SYMS "#Copy or symlink the following shell scripts into a standard area\n";

    my $envbuffer;
    my $env_hash = {'WORKFLOW_DOCS_DIR' => "$workflowdocsdir",
		    'SCHEMA_DOCS_DIR' => "$schemadocsdir",
		    'WORKFLOW_WRAPPERS_DIR'  => "$instdir/bin"
		    };
    
    foreach my $key (keys %$env_hash){
	$envbuffer .= "if [ -z \"\$$key\" ]\nthen\n    $key=$env_hash->{$key}\nexport $key\nfi\n";
    }


    while(my $line = <FILE>){
	chomp $line;
	if($line =~ m|julia/[\w-]+\.jl| ){
	    my($fname) = basename($line);
	    my($strip_fname) = ($fname =~ /(.*)\.jl$/);
	    open WRAPPER, "+>bin/$strip_fname" or die "Can't open file bin/$strip_fname\n";
	    my($shell_args)  = q/"$@"/;
	    my $addbuffer = $envbuffer;
            print WRAPPER <<_END_WRAPPER_;
#!/bin/sh
$addbuffer

umask 0000

unset PERL5LIB
unset LD_LIBRARY_PATH

LANG=C
export LANG
LC_ALL=C
export LC_ALL

PERL_MOD_DIR=$instdir/lib/5.8.8
export PERL_MOD_DIR

export PERL5LIB=$instdir/lib/perl5/

    $julia_path $instdir/bin/$fname $shell_args    

_END_WRAPPER_
   ;
	    close WRAPPER;
	    
	    print SYMS "$instdir/bin/$strip_fname\n";
	    
	    push @binfiles,"$line";
	    push @binfiles,"bin/$strip_fname";
	    $wrapper_str .= "bin/$strip_fname ";
	}
    }
    close SYMS;
    close FILE;
    return (\@binfiles,$wrapper_str);
}


return 1;
