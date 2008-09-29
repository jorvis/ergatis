#!/usr/bin/perl

#./gff2fasta.pl [featuretype] [keysfile]
#featuretype is stored in column 3 of gff3 
#For example, this is a line with featuretype cds
#nmpdr|158878.1.contig.NC_002758 NMPDR   cds ...
#
use strict;
use CGI qw/escape unescape/;

my $features = {};

my $fastadirective=0;

my $savefastaentry=0;
my $currfeature;
my @pepfastaoutbuffer;
my @seqfastaoutbuffer;

while(my $line=<STDIN>){
    $line = unescape($line);
    if($line =~ /^\#\#FASTA/){
	$fastadirective=1;
    }
    if($line !~ /^\#/){
	if(!$fastadirective){
	    #parse tab delimeted
	    chomp $line;
	    my (@elts) = split(/\t/,$line);
	    $elts[2] = lc($elts[2]);
	    if($elts[2] eq "$ARGV[0]"){
		my(%attrs) = split(/[;=]/,$elts[8]);
		my(%dbxrefs) = split(/[,:]/,$attrs{'Dbxref'});
		$features->{$attrs{'ID'}}->{'save'}=1;
		$features->{$attrs{'ID'}}->{'description'} = $attrs{'description'};
		$features->{$attrs{'ID'}}->{'center'} = $elts[1];
		$features->{$attrs{'ID'}}->{'genomic_source'} = $elts[0];
		$features->{$attrs{'ID'}}->{'taxon'} = $features->{$elts[0]}->{'taxon'};
		$features->{$attrs{'ID'}}->{'type'} = $elts[2];
	    }
	    elsif($elts[2] eq 'contig'){
		my(%attrs) = split(/[;=]/,$elts[8]);
		my(%dbxrefs) = split(/[,:]/,$attrs{'Dbxref'});
		$features->{$attrs{'ID'}}->{'save'}=1;
		$features->{$attrs{'ID'}}->{'type'} = $elts[2];
		$features->{$attrs{'ID'}}->{'taxon'} = $dbxrefs{'taxon'};
		$features->{$attrs{'ID'}}->{'organism_name'} = $attrs{'organism_name'};
		$features->{$attrs{'ID'}}->{'strain'} = $attrs{'strain'};
		#print STDERR "$dbxrefs{'taxon'}\t$attrs{'organism_name'} $attrs{'strain'}\t$attrs{'ID'}\n";
	    }
	}
	else{
	    #parse fasta
	    if($line =~ /^>(\S+)/){
		if($features->{$1}->{'save'}==1){
		    $savefastaentry=1;
		    $currfeature=$1;
		    #use counter to ensure we only print out each entry once.
		    $features->{$currfeature}->{'save'}=0;

		}
		else{
		    $savefastaentry=0;
		}
	    }
	    if($savefastaentry){
		if($features->{$currfeature}->{'type'} eq 'cds'){
		    push @pepfastaoutbuffer,$line;
		}
		elsif($features->{$currfeature}->{'type'} eq 'contig'){
		    push @seqfastaoutbuffer,$line;
		}
		if($line =~ /\*/){
		    print "ERROR:$currfeature $features->{$currfeature}->{'description'}\n";
		}

	    }
	}
    }
}



if(defined $ARGV[1]){
    open FILE,">>$ARGV[1]" or die;
    foreach my $id (keys %$features){
	if($features->{$id}->{'type'} eq $ARGV[0]){
	    print FILE "$id $features->{$id}->{'center'},$features->{$id}->{'taxon'},$features->{$id}->{'genomic_source'},$features->{$id}->{'description'}\n";
	}
    }
    close FILE;
}

if(defined $ARGV[2]){
    open FILE,">>$ARGV[2]" or die;
    print FILE @seqfastaoutbuffer;
    close FILE;
}

print @pepfastaoutbuffer;
