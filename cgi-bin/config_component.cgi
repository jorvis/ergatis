#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);

use Tree::DAG_Node;
use Config::IniFiles;

my $conffile = param('conffile');
my $outputfile = param('outputfile');
my $dowrite = param('dowrite');
my $limitsect = param('limitsect');
my $ignoresect = param('ignoresect');
my $sharedconf = param('sharedconf');
my $id = param('id');

if($dowrite){
    my $newcfg = new Config::IniFiles(-file => $conffile);
    print header();
    &write_config_from_params($newcfg,$outputfile);
    print "<html><body onLoad='window.parent.location.reload();'>Output file $outputfile written.<br><pre>\n";
    &print_config($newcfg);
    print "</pre><br><a href='javascript:window.close();'>[close]</a></body></html>";
    exit;
}

my $cfg = new Config::IniFiles(-file => $conffile);

#print all conf options

my $ignoresectlookup = {};
foreach my $elt (split (/,/,$ignoresect)){
    $ignoresectlookup->{$elt} = 1;
}
my $limitsectlookup = {};
foreach my $elt (split (/,/,$limitsect)){
    $limitsectlookup->{$elt} = 1;
}

my $currcfg = &import_includes($cfg,$sharedconf);

print header();

print "<html><body>";
print "<form name='myform'><input type=hidden name='dowrite' value=1><input type=hidden name='conffile' value='$conffile'><input type=hidden name='outputfile' value='$outputfile'>";
print "Configuration will be written to $outputfile&nbsp;<a href='remove_component.cgi?conffile=$conffile'>[remove]</a>\n";

if($id ne ""){
    $currcfg->setval('init','$;PIPELINE$;',$currcfg->val(($currcfg->GroupMembers("workflowdocs"))[0],'$;NAME$;')."_$id");
}

my @sections = $currcfg->Sections();

foreach my $section (@sections){
    if(($limitsect && $limitsectlookup->{$section}) || !$limitsect){
	if($ignoresectlookup->{$section}){
	    print "<h3>$section</h3><table>";
	}
	else{
	    print "<h3>$section</h3><table>";
	}
	my @parameters = $currcfg->Parameters ($section);
	foreach my $param (@parameters){
	    my $value = $currcfg->val($section,$param);
	    $value =~ s/\s//g;
	    if($ignoresectlookup->{$section}){
		print "<tr><td>$param</td><td>$value</td></tr>";
	    }
	    
	    else{ 
		print "<tr><td>$param</td><td><input type=text name='edit_",$section,"::","$param' value='$value' size=100></td></tr>";
	    }
	}
    }
    print "</table>";
}
print "<input type='submit' value='Write config'></form></html></body>";

sub import_includes{
    my($cfg,$sharedconf) = @_;
    my @includes = $cfg->GroupMembers("include");
    my $currcfg = $cfg;
    foreach my $member (@includes){
	my @parameters = $cfg->Parameters ($member);
	foreach my $param (@parameters){
	    if($param eq '$;SHARED_CONFIG$;'){
		my $includefile =  $sharedconf;#$cfg->val($member,$param);
		if(-e $includefile){
		    my $newcfg = new Config::IniFiles( -file => $includefile, 
						       -import => $currcfg);
		    $newcfg->setval($member,$param,$includefile);
		    $currcfg = $newcfg;
		}
	    }
	}
    }
    return $currcfg;
}

sub write_config_from_params{
    my($newcfg,$filename) = @_;
    my @allparams = param();
    foreach my $p (@allparams){
	if($p =~ /^edit_/){
	    #print "Parsing $p<br>\n";
	    my $pval = $p;
	    $pval =~ s/^edit_//;
	    my($section,$secp) = split(/\:\:/,$pval);
	    #print "Writing [$section] $secp=",param($p),"<br><\n>";
	    my $s = $newcfg->setval($section,$secp,param($p));
	    if(!$s){
		die "Can't write  [$section] $secp".param($p);
	    }
	}
    }
    $newcfg->WriteConfig($filename);
    `chmod 666 $filename`;
}

sub print_config{
    my($cfg) = @_;
    my @sections = $cfg->Sections();
    foreach my $section (@sections){
	print "[$section]\n";
	my @parameters = $cfg->Parameters ($section);
	foreach my $param (@parameters){
	    my $value = $cfg->val($section,$param);
	    print "$param=$value\n";
	}
    }
}


    
