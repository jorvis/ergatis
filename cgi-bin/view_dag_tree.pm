package view_dag_tree;

use strict;
use Tree::DAG_Node;
use Exporter;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT= qw(&define_dag_html &open_root_nodes &open_new_node);

sub new { 
    my $classname = shift;
    my $self = {};
    bless($self,$classname);
    $self->{DEBUG} = 0;
    $self->{dag_file} = undef;
    $self->{link_root} = undef;
    $self->_init(@_);
    if(! $self->{dag_obj}){
	print STDERR "No dag specified. Retrieving from $self->{dag_file}\n";
	$self->{dag_obj} = &LockRetrieve($self->{dag_file});
    }
    $self->{new_dag_file} = $self->{dag_file};
    $self->{link_url_root} = &ParseLinkRoot($self->{link_root});
    $self->{link_url} =  $self->{link_url_root}."?&xmltemplate=$self->{xmltemplate}";
    print STDERR "Writing to $self->{new_dag_file} ".scalar($self->{dag_obj}->descendants())."\n";
    return $self;
}
sub _init {
    my $self = shift;
    if(@_){
	my %extra = @_;
	@$self{keys %extra} = values %extra;
    }
}
sub ParseLinkRoot {
    my($link_tmp) = @_;
    my @a = split(/\//,$link_tmp);
    my $num = scalar(@a);

    return($a[$num - 1]);
}
sub define_dag_html {
    my($self) = @_;
    my $html = $self->make_all_property_html;
    &LockStore($self->{dag_obj},$self->{new_dag_file});
    return($html);
}
sub open_root_nodes {
    my($self) = @_;
    my @daughters = &get_daughters($self->{dag_obj});
    for(my $i = 0; $i < scalar(@daughters); $i++) {
	my $daughter = $daughters[$i];
	&SetNodeAsOpen($daughter);
    }
}

sub open_nodes{
    my($self,$node) = @_;
    my @daughters = &get_daughters($node);
    for(my $i = 0; $i < scalar(@daughters); $i++) {
	my $daughter = $daughters[$i];
	&SetNodeAsOpen($daughter);
	$self->open_nodes($daughter);
    }
}

sub close_nodes{
    my($self,$node) = @_;
    &SetNodeAsClosed($node);
    my @daughters = &get_daughters($node);
    for(my $i = 0; $i < scalar(@daughters); $i++) {
	my $daughter = $daughters[$i];
	&SetNodeAsClosed($daughter);
	$self->close_nodes($daughter);
    }
}

sub open_new_node {
    my($self,$open_node) = @_;

    my @daughters = &get_daughters($self->{dag_obj});
    for(my $i = 0; $i < scalar(@daughters); $i++) {
	my $daughter = $daughters[$i];
	my $key = $daughter->name();

	if($key eq $open_node) {
	    &open_close($daughter);
#	    $self->open_nodes($daughter);
	    last;
	}else{
	    ##Walk down the dag till find the _key node and open daughters
	    $self->find_node_in_daughters($daughter,$open_node);
	}
    }
}
sub find_node_in_daughters {
    my($self,$mother,$open_node) = @_;
    my @daughters = &get_daughters($mother);
    for(my $i = 0; $i < scalar(@daughters); $i++) {
	my $daughter = $daughters[$i];
	my $key = $daughter->name();

	if($key eq $open_node) {
	    &open_close($daughter);
#$self->open_nodes($daughter);
	    last;
	}else{
	    $self->find_node_in_daughters($daughter,$open_node);
	}
    }
}
sub open_close {
    my($node) = @_;
    my $daughters_open = &TestIfAlreadyOpen($node);
    if($daughters_open == 1) {
	&SetDaughterNodesAsClosed($node);
    }else{
	&SetDaughterNodesAsOpen($node);
    }
}
sub TestIfAlreadyOpen {
    my($node) = @_;
    my $already_open = 0;

    my @daughters = &get_daughters($node);
    for(my $i = 0; $i < scalar(@daughters); $i++) {
	my $daughter = $daughters[$i];
	my $open = $daughter->attributes->{'open'};
	$already_open = 1 if($open == 1);
    }
    return $already_open;
}
sub SetDaughterNodesAsClosed {
    my($node) = @_;
    my @daughters = &get_daughters($node);
    for(my $i = 0; $i < scalar(@daughters); $i++) {
	my $daughter = $daughters[$i];
	&SetNodeAsClosed($daughter);
	&SetDaughterNodesAsClosed($daughter);
    }
}
sub SetDaughterNodesAsOpen {
    my($node) = @_;
    my @daughters = &get_daughters($node);
    for(my $i = 0; $i < scalar(@daughters); $i++) {
	my $daughter = $daughters[$i];
	&SetNodeAsOpen($daughter);
	SetDaughterNodesAsOpen($daughter);
    }
}
sub SetNodeAsClosed {
    my($node) = @_;
    $node->attributes->{'open'} = 0;
}
sub SetNodeAsOpen {
    my($node) = @_;
    $node->attributes->{'open'} = 1;
}
sub check_have_daughters {
    my($node) = @_;
    my $has_daughters = 0;
    my @daughters = &get_daughters($node);
    $has_daughters = 1 if(scalar(@daughters) > 0);
    return($has_daughters);
}
sub make_all_property_html {
    my($self) = @_;
    my %html_hash;
    $html_hash{'max_y'} = 0;
    $html_hash{'max_x'} = 0;
    my $y_place = 0;

    #these are the major nodes below root
    my @daughters = &get_daughters($self->{dag_obj});

    for(my $i = 0; $i < scalar(@daughters); $i++) {
	my $x_place = 0;
	my $daughter = $daughters[$i];

	if($daughter->attributes->{'open'} == 1) {
	    $html_hash{$y_place}->{$x_place}->{'label'} = $self->get_label($daughter,$self->{link_url},&check_have_daughters($daughter));
	    $html_hash{'max_y'} = $y_place;
	    $html_hash{'max_x'} = $x_place if($x_place > $html_hash{'max_x'});
	    $y_place++;
	    $y_place = $self->add_daughters($x_place,$y_place,$daughter,\%html_hash, $self->{link_url});
	}
    }
    print "<table border=0>\n";
    for(my $y = 0; $y <= $html_hash{'max_y'}; $y++) {
	print "<tr border=0>\n";
	for(my $x = 0; $x <= $html_hash{'max_x'}; $x++) {
	    my $text = "\t";
	    if(exists $html_hash{$y}->{$x}->{'label'}) {
		$text = $html_hash{$y}->{$x}->{'label'};
	    }
	    print "<td border=0>$text</td>\n";
	    
	}
    	print "</tr>\n";
    }
    print "</table>";
    return 1;
}

sub get_label {
    my($self,$node,$url,$has_daughters) = @_;
    my $elt = $node->attributes->{'elt'};
    my $xmlfile = $node->attributes->{'xmlfile'};
    my $state = $node->attributes->{'state'};
    my $name = $node->attributes->{'label'};
    my $nodeid = $node->name();
    my $image;
    my $label;
    if($elt eq "commandSet"){
	my $type = $node->attributes->{'type'};
		
	$label = "$name<a href='show_commandset.cgi?xmltemplate=$xmlfile&node=$nodeid'>[info]</a><a href='show_pipeline.cgi?xmltemplate=$xmlfile'>[view]</a>";
	if($self->{'editmode'}){
	    my $left_sister = $node->left_sister;
	    my $right_sister = $node->right_sister;
	    my $prevnode;
	    my $nextnode;
	    if($left_sister){
		$prevnode = $left_sister->name();
	    }
	    if($right_sister){
		$nextnode = $right_sister->name();
	    }
	    
	    $label .= "<a href='edit_commandset.cgi?xmltemplate=$xmlfile&node=$nodeid'>[edit]</a><br>";
	    $label .= "<a href='move.cgi?xmltemplate=$xmlfile&node=$nodeid&target=$prevnode&location=before'>[<--]</a>" if($prevnode);
	    $label .= "<a href='move.cgi?xmltemplate=$xmlfile&node=$nodeid&target=$nextnode&location=after'>[-->]</a>" if($nextnode);
	    if($nodeid =~ /^component/){
		$label .= "<a href='remove_commandset.cgi?xmltemplate=$xmlfile&node=$nodeid'>[-]<a href='add.cgi?xmltemplate=$xmlfile&node=$nodeid&location=after'>[+]</a></a>";
	    }
	    else{
		$label .= "<a href='remove_commandset.cgi?xmltemplate=$xmlfile&node=$nodeid'>[-]</a><a href='add.cgi?xmltemplate=$xmlfile&node=$nodeid&location=last_child'>[+]</a>";
	    }
	}
	
	my $imaget;
	if($type =~ /serial/){
	    $imaget = &get_commandset_image($state);
	}
	elsif($type =~ /parallel/){
	    $imaget = &get_commandset_image($state);
	    
	}
	$image = "$imaget<BR>$imaget";
    }
    elsif($elt eq "command"){
	$label = "$name<a href='show_command.cgi?xmltemplate=$xmlfile&node=$nodeid'>[info]</a>";
	my $log = $node->attributes->{'log'};
	my $jobid = $node->attributes->{'jobid'};
	my $conf = $node->attributes->{'conf'};
	if(-e $log){
	    $label .= "<a href='show_file.cgi?&file=$log' target='_condorlog'>[log]</a><a href='http://htcworker1:8080/antware/htcservice/html/request_display.jsp?RequestID=$jobid' target='_condorlog'>[htcinfo]</a>";
	    $label .= "<a href='add.cgi?xmltemplate=$xmlfile&node=$nodeid&location=after'>[add]</a>" if($self->{'editmode'});
	}
	if(-e $conf){
	    $label .= "<a href='show_file.cgi?&file=$conf' target='_condorlog'>[conf]</a>";
	}
	
	my $imaget = &get_command_image($state);
	$image = "$imaget";
    }


    my $key = $node->name();
    my $link = $url . "&open=$key";
    my $click_link;
    if($has_daughters==1){
	$click_link = "<a href=\"$link\">$image</a>$label";
    }
    else{
	$click_link = "<a href=\"$link\">$image</a>$label";
    }
    return($click_link);
}
sub add_daughters {
    my($self,$x_place,$y_place,$mother,$html_ref,$url) = @_;
    $x_place++;
    if($self->{'editmode'} && ($mother->name() =~ /^component_/)){
	$self->close_nodes($mother);
    }
    my @daughters = &get_daughters($mother);
    for(my $i = 0; $i < scalar(@daughters); $i++) {
	my $daughter = $daughters[$i];
	if(($daughter->attributes->{'open'} == 1)){
	    if(($self->{'max_level'} eq '') || ($x_place < $self->{'max_level'})) {
		$html_ref->{$y_place}->{$x_place}->{'label'} = $self->get_label($daughter,$url,&check_have_daughters($daughter));
		$html_ref->{'max_y'} = $y_place;
		$html_ref->{'max_x'} = $x_place if($x_place > $html_ref->{'max_x'});
		$y_place++;
		$y_place = $self->add_daughters($x_place,$y_place,$daughter,$html_ref,$url);
	    }
	    else{
		$self->close_nodes($daughter);
	    }
	}
    }
    return($y_place);
}
sub get_daughters {
    my($cnode) = @_;
    my @daughters = $cnode->daughters;
    return(@daughters);
}
sub LockStore {
    my($data_ref,$file) = @_;
    print STDERR "Writing tree with daughters ".scalar($data_ref->descendants)." to $file\n";
    use Storable qw(lock_store);  #to serialize
    lock_store($data_ref,$file);
}

sub LockRetrieve {
    my($file) = @_;
    use Storable qw(lock_retrieve);  #to serialize
    my $data_ref = undef;
    $data_ref = lock_retrieve($file) if(-r $file);
    return($data_ref);
}


sub get_command_image{
    my ($status) = @_;
    if($status eq 'complete'){
	return "<img width=5 height=15 border=0 src = \"/papyrus/htdocs/greenpixel.gif\">";
    }
    if($status eq 'interrupted'){
	return "<img width=5 height=15 border=0 src = \"/papyrus/htdocs/purplepixel.gif\">";
    }
    if($status eq 'incomplete'){
	return "<img width=5 height=15 border=0 src = \"/papyrus/htdocs/greypixel.gif\">";
    }
    if($status eq 'pending'){
	return "<img width=5 height=15 border=0 src = \"/papyrus/htdocs/yellowpixel.gif\">";
    }
    if($status eq 'running'){
	return "<img width=5 height=15 border=0 src = \"/papyrus/htdocs/bluepixel.gif\">";
    }
    if($status eq 'running_distributed'){
	return "<img width=5 height=15 border=0 src = \"/papyrus/htdocs/ltbluepixel.gif\">";
    }
    if($status eq 'error' || $status eq 'failed'){
	return "<img width=5 height=15 border=0 src = \"/papyrus/htdocs/redpixel.gif\">";
    }
    return "<img width=5 height=15 border=0 src = \"/papyrus/htdocs/blackpixel.gif\">";
}

sub get_commandset_image{
    my ($status) = @_;
    if($status eq 'complete'){
	return "<img width=20 height=5 border=1 src = \"/papyrus/htdocs/greenpixel.gif\">";
    }
    if($status eq 'interrupted'){
	return "<img width=20 height=5 border=1 src = \"/papyrus/htdocs/purplepixel.gif\">";
    }
    if($status eq 'incomplete'){
	return "<img width=20 height=5 border=1 src = \"/papyrus/htdocs/greypixel.gif\">";
    }
    if($status eq 'pending'){
	return "<img width=20 height=5 border=1 src = \"/papyrus/htdocs/yellowpixel.gif\">";
    }
    if($status eq 'running'){
	return "<img width=20 height=5 border=1 src = \"/papyrus/htdocs/bluepixel.gif\">";
    }
    if($status eq 'error' || $status eq 'failed'){
	return "<img width=20 height=5 border=1 src = \"/papyrus/htdocs/redpixel.gif\">";
    }
    return "<img width=20 height=5 border=1 src = \"/papyrus/htdocs/blackpixel.gif\">";
}
