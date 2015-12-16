<?php
	$study_id = $_GET['study_id'];
	$client = $_GET['client'];

	$id = array();			
	$command = "/usr/bin/perl ./perl/get_study_stage.pl $client $study_id";
	echo exec($command, $id);
?>	
