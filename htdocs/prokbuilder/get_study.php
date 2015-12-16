<?php
	$sample_id = $_GET['sample_id'];

	$id = array();			
	$command = "/usr/bin/perl ./perl/get_study_from_sample.pl $sample_id";
	echo exec($command, $id);
?>	
