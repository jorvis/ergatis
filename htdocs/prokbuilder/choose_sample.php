<?php
	$proj_id = $_GET['proj_id'];
	$client = $_GET['client'];

	$id = array();			
	$command = "/usr/bin/perl ./perl/IPD_to_list.pl $client $proj_id";
	exec($command, $id);
	$string = '';
	$result = array();
	foreach ($id as $index => $key) {
		preg_match('/\s+\(ID:\s+(\d+)\)/', $key, $match);
		$value = $match[1];
		$result[$key] = $value;	#want to transform into array where IDs are values, and exec output is the keys
	}
	echo json_encode($result);		#return array in JSON format
?>					
