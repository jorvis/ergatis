<html>
	<head>
		
		<title><?php echo $page_title; ?> - Ergatis</title>
		<link rel="stylesheet" type="text/css" href="./css/common.css">
                <link rel="stylesheet" type="text/css" href="./css/header.css">
                <link rel="stylesheet" type="text/css" href="./css/monitor.css">
                <link rel="stylesheet" type="text/css" href="./css/index.css">
                <link rel="stylesheet" type="text/css" href="./css/forms.css">
                <link rel="stylesheet" type="text/css" href="./css/config_file.css">
                
                <?php
 			if ( isset($extra_head_tags)) {               
                		echo $extra_head_tags;
                	}
                     ?>
                
	</head>
	<body>
		<div id='page_header'>
			<div id='title_container'>
        			<span id='page_title'>Ergatis - Prokaryotic Pipeline Builder</span>
	 		</div>
		</div> 
		<div id='primary_nav_container'>
			<ul id='primary_nav_elements'>
				<li id='nav_pipelines'><a href="http://ergatis.igs.umaryland.edu">Ergatis</a></li>
				<li id='nav_documentation'><a href="help.php">Help</a></li>
			</ul>
		</div>
		<div id='sub_nav_container'>
			<ul id='sub_nav_elements'>
				<li><a href='prok_pipeline_step1.php'>Home</a></li>
			</ul>
		</div>
		
		<div id='page_container'>
			<div id='content_container'>
			

