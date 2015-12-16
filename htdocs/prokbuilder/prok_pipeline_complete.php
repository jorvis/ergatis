<?php
$page_title = 'Pipeline Builder - Complete';
include_once('header.php');
/*
<html>
	<head>
	</head>
	<body>
		<div id='page_container'>

			<div id='content_container'>
*/  ?>
                <?php
                    include_once('prok_pipeline_step2_submit.php');
					if ($errFlag == 0) {
						$pipeline_url = `/usr/bin/perl ./perl/run_prok_pipeline.pl --layout {$formValuesArr['layout_file']['default']} --config {$formValuesArr['config_file']['default']} --repository_root {$formValuesArr['rep_root']['default']}`;

						if (preg_match('/pipeline_id/',$pipeline_url,$matches)) {
							$res = array();
							$temp = preg_split("/\|/",$pipeline_url);
							foreach ($temp as $t) {
								list($label, $val) = preg_split("/\-\>/",$t);
								trim($label);
								$res[$label] = trim($val);
							}
							echo "<br>";
							echo "<h2>STEP FINAL : Pipeline {$res['pipeline_id ']} created successfully<sup><a href=\"./help.php\">?</a></sup></h2>";
							echo "<h3>Following files have been created :</h3>";
							echo "<ul>";
							echo "<li>Pipeline Configuration File : {$formValuesArr['config_file']['default']}</li>";
							echo "<li>Pipeline Layout File : {$formValuesArr['layout_file']['default']}</li>";
							echo "<li>Repository Root : {$formValuesArr['rep_root']['default']}";
							echo "</ul>";
							echo "<h3>Click the link below to view and run the pipeline in Ergatis. Hit rerun to start the pipeline</h3>";
							echo "<font color=\"blue\">Note : Open the link below in a browser. The pipeline is only <b>created</b> but <b>NOT</b> running. Click <b>rerun</b> in ergatis to start running the pipeline.</font><br><br>";
							echo "<a href=\"{$res[' pipeline_url ']}\" target=\"_blank\">{$res[' pipeline_url ']}</a>";
						} else {
							echo "<br>";
							echo "<h3><font color=\"red\">Error creating the pipeline !! Contact system administrator.</font></h3>";
							echo "<br>";
							echo "Error : $pipeline_url<br>";
						}
					} else {
						echo "<br>";
						echo "<h3>Hit the browser's Back button and resolve the following errors to proceed:</h3>";
						echo "<ul>";
						foreach ($formFieldsArr as $formField) {
							if ($formValuesArr[$formField]['error'] > 0) {
								echo "<li><font color=\"red\">ERROR !! {$formValuesArr[$formField]['msg']}</font></li>";
							}
						}
						echo "</ul>";
					}
				?>
<?php
include_once('footer.php');
/*			</div>
		</div>
	</body>
</html>
*/
?>
