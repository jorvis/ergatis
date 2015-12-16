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
					$formFieldsArr = Array("config_file", "layout_file", "metadata");
					if (isset($_POST['bsubmit'])) {
						foreach ($formFieldsArr as $formField) {
							if (isset($_POST[$formField])) {
								$formValuesArr[$formField]['default'] = trim($_POST[$formField]);
							} else {
								$formValuesArr[$formField]['default'] = "";
							}
							$formValuesArr[$formField]['error'] = 0;
							$formValuesArr[$formField]['msg'] = "";
						}
						if (strlen($formValuesArr['config_file']['default']) == 0) {
							$errFlag++;
							$formValuesArr['config_file']['error'] = $errFlag;
							$formValuesArr['config_file']['msg'] = "pipeline.config not available";
						} else {
							chmod($formValuesArr['config_file']['default'], 0777);
						}
						if (strlen($formValuesArr['layout_file']['default']) == 0) {
							$errFlag++;
							$formValuesArr['layout_file']['error'] = $errFlag;
							$formValuesArr['layout_file']['msg'] = "pipeline.layout not available";
						} else {
							chmod($formValuesArr['layout_file']['default'], 0777);
						}
						if (strlen($formValuesArr['metadata']['default']) == 0) {
							$errFlag++;
							$formValuesArr['metadata']['error'] = $errFlag;
							$formValuesArr['metadata']['msg'] = "Metadata file is a mandatory parameter. Could not proceed without entering it";
						}
					}


					if ($errFlag == 0) {
                        $stdout = `/usr/bin/perl ./perl/create_batch_prokpipes.pl --config_file {$formValuesArr['config_file']['default']} --input_file {$formValuesArr['metadata']['default']}`;
						if (preg_match("/=/", $stdout, $matches)){
							echo "<br>";
							echo "<h2>STEP FINAL : Pipelines created successfully<sup><a href=\"./help.php\">?</a></sup></h2>";
							echo "<h3>Following files have been created :</h3>";
							echo "<ul>";
							echo "<li>Pipeline Configuration File Template : {$formValuesArr['config_file']['default']}</li>";
							echo "<li>Pipeline Layout File Template: {$formValuesArr['layout_file']['default']}</li>";
							echo "</ul>";
							echo "<br>";
							echo "<h3>The following repositories have new pipelines waiting to start running:</h3>";
							$temp = preg_split("/\n(\r)?/", $stdout);	# could actually use explode("\n", $stdout); as well
                            $tmp = array_pop($temp);	# Since last line also has "\n", it adds a blank line
                            # Go through each entry and link the repo root to the page
							foreach ($temp as $t) {
								list($repo, $count) = preg_split("/=/",$t);
								if ($count == 1) {
									echo "<p>$count new pipeline created in <a href=\"http://ergatis.igs.umaryland.edu/cgi/pipeline_list.cgi?repository_root=$repo\" target='_blank'>$repo</a></p>";
								} else {
									echo "<p>$count new pipelines created in <a href=\"http://ergatis.igs.umaryland.edu/cgi/pipeline_list.cgi?repository_root=$repo\" target='_blank'>$repo</a></p>";
								}
							}
							echo "<p>NOTE:  The pipelines do not automatically start running.  To start a pipeline, click any repository root link above, click the numerical ID of an <b>incomplete</b> pipeline, and hit <b>Rerun</b></p>";
						} else {
							echo "<h3><font color='red'>ERROR !! No pipelines were created from metadata file</font></h3>";
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
