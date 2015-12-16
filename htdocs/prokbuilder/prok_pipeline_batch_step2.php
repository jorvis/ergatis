<?php
$page_title = 'Pipeline Builder - Step 2';
$extra_head_tags .= '<link rel="stylesheet" type="text/css" href="./css/jquery-ui.css " media="screen"></link>' . "\n";
$extra_head_tags .= '<script type="text/javascript" src="./js/jquery.js"></script>' . "\n";
$extra_head_tags .= '<script type="text/javascript" src="./js/jquery-ui-1.10.1.custom.js"></script>' . "\n";
$extra_head_tags .= '<script type="text/javascript" src="./js/DropDownList.js"></script>';
$extra_head_tags .= '<style>
						.ui-autocomplete {
							max-height: 100px;
							overflow-y: auto;
							overflow-x: hidden;
							padding-right: 20px;
						}
					</style>';
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
					include_once('prok_pipeline_step1_submit.php');
					if ($errFlag == 0) {
						$output = `/usr/bin/perl ./perl/create_prok_pipeline_config.pl $args --template_directory /local/projects/ergatis/package-latest/global_pipeline_templates`;
					}
					$content = array('cpseudomolecule' => 'Create pseudomolecule',
							 'cgenecalls' => 'Genecalls Pipeline',
							 'crna' => 'Predict RNA',
							 'cgene' => 'Gene Prediction Algorithm : '.$formValuesArr['rgene_algo']['default'],
							 'cannotate' => 'Annotate',
							 'cload_db' => 'Load into database',
							 'cmulti' => 'Multisequence output',
							 'cipd' => 'Connect to IPD',
							 'cbatch' => 'Batch pipeline creation mode');
				?>
				<br>
				<h2>STEP 2 : Configure the prokaryotic pipeline components<sup><a href="./help.php" target='_blank'>?</a></sup></h2>
				<h3>You selected the following major components to be included in the pipeline</h3>
				<ul>
					<?php
						foreach ($content as $cnt => $desc) {
							if (isset($formValuesArr[$cnt]['default']) && !(isset($formValuesArr[$cnt]['error']))) {
								echo "<li>$desc</li>";
							} elseif (isset($formValuesArr[$cnt]['error'])) {
								echo "<li>$desc : <font color=\"red\"><b>ERROR !! {$formValuesArr[$cnt]['msg']}</b></font></li>";
							}
						}
					?>
				</ul>
				<h3>Following directories and files have been created</h3>
				<ul>
					<?php
						foreach ($formFieldsArr as $formField => $val) {
							if ($formField != 'rgene_algo') {
								if (isset($formValuesArr[$formField]['error'])) {
									echo "<li>$val : <font color=\"red\"><b>ERROR !! {$formValuesArr[$formField]['msg']}</b></font></li>";
								} else {
									echo "<li>$val : {$formValuesArr[$formField]['default']}</li>";
								}
							}
						}
					?>
				</ul>
				<h3>Provide the path to a metadata file of relevant fields : </h3>
				<p>The fields in the metadata file can be either comma or tab delimited.  For an explanation of what order the fields are, and what each field represents, <a href='./help.php#batch' target='_blank'>please click here</a></p>
				<form name="prok_pipeline_batch_step2_form" method="post" action="prok_pipeline_batch_complete.php">
					<?php
						$pipeline_config = `find {$formValuesArr['output_dir']['default']} -name "*.config" -type f`;
						$pipeline_layout = `find {$formValuesArr['output_dir']['default']} -name "*.layout" -type f`;
						if (isset($pipeline_config) && isset($pipeline_layout)) {
							echo "<table class='config_table' cellspacing='10'>";
							echo "<tr><td>METADATA FILE PATH</td>";
							echo "<td><input type='text' name='metadata'></td></tr>";
							echo "</table>";
							echo "<br><br>";
							echo "<input type='hidden' name='config_file' value=\"$pipeline_config\">";
							echo "<input type='hidden' name='layout_file' value=\"$pipeline_layout\">";
							echo "<input type='reset' name='breset' value='Reset'>&nbsp&nbsp&nbsp";
							if ($errFlag) {
								echo "<input type='submit' name='bsubmit' value='Submit' disabled='disabled'>";
							} else {
								echo "<input type='submit' name='bsubmit' value='Submit' onclick=\"this.value='Processing ...'\">";
							}
						} else {
							echo "<font color='red'>Error creating pipeline.config and pipeline.layout files</font><br>";
						}
					?>
				</form>
<?php
include_once('footer.php');
/*			</div>
		</div>
	</body>
</html>
*/
?>
