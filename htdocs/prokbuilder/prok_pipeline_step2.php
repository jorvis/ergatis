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
							 'cipd' => 'Connect to IPD');
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
				<h3>Configure the following global parameters for the prokaryotic pipeline : </h3>
				<p>For an explanation of what goes in the form fields, <a href='./help.php#form' target='_blank'>please click here</a></p>
				<form id="prok_pipeline_step2_form" name="prok_pipeline_step2_form" method="post" action="prok_pipeline_complete.php">
					<?php
						$pipeline_config = `find {$formValuesArr['output_dir']['default']} -name "*.config" -type f`;
						$pipeline_layout = `find {$formValuesArr['output_dir']['default']} -name "*.layout" -type f`;
						if (isset($pipeline_config) && isset($pipeline_layout)) {
							$params = file(trim($pipeline_config), FILE_IGNORE_NEW_LINES);
							$line_num = 0;
							$host = "manatee-db";
							echo "<table class='config_table' cellspacing='10'>";
							while (!empty($params[$line_num])) { 
								if (!($params[$line_num] == "[global]")) {
									list($header,$val) = preg_split("/=/",$params[$line_num]);
									$temp = preg_replace('/\$\;/','',$header);
									$formatted_head = preg_replace('/\_/',' ',$temp);
									echo "<tr><td>$formatted_head</td>";
									if ($formatted_head == "HOST") {
										echo "<td>$host</td>";
										echo "<td><input type='hidden' name=\"$header\" value=\"$host\"></td></tr>";
									} else if ($formatted_head == "GENOME" && isset($formValuesArr['cipd']['default']) && !(isset($formValuesArr['cipd']['error']))) {
										# If IPD is enabled, Genome will be filled in automatically rather than input by the user
										echo "<td><p id='g_val'></p></td>";
										echo "<td><input type='hidden' id='genome_val' name=\"$header\" value=\"$val\"></td></tr>";
									} else if ($formatted_head == "IPD STUDY STAGE") {
										# Automatically is filled in when study stage is obtained or found to not exist
										echo "<td><p id='ss_val'></p></td>";
										#echo "<td>(Automatic annotation stage ID or -1 if none exist)</td>";
										echo "<td><input type='hidden' id='study_stage_val' name=\"$header\" value=\"$val\"></td></tr>";
									} else {
										# Default behavior is the let user fill in the field
										echo "<td><input type='text' name=\"$header\" value=\"$val\"></td></tr>";
									}
									echo "<input type='hidden' name=\"paramsArr[]\" value=\"$params[$line_num]\">";
								}
								$line_num++;
							}
							echo "<tr><td>REPOSITORY ROOT</td>";
							echo "<td><input type='text' name='rep_root'></td></tr>";
							echo "</table>";

							if (isset($formValuesArr['cipd']['default'])) {
								echo "<table class='dropdown_table' cellspacing='10'>";
								echo "<td><input type='hidden' name='project_id' id='project_id' value=''></td>";
								echo "<td><input type='hidden' name='study_id' id='study_id' value=''></td>";
								echo "<td><input type='hidden' name='study_stage_id' id='study_stage_id' value=''></td>";
								$id = array();		
								if ($_SERVER['SERVER_NAME'] == 'localhost' || $_SERVER['SERVER_NAME'] == 'sadkins-lx.igs.umaryland.edu') {
									$client = 'devel';
								} else {
									$client = 'production';
								}	
								echo "<td><input type='hidden' name='client' id='client' value=\"$client\"></td>";
								$command = "/usr/bin/perl ./perl/IPD_to_list.pl $client";
								exec($command, $id);	#get IPD project names for a drop-down list
								echo "<tr><td><select id='project'>";
								echo "<option value=0>-- Please select a project --</option>";
								IPD_DropDownList($id);
								echo "</select></td></tr>";
								echo "<tr><td><select id='sample'>";
								echo "</select></td></tr>";
								echo "<tr><td><div id='ssDiv'><input type='checkbox' name='ss_create' id='ss_create' value='c_ss'>Create a new automatic annotation study stage</div></td></tr>";
								echo "</table>";
							}							

							echo "<br><br>";
							echo "<input type='hidden' name='config_file' value=\"$pipeline_config\">";
							echo "<input type='hidden' name='layout_file' value=\"$pipeline_layout\">";
							echo "<input type='reset' name='breset' value='Reset'>&nbsp&nbsp&nbsp";
							if ($errFlag) {
								echo "<input type='submit' name='bsubmit' value='Submit' disabled='disabled'>";
							} else {
								echo "<input type='submit' name='bsubmit' value='Submit' onclick=\"this.value=\"Processing ...\"\">";
							} 
						} else {
							echo "<font color='red'>Error creating pipeline.config and pipeline.layout files</font><br>";
						}

						function IPD_DropDownList($list = array()) {
							$arr = array();
							foreach ($list as $index => $key) {
								preg_match('/\s+\(\w{5}\)\s+\-\-\-(\d+)\-\-\-/', $key, $match);
								$value = $match[1];
								$pattern = '/\-\-\-\d+\-\-\-/';
								$key2 = preg_replace($pattern, '', $key);	#remove the ID from the key (---1234---)
								preg_match('/\s+(\(\w{5}\))/', $key2, $m);
								$abbr = $m[0];
								$abbr_ptrn = '/\(\w{5}\)/';
								$full = preg_replace($abbr_ptrn, '', $key2);    #remove the LIMS ID from the key (PVCHO)
								$arr[$value] = (preg_match('/(12345)/', $key2) ) ? $full : $abbr. "\t" . $full;	# set key in form to sort by LIMS
								# if no LIMS was created for project, a dummy value (12345) that aided in regex matching will be removed

							}
							asort($arr);	# sort by name while keeping index association
							foreach ($arr as $val => $name) {
								echo "<option value=$val>$name</option>";
							}
							return;
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
