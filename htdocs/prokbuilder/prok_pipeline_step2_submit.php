<?php
	$formFieldsArr = Array("config_file", "layout_file", "rep_root");
	$errFlag = 0;
	$ipdFlag = 0;
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
		if (strlen($formValuesArr['config_file']['default']) > 0) {
			$formValuesArr['config_file']['values'] = file($formValuesArr['config_file']['default'], FILE_IGNORE_NEW_LINES);
		} else {
			$errFlag++;
			$formValuesArr['config_file']['error'] = $errFlag;
			$formValuesArr['config_file']['msg'] = "pipeline.config not available";
		}
		if (strlen($formValuesArr['layout_file']['default']) == 0) {
			$errFlag++;
			$formValuesArr['layout_file']['error'] = $errFlag;
			$formValuesArr['layout_file']['msg'] = "pipeline.layout not available";
		} else {
			chmod($formValuesArr['layout_file']['default'], 0777);
		}
		if (strlen($formValuesArr['rep_root']['default']) == 0) {
			$errFlag++;
			$formValuesArr['rep_root']['error'] = $errFlag;
			$formValuesArr['rep_root']['msg'] = "Repository root is a mandatory parameter. Could not proceed without entering it"; 
		} 
		if(isset($_POST['paramsArr'])) {
			$paramsArr = $_POST['paramsArr'];
			foreach ($paramsArr as $line) {
				list($header,$val) = preg_split("/=/",$line);
				if (isset($_POST[$header])) {
					if ($header == '$;IPD_STUDY_STAGE$;') {
						if ($_POST[$header] < 0 && isset($_POST['ss_create'])) {
							$ss_id = create_study_stage();
							$_POST[$header] = $ss_id;
						}
					}
					$paramVal[$header] = trim($_POST[$header]);
				}
			}
		}

		if ($errFlag == 0) {
			$fh = fopen($formValuesArr['config_file']['default'], 'w') or $errflag++ ;
			if ($errFlag == 0) {
				foreach ($formValuesArr['config_file']['values'] as $data) {
					$count = count(preg_split('/=/', $data));
					if ($count == 2) {	// some data lines like [global] only have 1 element in array
						list($k,$v) = preg_split("/=/",$data);
						if (isset($paramVal[$k])) {
							fwrite($fh, "$k=$paramVal[$k]\n");
						} else {
							fwrite($fh, "$data\n");
						}
					} else {
						fwrite($fh, "$data\n");
					}
				}
				if (isset($_POST['project_id']) > 0) {
					if (isset($_POST['study_id']) == 0) {
						$formValuesArr['config_file']['error'] = $errFlag;
						$formValuesArr['config_file']['msg'] = "A project was selected from IPD but a sample was not.  Will not be written to pipeline.config file";
					}
				} else {
					$formValuesArr['config_file']['error'] = $errFlag;
					$formValuesArr['config_file']['msg'] = "No project ID was selected.  No changes made to the pipeline.config file";
				}
				chmod($formValuesArr['config_file']['default'], 0777);
				fclose($fh);
			} else {
				$formValuesArr['config_file']['error'] = $errFlag;
				$formValuesArr['config_file']['msg'] = "Could not open pipeline.config file for writing. No changes made to the pipeline.config file";
			}	
		} else {
			$formValuesArr['config_file']['error'] = $errFlag;
			$formValuesArr['config_file']['msg'] = "No changes made to the pipeline.config file";
		}
	}

	function create_study_stage() {
		$study_id = $_POST['study_id'];
		$client = $_POST['client'];
		$command = "/usr/bin/perl ./perl/create_study_stage.pl $client $study_id";
		$ss_id = shell_exec($command);	#create a study stage and return the ID
	 	return $ss_id;
	}
?>
