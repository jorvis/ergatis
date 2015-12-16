<?php
	$formFieldsArr = Array("rgene_algo" => "", "output_dir" => "Output directory", "log_file" => "Log file");
	$args = "";
	$local_dir = "/local/scratch/sadkins_devel/prok_pipeline_runs/pipeline_dir";
	$errFlag = 0;
	if (isset($_POST['bsubmit'])) {
		if (isset($_POST['rgene_algo'])) {
			$formValuesArr['rgene_algo']['default'] = trim($_POST['rgene_algo']);
		} else {
			$formValuesArr['rgene_algo']['default'] = "";
		}
		$formValuesArr['rgene_algo']['error'] = 0;
		$formValuesArr['rgene_algo']['msg'] = "";
		if(isset($_POST['selections'])) {
			$options = $_POST['selections'];
			foreach ($options as $opt) {
				$formValuesArr[$opt]['default'] = 1;
				switch($opt) {
					case "cpseudomolecule"	:	$args .= "--pseudomolecule 1 ";
									break;
					case "cannotate"	:	$args .= "--annotation 1 ";
									break;
					case "cipd"	:	$args .= "--ipd 1 ";
									break;
					case "cmulti"	:	$args .= "--multifasta 1 ";
									break;
				}
			}
			if (isset($formValuesArr['crna']['default'])) {
			    $args .= "--rna_prediction 1 ";
			} else {
			 	$args .= "--rna_prediction 0 ";
			}
			if (isset($formValuesArr['cgene']['default'])) {
				if ($formValuesArr['rgene_algo']['default'] == "glimmer") {
					$args .= "--gene_prediction glimmer ";
				} else {
					$args .= "--gene_prediction prodigal ";
				}
			} else {
				$args .= "--gene_prediction none ";
			}
			if (isset($formValuesArr['cload_db']['default']) && !(isset($formValuesArr['cannotate']['default']))) {
				$errFlag++;
				$formValuesArr['cload_db']['error'] = 1;
				$formValuesArr['cload_db']['msg'] = "Without annotating the genome, loading into the database is not possible";
			} elseif (isset($formValuesArr['cload_db']['default']) && isset($formValuesArr['cannotate']['default'])) {
				$args .= "--load 1 ";
			}
			if (isset($formValuesArr['cgenecalls']['default'])) {
				$args .= "--genecalls 1 ";
			}
		}
		error_reporting(0);
		$dir = create_pipeline_dir($local_dir);
		$formValuesArr['output_dir']['default'] = $dir;
		if(!dir) {
			$errFlag++;
			$formValuesArr['output_dir']['error'] = 1;
			$formValuesArr['output_dir']['msg'] = "Default output directory could not be created";
		} else {
			$args .= "--output_directory $dir ";
		}
		error_reporting(-1);
		if (isset($dir)) {
			$formValuesArr['log_file']['default'] = $dir."/create_prok_pipeline_config.log";
			$args .= "--log $dir"."/create_prok_pipeline_config.log ";
		} else {
			$errFlag++;
			$formValuesArr['log_file']['error'] = 1;
			$formValuesArr['log_file']['msg'] = "Could not create log file";
		}

#		echo "$args<br>";
	}

	function create_pipeline_dir ($local_dir) {
		$dir_num = mt_rand(1, 999999);
		$temp_num = str_pad($dir_num, 6, "0", STR_PAD_LEFT); 
		$dir = $local_dir.$temp_num;
		if (!file_exists($dir)) {
			if(!(mkdir($dir, 0777, true))) {
				return(0);
			} else {
				return($dir);
			}
		} else {
			create_pipeline_dir();
		}
	}
?>
