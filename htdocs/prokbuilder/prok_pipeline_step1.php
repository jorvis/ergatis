<?php
$page_title = 'Pipeline Builder - Step 1';
$extra_head_tags = '<script type="text/javascript" src="./js/ProkPipe1.js"></script>';
include_once('header.php');
/*
<html>
	<head>
	</head>
	<body>
		<div id='page_container'>
			
			<div id='content_container'>
*/  ?>
				<form id='prok_pipeline_step1_form' name='prok_pipeline_step1_form' method='post' action='prok_pipeline_step2.php'>
					<br>
					<h2>STEP 1 : Configure the prokaryotic pipeline<sup><a href='./help.php#form' target='_blank'>?</a></sup></h2>
					<h3>Select the components to be included in the pipeline from the following list : (Default selections are shown)</h3>
					<table cellspacing="10">
						<tr>
							<td><input type="checkbox" name="selections[]" id="selections" value="cpseudomolecule" onclick="CheckSubmit()"></td> 
							<td>Create pseudomolecule</td>
							<td></td>
						</tr>
						<tr>
							<td><input type="checkbox" name="selections[]" id="selections" value="cgenecalls" onclick="CheckSubmit(); SelectedGenecalls()"></td>
							<td>Annotate on established genecalls</td>
							<td></td>
						</tr>
						<tr>
							<td><input type="checkbox" name="selections[]" id="selections" value="crna" checked onclick="CheckSubmit()"></td>
							<td>Predict RNA [tRNA and rRNA]</td>
							<td></td>
						</tr>
						<tr>
							<td valign="top"><input type="checkbox" name="selections[]" id="selections" value="cgene" checked onclick="CheckSubmit()"></td>
							<td>
								Gene Prediction Algorithm. Select one
								<table cellspacing="10">
									<tr>
										<td><input type="radio" name="rgene_algo" value="glimmer" checked></td>
										<td>Glimmer</td>
									</tr>
									<tr>
										<td><input type="radio" name="rgene_algo" value="prodigal"></td>
										<td>Prodigal</td>
									</tr>
								</table>
							</td>
							<td></td>
						</tr>
						<tr>
							<td><input type="checkbox" name="selections[]" id="selections" value="cannotate" checked onclick="CheckSubmit()"></td>
							<td>Annotate</td>
							<td></td>
						</tr>
						<tr>
							<td><input type="checkbox" name="selections[]" id="selections" value="cload_db" onclick="CheckSubmit()"></td>
							<td>Load into Database</td>
							<td></td>
						</tr>
						<tr>
							<td><input type="checkbox" name="selections[]" id="selections" value="cmulti" checked disabled onclick="CheckSubmit()"></tdd>
							<td>Multisequence Output (reduces number of files created)</td>
							<td></td>
						</tr>
						<tr>
							<td><input type="checkbox" name="selections[]" id="selections" value="cipd" onclick="CheckSubmit()"></td>
							<td>Connect to IPD Study Stage</td>
							<td></td>
						</tr>
						<tr>
							<td><input type="checkbox" name="selections[]" id="selections" value="cbatch" onclick="CheckSubmit();SendToBatch()"></td>
							<td>Create multiple pipelines</td>
							<td></td>
						</tr>								
					</table>
					<br>
						<tr>
							<td><input type="button" name="bselect_all" onclick="SelectUnselectAll()" value="Select All"></td>
							<td><input type="reset" name="breset" onclick="ChangeLabel()" value="Reset"></td>
							<td><input type="submit" name="bsubmit" value="Submit"></td>
						</tr>
					</table>
				</form>
				
<?php
include_once('footer.php');
/*			</div>
		</div>
	</body>
</html>
*/
?>
