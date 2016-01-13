
<?php
$page_title = 'Pipeline Builder - Help';
$extra_head_tags = '<meta http-equiv="Content-Language" content="en-us">' . "\n";
$extra_head_tags .= '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">' . "\n";
$extra_head_tags .= '<link rel="stylesheet" type="text/css" href="./css/documentation.css">';
include_once('header.php');
/*
<html>
	<head>
	</head>
	<body>
		<div id='page_container'>
			
			<div id='content_container'>
*/  ?>
				<h2><a href="help.php">Documentation</a></h2>
				<p>
					This is a simple interface build to create and configure prokaryotic annotation pipeline in 
					<a href="http://ergatis.igs.umaryland.edu" target="_blank"><strong>Ergatis</strong></a>. It involves the following <a href="images/documentation/Prok_Pipeline_Flowchart.pdf" target="_blank">steps</a> :
				</p>
				<ol>
					<li>Selecting from a list of available sections to be included in the pipeline. The sections are :
						<ul>
							<li><a href='#create_pseudomolecule'>Create pseudomolecule</a></li>
							<li><a href='#predict_rna'>Predict RNA</a></li>
							<li><a href='#gene'>Predict protein coding genes</a></li>
							<li><a href='#annotate'>Annotate the predicted genes</a></li>
							<li><a href='#load'>Load the genome with annotation in the database</a></li>
							<li><a href='#multi'>Multisequence output</a></li>
							<li><a href='#ipd'>Associate the pipeline with a sample from the IGS Projects Database (IPD)</a></li>
						</ul>
					    An output directory will be created and listed on the next page after Submit is clicked.  This output directory will contain the requisite pipeline.config and pipeline.layout files needed to create the pipeline in addition to a log file.
					</li>
					<li>Configuring the pipeline.config file based on the selections made in step 1.<a href='#form'> Here is a list of form fields that you may encounter.</a></li>
					<li>The last page will display a link to the newly created Ergatis pipeline. Clicking on the link will take you directly to the Ergatis interface to run the and monitor the pipeline.</li>
				</ol>
				<p>
					If there are any issues with the Pipeline Builder interface, or if you have any questions, comments, etc. then feel free to send an e-mail to <strong>sadkins [at] som.umaryland.edu</strong>
				</p>
				<a name='create_pseudomolecule'></a>
				<h2>Create pseudomolecule</h2>
				<p>
					<dl>					
						<dt>Required Inputs</dt>
						<dd>Contigs or sequence file in FASTA format</dd>
						<dt>Optional Inputs</dt>
						<dd>Reference accession ids file &#45 This is a tab-delimited file with the following 2 columns :
							<ol>
								<li>Reference accession id &#45 GenBank id of the genome to be used as reference for the input contigs or sequence file.</li>
								<li>Group Number &#45 This is an arbitrary sequential number assigned to the reference accession ids in cases where there are more 
								    than one reference genomes, separate for chromosome and plasmids. All the chromosome reference genomes ids should be assigned 
								    the same group and all the reference plasmids ids another group.
								</li>
							</ol>
						</dd>
					</dl>
					During the pseudomolecule creation process, the reference genome sequences [chromosome references] in the first group are concatenated 
					together and aligned with the input contigs. The left over contigs which do not align in this first pass are then aligned to the 
					concatenated reference sequences [plasmid references] in the next group and so on till all the groups are completed or there are no more 
					contigs left for alignment whichever happens first. The last set of unaligned contigs are concatenated together to form a pseudomolecule.  
					If no reference accession ids file is passed, the contigs will be arranged in decreasing order of their lengths to form a single 
					pseudomolecule.
				</p>
				<a name='predict_rna'></a>
				<h2>Predict RNA</h2>
				<p>
					Including this component will predict <strong>tRNAs</strong> and <strong>rRNAs</strong> in the genomic sequence.
				</p>  
				<a name='gene'></a>
				<h2>Predict protein coding genes</h2>
				<p>
					One of the two algorithms/programs could be selected to predict protein coding genes in the genomic sequence. The algorithms are : 
					<dl>
						<dt>Glimmer</dt>
						<dd>This is the default algorithm used for gene prediction in the pipeline.</dd>
						<dt>Prodigal</dt>
					</dl>
					If you have a tab delimited or gff file containing gene calls, convert it into a bsml file using 
					<a href="/home/sagrawal/ergatis_dev/prodigal/genecalls2bsml.pl">genecalls2bsml.pl</a> and select annotation section to annotate on these gene calls.
				</p>
				<a name='annotate'></a>
				<h2>Annotate the predicted genes</h2>
				<p>
					Several softwares are used to perfom gene annotation like BLASTX, BER, HMMER, SIGNALP, LIPOP, TMHMM, etc.  
				</p>
				<a name='load'></a>
				<h2>Load the genome with annotation in the database</h2>
				<p>
					The annotated genome could be stored into a <strong>CHADO</strong> database. Without selecting the annotate section above this section would not work. 
					Only annotated genomes could be loaded into the database. The annotated genome could be viewed in a visualization tool called 
					<a href="http://manatee.igs.umaryland.edu"><strong>Manatee</strong></a>.
				</p>
				<a name='multi'></a>
				<h2>Change files to multisequence output</h2>
				<p>
					Many pipeline components, such as translate_sequence and bsml2fasta are set up where only one sequence per file is written, meaning a lot of smaller files can accumulate in the output repository per pipeline.  Enabling this option will reduce the file volume by writing 100 sequence entry outputs per file instead of 1. 
				</p>
				<a name='ipd'></a>
				<h2>Associate the pipeline with a sample from the IGS Projects Database (IPD)</h2>
				<p>
					This option will add a drop-down list of projects sorted by LIMS ID found from <a href="http://projects.igs.umaryland.edu" target="_blank"><strong>IPD</strong></a> on the step 2 page.  Selecting a project will reveal another box where you can search for any sample whose study is listed as currently active.  Selecting one of these samples will automatically fill in the GENOME section of the config file.  
					In addition the automatic annotation study stage ID for that study will be collected if it exists, and will be used to fill in the IPD STUDY STAGE section of the config file.  If there is no automatic annotation study stage for this study, then IPD STUDY STAGE will be listed as -1, and the checkbox that will allow you to create an automatic annotation study stage ID will be present.  If checked, the new study stage ID will be used in the config file.
					A component will be added to the end of the pipeline that will mark the automatic annotation study stage ID as complete.  If no automatic annotation study stage ID existed and the user does not wish to create one, then the component will finish without doing anything.
				</p>
				<a name='form'></a>
				<h2>Form fields used in step 2</h2>
					<dl>
						<dt>Abbreviation</dt>
							<dd>Used in id generation when generating sequences in the pipeline.  Locus prefixes make good abbreviations.  Try to avoid underscores or other random symbols if possible.</dd>
						<dt>Accession File</dt>
							<dd>A tab-delimited file with each line consisting of the Genbank accession ID and a designated group number.  Will be used to create pseudomolecules for each group number.</dd>
						<dt>DB</dt>
							<dd>The name of the Chado database that the analyses from this pipeline should be loaded into.</dd>
						<dt>Genome</dt>
							<dd>This is where you specify the species you wish to use.  If the 'Connect to IPD' option is selected, then this will be filled in automatically upon choosing a sample from an IPD project.</dd>
						<dt>Host</dt>
							<dd>Server where the Chado database is held.  Automatically set to 'manatee-db'.</dd>
						<dt>Input BSML List</dt>
							<dd>This is a BSML file typically based on the output from the promote_gene_prediction component in Ergatis.  This field will typically show up if the user selected for no gene prediction to be performed.</dd>
						<dt>Input FSA File</dt>
							<dd>Specify the location of your FASTA file here.</dd>
						<dt>Input Fasta</dt>
							<dd>See Input FSA File.</dd>
						<dt>Input FSA List</dt>
							<dd>If you have a list of FASTA-formatted files that you wish to use as input, then use this field instead of Input FSA File.</dd>
						<dt>IPD Study Stage</dt>
							<dd>Automatic annotation study stage ID acquired from IPD.  Will be automatically filled in after a sample is chosen.  If no automatic annotation study stage exists for the corresponding study associated with this sample, then it will be labeled as '-1' and a checkbox will show up below asking if a new study stage should be created for this sample.</dd>
						<dt>Pass</dt>
							<dd>Password for the provided username for the database.</dd>
						<dt>Repository Root</dt>
							<dd>The directory where the output and workflow process files should be placed.  The directory structure for your repository root project area must have already been created before a new pipeline is created in the area.  If a new project area needs to be created, <a href="http://ergatis.igs.umaryland.edu/cgi/create_project.cgi" target="_blank">please click here</a>.</dd>
						<dt>TaxID</dt>
							<dd>Taxon ID of the genome provided.  Most likely, the ID you will want to provide is the taxonomy ID taken from NCBI's <a href="http://www.ncbi.nlm.nih.gov/taxonomy" target="_blank">Taxonomy Browser</a>.</dd>
						<dt>Trans Table</dt>
							<dd>The ID of the codon translation table you want to use.  The default (which is most likely what you want) is '11' since this is the code for prokaryotic proteins.</dd>
						<dt>User</dt>
							<dd>Username for the Chado database specified.</dd>
					</dl>
				<a name='batch'></a>
				<h2>Metadata File Fields</h2>
				<p>A tab-separated or comma-separated file containing the following columns:</p>
				<ol>
					<li>Repository Root path</li>
					<li>Abbreviation (name that will be typically tagged onto output files and headers)</li>
					<li>Genome (Two space-separated names)</li>
					<li>Taxon ID</li>
					<li>Input Fasta File</li>
					<li>Input Fasta List (if you want to pass that in instead of a file)</li>
					<li>Input BSML List</li>
					<li>Accession File Path</li>
					<li>Chado DB name</li>
					<li>Chado Host server (defaults to manatee-db)</li>
					<li>Chado server username</li>
					<li>Chado server password</li>
					<li>IPD Study Stage (put Study Stage ID)</li>
 				</ol>
<?php
include_once('footer.php');
/*			</div>
		</div>
	</body>
</html>
*/
?>
