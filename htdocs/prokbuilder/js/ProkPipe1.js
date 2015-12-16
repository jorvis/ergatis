
// Operates the Select All button
function SelectUnselectAll() {		
    var FormName = document.forms["prok_pipeline_step1_form"];
    var Field = FormName.selections
	if (FormName.elements['bselect_all'].value == 'Deselect All') {
        for (i = 0; i < Field.length; i++) {
        	Field[i].checked = false;
        	Field[i].disabled = false;
		}
		FormName.elements['bsubmit'].disabled = true;
		FormName.elements['bselect_all'].value = "Select All";
	} else {
		for (i = 0; i < Field.length; i++) {
			if (Field[i].value !== "cgenecalls") {
				Field[i].checked = true;
				Field[i].disabled = false;
			}
		}
		FormName.elements['bsubmit'].disabled = false;
		FormName.elements['bselect_all'].value = "Deselect All";
	}
}

// Hitting the Reset button will toggle the Select All button to select all checkboxes		
function ChangeLabel() {		
	var FormName = document.forms["prok_pipeline_step1_form"];
	var Field = FormName.selections;	
	FormName.elements['bsubmit'].disabled = false;
    for (i = 0; i < Field.length; i++) {
       	Field[i].disabled = false;
	}	
	if (FormName.elements['bselect_all'].value == 'Deselect All') {
		FormName.elements['bselect_all'].value = "Select All";
	} 
}

// Keeps track of the status of each checkbox
function CheckSubmit() {		
	var FormName = document.forms["prok_pipeline_step1_form"];
	var Field = FormName.selections;
	var box_checked = false;
	for (i = 0; i < Field.length; i++) {
		if (Field[i].checked) {
			box_checked = true;
		}
	}
	if(box_checked == true){
		FormName.elements['bsubmit'].disabled = false;
	} else {
		FormName.elements['bsubmit'].disabled = true;
	}
}	

// If batch pipeline creation is set, send to a different Step 2 form
function SendToBatch() {	
	var FormName = document.forms["prok_pipeline_step1_form"];
	var Field = FormName.selections;
	var box_checked = false;
	for (i = 0; i < Field.length; i++) {
		if (Field[i].value == 'cbatch' && Field[i].checked) {
			box_checked = true;
		}
	}	
	
	if (box_checked == true) {
		FormName.action = "prok_pipeline_batch_step2.php";
	} else {
		FormName.action = "prok_pipeline_step2.php";
	}
	return true;
}

// Uncheck and disable certain buttons if genecalls option is selected
function SelectedGenecalls() {	
	var FormName = document.forms["prok_pipeline_step1_form"];
	var Field = FormName.selections;
	var genecalls;
	var pseudo;
	var rna;
	var gene;
	var change = ["cpseudomolecule", "crna", "cgene"];
	
	// Get genecalls field and other fields
	for (i=0; i<Field.length; i++) {
		switch (Field[i].value) {
			case "cgenecalls":
				genecalls = Field[i];
				break;
			case "cpseudomolecule":
				pseudo = Field[i];
				break;
			case "crna":
				rna = Field[i];
				break;
			case "cgene":
				gene = Field[i];
				break;
		}
	}
	// If Genecalls is checked, then disable certain options.  If unchecked then enable
	if (genecalls.checked == true) {
		pseudo.checked = false;
		pseudo.disabled = true;
		rna.checked = false;
		rna.disabled = true;
		gene.checked = false;
		gene.disabled = true;
	} else {
		pseudo.disabled = false;
		rna.disabled = false;
		gene.disabled = false;	
	}
}
