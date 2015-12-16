$(document).ready(function (){
	$('#sample').hide();	//hide the sample and study drop-down-list at startup
	$('#ssDiv').hide();	//hide the study stage checkbox

	$('#project_id').val(0);
	$('#study_id').val(0);
	$('#study_stage_id').val(0);
	
	if ($("#project").length) {
		$("input[name='bsubmit']").prop("disabled", true);
	} else {
		$("input[name='bsubmit']").prop("disabled", false);	//Make sure Submit button is enabled initially in case IPD is not needed
	}

	$('#project').change(function () {
		$('#project_id').val($('#project').val());	//setting project_id value
		$('#study_id').val(0);	
		reset_genome_ss();
		
		if ( $('#sample').combobox() != "undefined") {	// to avoid errors when initially loading page (no sample list yet)
			$('#sample').combobox("destroy").hide();	
		}
		
		if ( $('#project_id').val() > 0 ) {	// Only load sample list if an actual project was chosen
			//alert("Please be patient while the sample list is loading for this project...");
			$(".dropdown_table").append("<tr><td><strong>Please be patient while the sample list is loading for this project...</strong></td></tr>");
			$.ajax ({
				url: './choose_sample.php',
				data: { 'proj_id' : $('#project_id').val() ,
					    'client' : $('#client').val()},
				cache: false,
				dataType: "text"
			}).done( function (data) { 
				$("strong").remove(":contains('Please be patient')");
				array_to_sample(data);
			}).fail( function(http_r, opts, thrown) {
				alert("status=" + opts + " "  + http_r.status + ", error=" + thrown);
			});
		}
	} );

	$("input[name='breset']").click(function() {
		reset_genome_ss();
		$('#sample').combobox("destroy").hide();
	});
} );

// Taking JSON sample data for the selected project and using it to create a Select HTML element
function array_to_sample (result) {
	var json = jQuery.parseJSON(result);
	var sample = $('#sample');

	var options = sample.prop('options');	// store existing options in variable
	$('option', sample).remove();	//remove existing options to replace
	options[options.length] = new Option("-- Please type an active sample name --", 0); 
	if ($.isEmptyObject(json)){
		alert ("This project does not appear to have any samples with an 'active' status");
	} else {
		$.each(json, function(key, val) {
			options[options.length] = new Option (key, val);	// loop through new values and add to Sample dropdown list
		});
	}	
	sample.combobox({
		// The trigger response when a sample option is selected
		selected: function(event, ui) {
			//console.log("Selected function " + $("#sample").val());
			$('#study_id').val( $('#sample').val() );	// values used are study ID keys associated with the sample
			var sample_text = $('#sample option:selected').text();
			var formatted = sample_text.replace(/\s\(ID:\s+\d+\)/, '');
			$('#g_val').text(formatted);	//autofills text area
			$("#genome_val").val($('#g_val').text());
			//$("input[name='bsubmit']").prop("disabled", true);	//Disable Submit until StudyStageID has been acquired

			if ( $('#study_id').val() <= 0 ) {	// For times when the default "please select a sample" option is checked
				reset_genome_ss();				
			} else {
				$.ajax ({	//get_ss_id (since we now already have study ID)
					url: './get_study_stage.php',
					data: { 'study_id' : $('#study_id').val() ,
						    'client' : $('#client').val() },
					cache: false,
					dataType: "text"
				}).done(function (ss_id) { 
					$('#study_stage_id').val(ss_id);
					$('#ss_val').text(ss_id);
					$("#study_stage_val").val($('#ss_val').text() );		
					if ($('#study_stage_id').val() <= 0) {
						$('#ssDiv').show();	//show checkbox if ID has not been found
					} else {
						$('#ssDiv').hide({
							complete: function () {
								$('#ss_create').prop('checked', false);
							}
						});
					} 
					$("input[name='bsubmit']").prop("disabled", false);
				}).fail( function(http_r, opts, thrown) {
					alert("status=" + opts + " "  + http_r.status + ", error=" + thrown);
				});
			}
		}
	});
}

// Creates the combobox widget
(function ($) {
	$.widget( "ui.combobox", {
		_create: function () {
			var self = this;
			var select = this.element,
				theWidth = select.width(),	// Width of combobox will be that of the dropdown menu underneath
				selected = select.children(":selected"),	// Initialize 
				value = selected.val() ? selected.text() : "";
			select.hide();	// Hide dropdown element to replace with text input
			var wrapper = this.wrapper = $( "<span>" )
				.insertAfter( select );
				
			var input =$("<input style=\"width:" + theWidth + "px\">")
				.appendTo( wrapper )
				.val(value)
				.attr( "title", "")
				.autocomplete({
					delay: 0,
					minLength: 0,
					source: function (request, response) {
						var matcher = new RegExp( $.ui.autocomplete.escapeRegex(request.term), "i" );	// defines RegEx Object based on what is entered
						response(select.children("option").map(function () {	// Maps each option from the select box
							var text = $(this).text();
							var text_f = text.replace(/\s\(ID:\s+\d+\)/, '');
							// if value of option has a match or is a blank search, return the option 
							if (this.value && (!request.term || matcher.test(text) )) return {
								label: text_f.replace(	new RegExp(
                   							"(?![^&;]+;)(?!<[^<>]*)(" +
                  							$.ui.autocomplete.escapeRegex(request.term) +
                   							")(?![^<>]*>)(?![^&;]+;)", "gi" ),
                 							"<strong>$1</strong>" ),  // Makes typed text bold in the results
								value: text,
								option: this
							};
						}));
					},	// END source option
					select: function(event, ui) {
						ui.item.option.selected = true;	// The select item that is chosen in the Autocomplete menu
						self._trigger("selected", event, {
							item: ui.item.option	// Update select to be this option via trigger
						});
					},
					change: function(event, ui) {
						if ( !ui.item) {
							var matcher = new RegExp( "^" + $.ui.autocomplete.escapeRegex( $(this).val() ) + "$", "i"),
							valid = false;	// Initializes RegExp object based on a changed entered value
							select.children("option").each(function () {
								if(this.value.match(matcher) ){
									this.selected = valid = true;	// keep selected option's value if new RegExp still matches
									return false;
								}
							});
							if (!valid) {	// Select's value is null if there is no match to the selected option
								$(this).val("");
								select.val("");
								return false;
							}
						}
					}	// END change event
				})	// END autocomplete
				.addClass("ui-state-default ui-widget ui-widget-content ui-corner-left");

			//render autocomplete results with the bold label
   			 input.data( "ui-autocomplete" )._renderItem = function( ul, item ) {	
     				 return $( "<li>" )
        				.data( "item.autocomplete", item )
        				.append( "<a>" + item.label + "</a>" )
        				.appendTo( ul );
    			};

			/*	USING A SCROLLBAR FOR THE AUTOCOMPLETE FOR NOW
			
			//creates button for opening up all options from the select menu
    			$( "<button>" )	
    				.attr( "tabIndex", -1 )
    				.attr( "title", "Show All Items" )
    				.tooltip()
    				.appendTo( wrapper )
    				.button({
      					icons: {
        					primary: "ui-icon-triangle-1-s"
      					},
      					text: false
    				})
    				.removeClass( "ui-corner-all" )
    				.addClass( "ui-corner-right ui-button-icon" )
    				.click(function() {
    				      	console.log("a Input " + input.autocomplete("widget").is(":visible"));
      					if (input.autocomplete("widget").is(":visible")) {      // close if already visible
        					input.autocomplete("close");
        					return false;		//adding a false return will ensure the button does not redirect
      					}
      					console.log("b Input " + input.autocomplete("widget").is(":visible"));
      					input.autocomplete("search", "");	// pass empty string as value to search for, displaying all results
      					input.focus();     
      					return false;
    				});
    				
    			input.tooltip({
    				tooltipClass: "ui-state-highlight"
    			});
    			*/
		},	// END _create function
		_destroy: function() {
			this.wrapper.remove();
			this.element.show();
		}
	});
}) (jQuery);

// Resets text assigned to Genome and IPD Study Stage fields
function reset_genome_ss () {	
	$('#g_val').text('');
	$("#genome_val").val($('#g_val').text());
	$('#ss_val').text('');
	$("#study_stage_val").val($('#ss_val').text() );	
	$('#ssDiv').hide({
		complete: function() {
			$('#ss_create').prop('checked',false);	// uncheck the create_study_stage box if sample is changed
		}
	});	// hide the study stage checkbox if no sample has been selected
	$("input[name='bsubmit']").prop("disabled", true);
}

/*	TOOLTIP FUNCTIONALITY IS NOT CURRENTLY BEING USED
$(function() {	// Enables tooltips for HTML elements that have 'title' attributes present within
	$( document ).tooltip();
});
*/
