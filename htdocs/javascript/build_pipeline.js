// these colors correspond to values in build_pipeline.css
var color_invalid    = 'rgb(255,150,150)';
var color_incomplete = 'rgb(230,230,230)';
var color_error      = 'rgb(200,0,0)';

var next_set_num = 0;
var new_input_fields_row;
var no_input_message_row;
var repository_root;
var workflowdocs_dir;

// component config fetches and saves sometimes need to be throttled
var componentSaveQueue = new AjaxQueue( 3 );

// set in the ergatis.ini file
var build_directory;
var builder_animations = true;

var pipeline_root_node;
var pipeline_root_panel_node;

var component_being_configured;
var pipeline_insert_location;

var email_address;
var emailNotify = false;
var saveClicked = false;

jQuery(document).ready(function() {
    // make sure the form is reset here.
    document.pipeline.reset();

    // all pipelines start with at least these two nodes
    pipeline_root_node = new TreeNode('pipeline_root', 'root');
    pipeline_root_panel_node = new TreeNode('pipeline_root_panel', 'panel');
    pipeline_root_node.setNodeBelow(pipeline_root_panel_node);
    pipeline_root_panel_node.setNodeAbove(pipeline_root_node);

    if ( getObject('builder_animations').value == 0 ) {
        builder_animations = false;
    }

    new_input_fields_row = new fx.Combo('new_input_fields_row', {duration: 400});
    new_input_fields_row.hide();
    
    no_input_message_row = new fx.Combo('no_input_message', {duration: 400});
    
    repository_root = getObject('repository_root').value;
    build_directory = getObject('build_directory').value;
    workflowdocs_dir = getObject('workflowdocs_dir').value;
    
    // check if there is anything to autoload
    if ( ! isEmpty( getObject('autoload_template').value ) ) {
        // remember where this pipeline will be inserted
        pipeline_insert_location = 'pipeline_root';
        selectPipelineTemplate( getObject('autoload_template').value );
    }

    emailNotify = jQuery('#email_on_default').val();

    if (emailNotify == "1") {
        jQuery('#add_email').hide();
        jQuery('#email_input').show();
    }

    // If someone clicks the 'click to add' button for the email notification 
    // label we'll need to pop-up our input box
    jQuery('#add_email').click(function() {
        jQuery('#email_input').show();
        jQuery('#save_email').show();
        jQuery('#add_email').hide();
        jQuery('#add_email').css('visibility', 'hidden');
        emailNotify = true;
    });

    jQuery('#email_input').click(function() {
        email_address = jQuery(this).val();
        jQuery('#save_email').show();
    });

    jQuery('#email_input').keypress(function(e) {
        if (e.keyCode == 13) {
            jQuery('#save_email').mousedown();
            jQuery('#email_input').blur();
        }
    });

    jQuery('#email_input').mouseover(function() {
        inputDisplayEditable(this);
    });

    jQuery('#email_input').mouseout(function() {
        if (! jQuery(this).is(":focus")) {
            inputDisplayUneditable(this);
            jQuery('#save_email').hide();
        } 
    });

    jQuery('#save_email').mousedown(function() {
        saveClicked = true;
    });

    jQuery('#email_input').blur(function() {
        inputDisplayUneditable(this)

        if (! jQuery(this).val() ) {
            jQuery(this).hide();
            jQuery('#add_email').show();
            jQuery('#add_email').css('visibility', 'visible');
            emailNotify = false;
         
        } else{
            if (! saveClicked ) {
                jQuery(this).val(email_address);
            } else {
                emailNotify = true;
                saveClicked = false;
            }
        }

        jQuery('#save_email').hide();
    });

});

// we need to pass some linking information about each numbered component ID
//  that will get passed in the pipeline form submission
function addComponentInfo( component_id, name, token ) {
    var info_elm = document.createElement('input');
    info_elm.setAttribute('type', 'text');
    info_elm.setAttribute('name', component_id + '_name.token');
    info_elm.setAttribute('value', name + '.' + token);
    
    getObject('variables').appendChild( info_elm );
}

function addComponentStub( set_root ) {
    var component_id = 'c' + components.length;

    // this will be inserted above the panel clicked

    var node_below = set_root + '_panel';
    var node_above = getObject(set_root + '_panel_up').value;
    
    var component_html = Component.templateHtml( component_id );
    
    getObject(set_root + '_panel').insertAdjacentHTML('BeforeBegin', component_html);

    components[component_id] = new Component( component_id );
    components[component_id].setPosition( node_below, node_above );
    
    // make sure the component is in basic view
    components[component_id].clearConfig();
    
    // the copy icon should be disabled until the component is saved
    getObject(component_id + '_copy').style.display = 'none';
    
    // set which component is being configured
    component_being_configured = component_id;
    
    // display the add menu
    getObject( 'pipeline_choices' ).style.display = 'none';
    getObject( 'component_choices' ).style.display = 'block';
    getObject( 'content_container' ).style.display = 'none';
    getObject( 'add_menu_container' ).style.display = 'block';
    
    // scroll to the top of the page
    window.scrollTo(0,0);
}

/* a 'locator' can be considered an edge in a doubly-linked list used
   here to keep track of each element's position in the pipeline layout.
*/
function addLocator(locator_id, val) {
    var new_input = document.createElement('input');
    new_input.setAttribute('type', 'text');
    new_input.setAttribute('class', 'locator');
    new_input.setAttribute('id', locator_id);
    new_input.setAttribute('name', locator_id);
    new_input.setAttribute('value', val);
    
    getObject('locators').appendChild( new_input );
}

function addParallelSet( set_root ) {
    var set_id = 's' + next_set_num++;

    // this will be inserted above the panel clicked

    var node_below = set_root + '_panel';
    var node_above = getObject(set_root + '_panel_up').value;

    var set_html = 
        "<div class='parallel' id='" + set_id + "'>" +
            "<h3>parallel group<span class='locator'> (" + set_id + ")</span></h3>" +
            "<ul class='add_panel' id='" + set_id + "_panel'>" +
                "<li class='comp'><a onclick=" + '"addComponentStub(' + "'" + set_id + "')" + '">add component</a></li>' +
                "<li class='pipe'><a onclick=" + '"viewPipelineAddMenu(' + "'" + set_id + "')" + '">add pipeline</a></li>' +
                "<li class='serial'><a onclick=" + '"addSerialSet(' + "'" + set_id + "')" + '">add serial group</a></li>' +
                "<li class='parallel'><a onclick=" + '"addParallelSet(' + "'" + set_id + "')" + '">add parallel group</a></li>' +
            "</ul>" +
        "</div>";
    getObject(set_root + '_panel').insertAdjacentHTML('BeforeBegin', set_html);

    addSetType( set_id, 'parallel' );
    
    // add the locators
    addLocator( set_id + '_up', node_above );
    addLocator( set_id + '_down', node_below );
    addLocator( set_id + '_panel_up', set_id );
    addLocator( set_id + '_panel_down', '' );
    
    // create nodes for the set and its panel
    nodes[set_id] = new TreeNode( set_id, 'set' );
    nodes[set_id + '_panel'] = new TreeNode( set_id + '_panel', 'panel');
    
    var set = nodes[set_id];
    var panel = nodes[set_id + '_panel'];
    
    // set the nodes on either side of the set
    set.setNodeNeighbors( nodes[node_above], panel );
    
    // set the nodes on either side of the panel
    panel.setNodeNeighbors( set, nodes[node_below] );
    
    // the node above should point down to this set
    set.up.setNodeBelow( set );
    
    // the node below should point up to this set's panel
    panel.down.setNodeAbove( panel );
    
    // need to check arrows if a component is above.
    if (set.up.type == 'component') {
        components[ set.up.id ].checkArrows();
    }
    
    return set_id;
}

function addSerialSet( set_root ) {
    var set_id = 's' + next_set_num++;

    // this will be inserted above the panel clicked

    var node_below = set_root + '_panel';
    var node_above = getObject(set_root + '_panel_up').value;

    var set_html = 
        "<div class='serial' id='" + set_id + "'>" +
            "<h3>serial group<span class='locator'> (" + set_id + ")</span></h3>" +
            "<ul class='add_panel' id='" + set_id + "_panel'>" +
                "<li class='comp'><a onclick=" + '"addComponentStub(' + "'" + set_id + "')" + '">add component</a></li>' +
                "<li class='pipe'><a onclick=" + '"viewPipelineAddMenu(' + "'" + set_id + "')" + '">add pipeline</a></li>' +
                "<li class='serial'><a onclick=" + '"addSerialSet(' + "'" + set_id + "')" + '">add serial group</a></li>' +
                "<li class='parallel'><a onclick=" + '"addParallelSet(' + "'" + set_id + "')" + '">add parallel group</a></li>' +
            "</ul>" +
        "</div>";
    getObject(set_root + '_panel').insertAdjacentHTML('BeforeBegin', set_html);

    addSetType( set_id, 'serial' );

    // add the locators
    addLocator( set_id + '_up', node_above );
    addLocator( set_id + '_down', node_below );
    addLocator( set_id + '_panel_up', set_id );
    addLocator( set_id + '_panel_down', '' );
    
    // create nodes for the set and its panel
    nodes[set_id] = new TreeNode( set_id, 'set' );
    nodes[set_id + '_panel'] = new TreeNode( set_id + '_panel', 'panel');
    
    var set = nodes[set_id];
    var panel = nodes[set_id + '_panel'];
    
    // set the nodes on either side of the set
    set.setNodeNeighbors( nodes[node_above], panel );
    
    // set the nodes on either side of the panel
    panel.setNodeNeighbors( set, nodes[node_below] );
    
    // the node above should point down to this set
    set.up.setNodeBelow( set );
    
    // the node below should point up to this set's panel
    panel.down.setNodeAbove( panel );
    
    // need to check arrows if a component is above.
    if (set.up.type == 'component') {
        components[ set.up.id ].checkArrows();
    }
    
    return set_id;
}

function addSetType(set_id, type) {
    var set_elm = document.createElement('input');
    set_elm.setAttribute('type', 'text');
    set_elm.setAttribute('name', set_id + '_type' );
    set_elm.setAttribute('value', type);
    
    getObject('variables').appendChild( set_elm );
}

function cancelAddMenu(cancel_type) {
    if ( cancel_type == 'component' ) {
        // delete the component stub
        components[ component_being_configured ].remove();
        component_being_configured = '';
    }
    
    getObject( 'pipeline_choices' ).style.display = 'none';
    getObject( 'component_choices' ).style.display = 'none';
    getObject( 'content_container' ).style.display = 'block';
    getObject( 'add_menu_container' ).style.display = 'none';
}

function cancelNewInput() {
    clearNewInput();
    
    // if no inputs are defined, display the message row again
    if ( inputs.length < 1 ) {
        getObject('no_input_message').style.display = 'table-row';
        no_input_message_row.toggle();
    }
}

function checkAndRunPipeline() {
    if ( checkPipeline() ) {
        // Quick little hack regarding our email notification. 
        // If our input is not empty go ahead and replicate it onto our
        // form
        var emailInput = jQuery("#email_input").val();
        if (emailInput && emailNotify) {
            jQuery("#email_notify").val(emailInput);
        }

        getObject('instantiate').value = 1;
        document.pipeline.submit();
    } else {
        // TODO: redirect the user to a div containing the error messages
    }
}

function checkPipeline() {
    // TODO: error checking code will go here

    return 1;
}

function clearNewInput() {
    getObject('start_control').style.display = 'inline';
    getObject('save_control').style.display = 'none';
    getObject('cancel_control').style.display = 'none';
    
    getObject('new_label').value = '';
    getObject('new_label').style.backgroundColor = color_incomplete;
    getObject('new_input').value = '';
    getObject('new_input').style.backgroundColor = color_incomplete;
    
    new_input_fields_row.toggle();
    getObject('new_input_fields_row').style.display = 'none';
}

/*
    this function fetches the HTML for a component configuration, calling
    the get_component_template.cgi script.
    
    called by: selectComponentConfig
*/
function getComponentConfig( mode, component_id, path, display_config, auto_save ) {
    function ajaxBindCallback() {
        // progressive transitions are from 0 .. 4
        if (ajaxRequest.readyState == 4) {
            // 200 is the successful response code
            if (ajaxRequest.status == 200) {
                // could add the component name here too if we later allow multiple configuration windows.
                componentSaveQueue._calls_out--;
                
                // 'exact' mode is used when we're pulling in a copy of a pipeline.  layout shouldn't
                //  be saved after each component config pulled
                if ( mode == 'exact' || mode == 'existing' ) {
                    ajaxCallback( component_id, ajaxRequest.responseText, display_config, auto_save, false );
                } else {
                    ajaxCallback( component_id, ajaxRequest.responseText, display_config, auto_save, true );                
                }
            } else {
                // error handling here
                alert("there was a problem fetching the component template");
            }
        }
    }

    var ajaxRequest = null;
    var ajaxCallback = updateConfigContainer;
    var url = './get_component_template.cgi?mode=' + mode +
              '&path=' + path +
              '&component_id=' + component_id;

    // bind the call back, then do the request
    if (window.XMLHttpRequest) {
        // mozilla, firefox, etc will get here
        ajaxRequest = new XMLHttpRequest();
        ajaxRequest.onreadystatechange = ajaxBindCallback;
        ajaxRequest.open("GET", url , true);
        
        // requests of mode 'exact' are queued
        if ( mode == 'exact' ) {
            componentSaveQueue.addCall( function(){ajaxRequest.send(null);} );
            
        } else {
            ajaxRequest.send(null);
        }
        
    } else if (window.ActiveXObject) {
        // IE, of course, has its own way
        ajaxRequest = new ActiveXObject("Microsoft.XMLHTTP");

        if (ajaxRequest) {
            ajaxRequest.onreadystatechange = ajaxBindCallback;
            ajaxRequest.open("GET", url, true);
            
            // requests of mode 'exact' are queued
            if ( mode == 'exact' ) {
                componentSaveQueue.addCall( function(){ajaxRequest.send();} );
                
            } else {
                ajaxRequest.send();
            }
        }
    }

}

function makeComponentEditable( component_id ) {
    // display the save button
    getObject(component_id + '_save').style.display = 'block';
    
    // set focus on the save button
    getObject(component_id + '_save').focus()
    
    // hide the edit button
    getObject(component_id + '_edit').style.display = 'none';
    
    // mark it as unconfigured
    components[component_id].setConfigured(false);
    
    var config_table = getObject(component_id + '_config_table');
    
    // make all the input boxes editable and add the background color
    var input_boxes = config_table.getElementsByTagName('input');
    for ( var i=0; i<input_boxes.length; i++ ) {
        // output token is the only one that isn't editable
        if ( input_boxes[i].name.search( /_OUTPUT_TOKEN/ ) == -1 ) {
            //  this color must correspond to the table.config_table input entry in build_pipeline.css
            input_boxes[i].style.backgroundColor = 'rgb(230,230,230)';
            input_boxes[i].removeAttribute('readonly');
        }
    }
    
    // redisplay any select boxes and hide the shadow input boxes
    var select_boxes = config_table.getElementsByTagName('select');
    for ( var i=0; i<select_boxes.length; i++ ) {
        select_boxes[i].style.display = 'inline';
        getObject( select_boxes[i].name + '_shadow' ).style.display = 'none';
    }
}

function makeComponentUneditable( component_id ) {

    // hide the save button
    getObject(component_id + '_save').style.display = 'none';
    
    // display the edit button
    getObject(component_id + '_edit').style.display = 'block';
    
    var config_table = getObject(component_id + '_config_table');
    
    // make all input boxes read-only and remove their background color
    var input_boxes = config_table.getElementsByTagName('input');
    for ( var i=0; i<input_boxes.length; i++ ) {
        input_boxes[i].style.backgroundColor = 'rgb(255,255,255)';
        input_boxes[i].setAttribute('readonly', true);
    }
    
    // hide select boxes by cloning the value into a shadow input box and swapping visibility (hack, hack, hack)
    var select_boxes = config_table.getElementsByTagName('select');
    for ( var i=0; i<select_boxes.length; i++ ) {
        select_boxes[i].style.display = 'none';
        getObject( select_boxes[i].name + '_shadow' ).value = select_boxes[i].options[ select_boxes[i].selectedIndex ].innerHTML;
        getObject( select_boxes[i].name + '_shadow' ).style.display = 'inline';
    }
}

/*  this is a debugging method and shouldn't normally be called */
function printNodes() {
    for (i in nodes) {
        var node = nodes[i];
        debug('<br>node: ' + i);

        for ( prop in node ) {
            debug("&nbsp; &nbsp; &nbsp;prop: " + prop + ', value: ' + nodes[i][prop]);
        }
    }
}

function removeInput( input_id ) {
    inputs[input_id].remove();
    updateInputLists();
}

function saveAsProjectTemplate() {
    if ( checkPipeline() ) {
       
        // check that the user entered a name.
        if ( isEmpty( getObject('pipeline_name').value ) ) {
            getObject('pipeline_name').style.backgroundColor = color_invalid;
            getObject('pipeline_name').focus();
            return false;
        }

        // grab the name
        document.save_project_template.template_name.value = getObject('pipeline_name').value;
        
        // save the pipeline template
        document.save_project_template.submit();
        
    } else {
        // TODO: redirect the user to a div containing the error messages
    }
}

function showPipelineNameForm( check_pipeline_first ) {
    if ( check_pipeline_first ) {
        if (! checkPipeline() ) {
        
            return false;
        }
    }
    
    getObject('pipeline_name_container').style.display = 'block';
    getObject('pipeline_name').focus();
}

function hidePipelineNameForm() {
    getObject('pipeline_name').style.backgroundColor = 'rgb(255,255,255)';
    getObject('pipeline_name_container').style.display = 'none';
}

/*
    component_id - id of the component to be saved
    save_layout - true/false, should the pipeline layout be saved?
*/
function saveComponentConfig( component_id, save_layout ) {

    // make sure no other components of the same name in this pipeline have the same output token
    var token = document.getElementsByName(component_id + '_OUTPUT_TOKEN')[0].value;
    var component_name = components[component_id].name;
    
    if (! components[component_id].configured_before ) { 
        for (cid in components) {
            if ( components[cid].name == component_name && components[cid].token == token ) {
                alert("another component of the same name and output token already exists.");
                return 1;
            }
        }
    }

    var outputs = getElementsByClassName(getObject( component_id + '_config_table' ), 'td', 'output')

    if (! components[component_id].configured_before ) {
        // set some component properties
        components[component_id].token = token;

        addComponentInfo(component_id, component_name, token);
    } 
    
    var output_directory = '';
    if ( getObject( component_id + '_OUTPUT_DIRECTORY' ) ) {
        output_directory = getObject( component_id + '_OUTPUT_DIRECTORY' ).value;
        components[component_id].output_directory = output_directory;
    }

    // get each TD of the output class and add each as an input
    //  or update it if it already exists
    for ( var i = 0; i < outputs.length; i++ ) {
        var parameter_name  = outputs[i].firstChild.name;
       
        var parameter_label = getObject(parameter_name + '_label').innerHTML;
        var parameter_type  = 'unknown';
        var parameter_value = outputs[i].firstChild.value;

        // skip this one if it isn't assigned a value
        if ( isEmpty( parameter_value ) ) {
            continue;
        }

        // skip the output token
        if ( parameter_label == 'output token' ) {
            continue;
        }

        // have to determine the parameter type, if possible, using a regex on the label
        if ( parameter_label.search( /list/ ) > -1 ) {
            parameter_type = 'list';
        } else if ( parameter_label.search( /directory/ ) > -1 ) {
            parameter_type = 'directory';
        } else if ( parameter_label.search( /file/ ) > -1 ) {
            parameter_type = 'file';
        }

        // some variables need to be replaced in the value if present
        //  this makes them cross-component compatable
        parameter_value = parameter_value.replace(/\$\;OUTPUT_DIRECTORY\$\;/g, output_directory);
        parameter_value = parameter_value.replace(/\$\;COMPONENT_NAME\$\;/g, components[component_id].name);
        parameter_value = parameter_value.replace(/\$\;OUTPUT_TOKEN\$\;/g, components[component_id].token);
        
        // if the input hasn't been deleted and this component has been configured before
        if ( inputs[parameter_name] && components[component_id].configured_before ) {
            // only value can change here.
            inputs[parameter_name].updateValue( 'input_list', parameter_value );
            
        } else {
            var new_input = new Input( parameter_name );
            
            new_input.input_label = components[component_id].name + ' (' + token + ') ' + parameter_label;
            new_input.input_value = parameter_value;
            new_input.input_type  = parameter_type;
            new_input.input_source = component_name + '.' + token;

            new_input.addTo( 'input_list' );
        }
    }

    updateInputLists();

    // shrink it
    components[component_id].config_view.hide();

    // jump to the top of the saved component
    // window.scrollTo(0,0);
    window.location.hash = component_id + '_marker';

    // make sure the label has both the name and the token
    getObject(component_id + '_name').innerHTML = components[component_id].name + ' (' + components[component_id].token + ')';

    // change the UI for this component so it's not editable.
    makeComponentUneditable( component_id );
    
    // saved components are clonable
    getObject(component_id + '_copy').style.display = '';
    getObject(component_id + '_copy_disabled').style.display = 'none';
    
    components[component_id].saveToDisk();
    
    // should we save the current pipeline build
    if ( save_layout ) {
        savePipelineLayout();
    }
}


function saveNewInput() {
    // create new row
    var new_cell;

    var new_input = new Input();
    new_input.input_label = getObject('new_label').value;
    new_input.input_value = getObject('new_input').value;
    new_input.input_type  = getObject('new_input_type').value;
    new_input.input_source = 'manual';
    
    // make sure both a label and value were defined
    var form_ok = 1;
    
    if ( ! isEmpty(new_input.input_value) ) {
        getObject('new_input').style.backgroundColor = color_incomplete;
    } else {
        getObject('new_input').style.backgroundColor = color_invalid;
        getObject('new_input').focus();
        form_ok = 0;
    }
    
    if ( ! isEmpty(new_input.input_label) ) {
        getObject('new_label').style.backgroundColor = color_incomplete;
    } else {
        getObject('new_label').style.backgroundColor = color_invalid;
        getObject('new_label').focus();
        form_ok = 0;
    }

    if (! form_ok ) {
        // reset for next attempt and return
        form_ok = 1;
        return 1;
    }
    
    new_input.addTo('input_list');
    updateInputLists();
    clearNewInput();
}

function _save_pipeline_layout_handler( ajax_resp_text ) {
    // error handling code needs to go here
}

function savePipelineLayout() {
    getObject('save_progress_label').innerHTML = 'saving pipeline layout';
    getObject('save_progress_label').style.display = 'block';

    function ajaxBindCallback() {
        // progressive transitions are from 0 .. 4
        if (ajaxRequest.readyState == 4) {
            // 200 is the successful response code.
            if ( ajaxRequest.responseText ) {
                ajaxCallback ( ajaxRequest.responseText );
            } else {
                // error handling here
                alert("there was a problem saving the pipeline layout");
            }
        }
    }

    var ajaxRequest = null;
    var ajaxCallback = _save_pipeline_layout_handler;
    var url = './run_pipeline.cgi';
    var skip_run_old_val = document.pipeline.skip_run.value; 
    document.pipeline.skip_run.value = 1;
    document.pipeline.skip_instantiation.value = 1;
    var form_string = formData2QueryString( document.forms['pipeline'] );
    document.pipeline.skip_run.value = skip_run_old_val;
    document.pipeline.skip_instantiation.value = 0;

    // bind the call back, then do the request
    if (window.XMLHttpRequest) {
        // mozilla, firefox, etc will get here
        ajaxRequest = new XMLHttpRequest();
        ajaxRequest.onreadystatechange = ajaxBindCallback;
        ajaxRequest.open("POST", url , true);
        ajaxRequest.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
        ajaxRequest.setRequestHeader("Connection", "close");
        ajaxRequest.send( form_string );
        
    } else if (window.ActiveXObject) {
        // IE, of course, has its own way
        ajaxRequest = new ActiveXObject("Microsoft.XMLHTTP");

        if (ajaxRequest) {
            ajaxRequest.onreadystatechange = ajaxBindCallback;
            ajaxRequest.open("POST", url , true);
            ajaxRequest.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
            ajaxRequest.setRequestHeader("Connection", "close");
            ajaxRequest.send( form_string );
        }
    }

    // display 'pipeline build progress saved' label
    getObject('save_progress_label').innerHTML = 'pipeline build progress saved';
    window.setTimeout( "getObject('save_progress_label').style.display = 'none'", 3000);
}

function selectComponentConfig( component_name ) {
    var component_id = component_being_configured;

    // set the name
    getObject( component_id + '_name' ).innerHTML = component_name;

    // close the add menu
    getObject( 'add_menu_container' ).style.display = 'none';
    getObject( 'content_container' ).style.display = 'block';
    
    // set the name of the component
    components[component_id].name = component_name;

    // pull a new component config
    getComponentConfig( 'new', component_id, 
                        workflowdocs_dir + '/' + component_name + '.config',
                        true, false );
}


function startNewInput() {
    no_input_message_row.toggle();
    getObject('no_input_message').style.display = 'none';
    getObject('new_input_fields_row').style.display = 'table-row';
    new_input_fields_row.toggle();
    getObject('start_control').style.display = 'none';
    getObject('save_control').style.display = 'inline';
    getObject('cancel_control').style.display = 'inline';
    getObject('new_label').focus();
}


// toggles the visibility of the component configuration
function toggleConfigVisibility( component_id ) {
   
    // if it's configured, we need to fetch the configuration which is saved on disk.
    if ( components[component_id].configured == true && 
         components[component_id].config_view.visible == false ) {
        
        // pull an existing component config
        // alert('TODO: make this not editable by default.');
        getComponentConfig( 'existing', component_id, 
                            build_directory + '/' + components[component_id].name + '.' + components[component_id].token + '.config',
                            true, false );
    }
    
    components[component_id].config_view.toggle();
}


// this should be called once the background request finishes parsing the template.
function updateConfigContainer( component_id, config_html, display_config, auto_save, save_layout ) {
    getObject(component_id + '_config').innerHTML = config_html;
    
    // update all input lists
    updateInputLists();
    
    components[component_id].config_view.hide();
    
    if ( display_config ) {
        components[component_id].config_view.show();
    }
    
    // make sure the component is in basic view
    components[component_id].basicView();
    
    // put the view window here
    window.location.hash = component_id + '_marker';
    
    if ( auto_save ) {
        saveComponentConfig( component_id, save_layout );
    }
}


// makes sure all component select menus reflect the current options in the input list
function updateInputLists() {
    var input_list;
    
    var input_lists = getObject('build_area').getElementsByTagName('select');
    for ( i = 0; i < input_lists.length; i++ ) {
        if ( input_lists[i].getAttribute('class') == 'input_selector' ) {
            // save the current value of the selected index, if any (-1 if none)
            var current_index = 0;
            
            if ( input_lists[i].selectedIndex > -1 ) {
                current_index = input_lists[i].selectedIndex;
            }
        
            // record current value
            var current_value = input_lists[i].options[current_index].value;
            var current_label = input_lists[i].options[current_index].text;

            // add 'please choose' unless current value IS null (value for 'please choose')
            input_lists[i].innerHTML = '<option value="">please choose</option>';

            // as we add the input elements, was the previously selected one found among the input elements?
            var input_was_found = false;
            var selected_index = 0;
            
            // add inlist iteratively
            for ( id in inputs ) {
                if ( inputs[id] instanceof Input ) {
                    input_lists[i].innerHTML += '<option value="' + inputs[id].input_value + '">' +
                                                inputs[id].input_label + '</option>';
                    
                    if ( inputs[id].input_value == current_value ) {
                        input_was_found = true;
                        selected_index = input_lists[i].length - 1;
                    }
                }
            }

            // if it wasn't found, add it to the input list
            if ( input_was_found ) {
                // set the index
                input_lists[i].selectedIndex = selected_index;            
            } else if ( ! isEmpty(current_value) ) {
            
                input_lists[i].innerHTML += '<option value="' + current_value + '">' +
                                            current_label + '</option>';

                input_lists[i].selectedIndex = input_lists[i].length - 1;
                
                var new_input = new Input();
                new_input.input_label = current_label;
                new_input.input_value = current_value;
                new_input.input_source = 'manual';
                // WARNING: this next one should not be hard coded, but it's not respected right now anyway
                new_input.input_type  = 'list';
                new_input.addTo( 'input_list' );
            }
        }
    }
}


function viewComponentsByCategory() {
    // hide the 'by name' view
    getObject('components_by_name').style.display = 'none';
    
    // show the 'by category' view
    getObject('components_by_category').style.display = 'block';
    
    getObject('view_by_name').style.display = 'inline';
    getObject('view_by_category').style.display = 'none';
}

function viewComponentsByName() {
    // hide the 'by category' view
    getObject('components_by_category').style.display = 'none';

    // show the 'by name' view
    getObject('components_by_name').style.display = 'block';
    
    getObject('view_by_category').style.display = 'inline';
    getObject('view_by_name').style.display = 'none';
}

/*  viewPipelineAddMenu
    displays a menu with pipeline templates to add to the current pipeline.
    called when the user clicks 'add pipeline' on the builder
*/
function viewPipelineAddMenu( set_root ) {
    // display the add menu
    getObject( 'pipeline_choices' ).style.display = 'block';
    getObject( 'component_choices' ).style.display = 'none';
    getObject( 'content_container' ).style.display = 'none';
    getObject( 'add_menu_container' ).style.display = 'block';
    
    // remember where this pipeline will be inserted
    pipeline_insert_location = set_root;
    
    // make sure we're at the top of the page
    window.scrollTo(0,0);
}

/*  selectPipelineTemplate
    called when the user clicks on a pipeline template label on the pipeline choices menu
*/
function selectPipelineTemplate( template_path ) {
    // hide the pipeline add menu
    getObject( 'pipeline_choices' ).style.display = 'none';
    getObject( 'component_choices' ).style.display = 'none';
    getObject( 'content_container' ).style.display = 'block';
    getObject( 'add_menu_container' ).style.display = 'none';

    function ajaxBindCallback() {
        // progressive transitions are from 0 .. 4
        if (ajaxRequest.readyState == 4) {
            // 200 is the successful response code
            if (ajaxRequest.status == 200) {
                ajaxCallback ( ajaxRequest.responseText, template_path );
            } else {
                // error handling here
                alert("there was a problem getting the pipeline layout");
            }
        }
    }

    var ajaxRequest = null;
    var ajaxCallback = parsePipelineLayout;
    var url = './get_pipeline_layout.cgi?path=' + template_path;

    // bind the call back, then do the request
    if (window.XMLHttpRequest) {
        // mozilla, firefox, etc will get here
        ajaxRequest = new XMLHttpRequest();
        ajaxRequest.onreadystatechange = ajaxBindCallback;
        ajaxRequest.open("GET", url , true);
        ajaxRequest.send( null );
        
    } else if (window.ActiveXObject) {
        // IE, of course, has its own way
        ajaxRequest = new ActiveXObject("Microsoft.XMLHTTP");

        if (ajaxRequest) {
            ajaxRequest.onreadystatechange = ajaxBindCallback;
            ajaxRequest.open("GET", url , true);
            ajaxRequest.send( );
        }
    }
}


function parsePipelineLayout( layout_string, template_path ) {
    // parse pipeline and add pieces
    var doc;
    
    // we don't need to close this message because it will be done by the post-processor
    getObject('save_progress_label').innerHTML = 'parsing pipeline layout';
    getObject('save_progress_label').style.display = 'block';
    
    // check for IE
    if ( window.ActiveXObject ) {
        doc = new ActiveXObject( "Microsoft.XMLDOM" );
        doc.async = "false";
        doc.loadXML( text );
    
    // then everything else
    } else {
        var parser = new DOMParser();
        doc = parser.parseFromString( layout_string, "text/xml" );
    }
    
    var commandSetRoot = doc.getElementsByTagName('commandSetRoot')[0];

    parsePipelineNode( commandSetRoot, pipeline_insert_location, template_path );

    // add a post processor to save the layout?
    componentSaveQueue.addPostProcessor( savePipelineLayout );
    
    // purge the queue
    componentSaveQueue.start();
}

// recursive function used to add a pipeline template to the current pipeline
function parsePipelineNode( pipeline_node, insert_loc, template_path ) {
    var children = pipeline_node.childNodes;
    
    // regex used to match component name attributes
    var regEx = /\S+\.\S+/;
    
    for (var i=0; i < children.length; i++ ) {
    
        if ( children[i].nodeName == 'commandSet' ) {
            
            if ( children[i].getAttribute('type') == 'serial' ) {
                var possible_names = children[i].childNodes;
                var is_component = false;
                var name_token;
                
                // why, oh why doesn't javascript have hasChild('foo') or a depth argument 
                //  to getElementsByTagName?  Since it doesn't, we do this loopy nonsense.
                for ( var j=0; j < possible_names.length; j++ ) {
                    if ( possible_names[j].nodeType == 1 && possible_names[j].nodeName == 'name' && regEx.test(possible_names[j].firstChild.nodeValue) ) {
                        is_component = true;
                        name_token = possible_names[j].firstChild.nodeValue;
                        break;
                    }
                }
            
                // components are represented as serial commandSet elements with a child
                // like <name>jaccard.default</name>
                // all others are just serial sets
                if ( is_component ) {
                
                    // add the component here
                    var component_id = 'c' + components.length;

                    // this will be inserted above the panel clicked

                    var node_below = insert_loc + '_panel';
                    var node_above = getObject(insert_loc + '_panel_up').value;

                    var component_html = Component.templateHtml( component_id );

                    getObject(insert_loc + '_panel').insertAdjacentHTML('BeforeBegin', component_html);

                    components[component_id] = new Component( component_id );
                    components[component_id].setPosition( node_below, node_above );

                    // make sure the component is in basic view
                    components[component_id].clearConfig();

                    // the copy icon should be disabled until the component is saved
                    getObject(component_id + '_copy').style.display = 'none';

                    // set which component is being configured
                    component_being_configured = component_id;
                    var regSplit = /(.+?)\.(.+)/
                    var parts = name_token.match( regSplit );
                    
                    // set the name
                    components[component_id].name = parts[1];
                    getObject(component_id + '_name').innerHTML = parts[1];

                    getComponentConfig( 'exact', component_id, 
                                         template_path + '/' + name_token + '.config',
                                        false, true );
                    
                } else {
                    var new_loc = addSerialSet( insert_loc );
                    parsePipelineNode( children[i], new_loc, template_path );
                }
            
            } else if ( children[i].getAttribute('type') == 'parallel' ) {
                var new_loc = addParallelSet( insert_loc );
                parsePipelineNode( children[i], new_loc, template_path );

            } else {
                alert('unhandled commandSet type: ' + children[i].getAttribute('type') );
            }

        }
    }
}

/**
* A quick-link shortcut to save and run a pipeline in the build area of 
* the build pipeline CGI script. This allows the user to quickly run the pipeline
* without having to scroll to the top of the page and click the run pipeline link.
*/
function saveAndRunPipeline() {
    // Check to make sure we have some components in place already.
    components = getElementsByClassName(document.getElementById('build_area'), 'div', 'component');
    if (components.length == 0) {
        alert("Please add one or more components before attempting to run the pipeline");
        return false;
    } 
   
    checkAndRunPipeline();   
}
