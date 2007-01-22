// these colors correspond to values in build_pipeline.css
var color_invalid    = 'rgb(255,150,150)';
var color_incomplete = 'rgb(230,230,230)';

var next_set_num = 0;
var new_input_fields_row;
var no_input_message_row;
var repository_root;

// set in the ergatis.ini file
var build_directory;
var builder_animations = true;

var pipeline_root_node;
var pipeline_root_panel_node;

window.onload = function() {
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
}

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
    var component_num = components.length;
    var component_id = 'c' + component_num;

    // this will be inserted above the panel clicked

    var node_below = set_root + '_panel';
    var node_above = getObject(set_root + '_panel_up').value;
    
    var component_html = 
        "<div class='component' id='" + component_id + "'>" +
            "<form method='post' id='" + component_id + "_form' name='" + component_id + "_form'>" +
            "<div class='component_nav_buttons'>" +
                "<img id='" + component_id + "_arrow_up_disabled' src='/ergatis/images/icon_arrow_up_disabled.png' alt='cannot move component up'>" +
                "<img id='" + component_id + "_arrow_up' src='/ergatis/images/icon_arrow_up.png' alt='move component up' onClick='components[" + '"' + component_id + '"' + "].moveUp()'>" +
                "<img id='" + component_id + "_arrow_down_disabled' src='/ergatis/images/icon_arrow_down_disabled.png' alt='cannot move component down'>" +
                "<img id='" + component_id + "_arrow_down'src='/ergatis/images/icon_arrow_down.png' alt='move component down' onClick='components[" + '"' + component_id + '"' + "].moveDown()'>" +
            "</div>" +
            "<div class='component_action_buttons'>" +
                "<img src='/ergatis/images/trashcan.png' onClick='components[" + '"' + component_id + '"' + "].remove()'>" +
            "</div>" +
            "component<span class='locator'> (" + component_id + ")</span>: <select name='available_components' onChange='selectComponentConfig(" + '"' + component_num + '"' + ")' id='" + component_id + "_selector'></select>" +
            "<div id='" + component_id + "_expander' class='config_expander' onClick='toggleConfigVisibility(" + '"' + component_id + '"' + ")'>" +
                "<div class='component_status' id='" + component_id + "_status'>not configured</div>" +
                "<img src='/ergatis/images/arrow_right.gif' alt='toggle configuration'>configuration" +
            "</div>" +
            "<div id='" + component_id + "_config>component not yet chosen</div>" +
            "</form>" +
        "</div>";
    
    getObject(set_root + '_panel').insertAdjacentHTML('BeforeBegin', component_html);

    buildComponentSelector(component_id);

    // add the locators
    addLocator( component_id + '_up', node_above );
    addLocator( component_id + '_down', node_below );

    components[component_id] = new Component( component_id );

    var component = components[component_id];
    
    // set the nodes on either side of this one
    component.node.setNodeNeighbors( nodes[node_above], nodes[node_below] );

    // the node below should now point up to this component, and the node above should point down to it
    component.node.down.setNodeAbove( component.node );
    component.node.up.setNodeBelow( component.node );
    
    component.checkArrows();
    
    if (component.node.up.type == 'component') {
        //debug("entered if with " + component.node.up.id);
        components[component.node.up.id].checkArrows();
    }
    
    if (component.node.down.type == 'component') {
        components[component.node.down.id].checkArrows();
    }
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
}

function addSetType(set_id, type) {
    var set_elm = document.createElement('input');
    set_elm.setAttribute('type', 'text');
    set_elm.setAttribute('name', set_id + '_type' );
    set_elm.setAttribute('value', type);
    
    getObject('variables').appendChild( set_elm );
}

// gets the template list of available components and populates it
//  within the select for the passed component.
function buildComponentSelector( component_id ) {
    getObject(component_id + '_selector').innerHTML = getObject('available_component_list').innerHTML;
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

function removeInput( input_id ) {
    inputs[input_id].remove();
    updateInputLists();
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

    // remove the component select box so the user can't change it.
    if (! components[component_id].configured_before ) {
        var name_span = document.createElement("span");
        name_span.setAttribute('class', 'component_name');
        name_span.innerHTML = components[component_id].name + ' (' + components[component_id].token + ')';
        getObject(component_id + '_selector').replaceNode( name_span );
    }
    
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

function saveBuildProgress(  ) {
    // display 'saving' label
    getObject('save_progress_label').innerHTML = 'saving pipeline build progress';
    getObject('save_progress_label').style.display = 'block';

    // save every component
    for (cid in components) {
        if ( components[cid].id ) {
            components[cid].saveToDisk();        
        }
    }
    
    // save pipeline layout
    savePipelineLayout();
    
    // display 'pipeline build progress saved' label
    getObject('save_progress_label').innerHTML = 'pipeline build progress saved';
    window.setTimeout( "getObject('save_progress_label').style.display = 'none'", 3000);
}

function saveComponentConfig( component_id ) {

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
    toggleConfigVisibility( component_id );

    // jump to the top of the page
    window.scrollTo(0,0);

    // change the UI for this component so it's not editable.
    makeComponentUneditable( component_id );
    
    // save the current pipeline build
    saveBuildProgress();
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

    function ajaxBindCallback() {
        // progressive transitions are from 0 .. 4
        if (ajaxRequest.readyState == 4) {
            // 200 is the successful response code
            if (ajaxRequest.status == 200) {
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
    var form_string = formData2QueryString( document.forms['pipeline'] );

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
}

function selectComponentConfig( component_num ) {
    var component_id = 'c' + component_num;

    // inner functions below.  using these, reassigning the onreadystatechange
    // function won't stomp over earlier requests
    var component_name = getObject( 'c' + component_num + '_selector' ).value;
    
    // if the name is '' it's because they chose the 'please choose' option.  
    //  don't attempt to fetch anything.  just shrink the display
    if ( isEmpty( component_name ) ) {
        components[component_id].config_view.hide();
        return;
    }
    
    // set the name of the component
    components[component_id].name = component_name;

    function ajaxBindCallback() {
        // progressive transitions are from 0 .. 4
        if (ajaxRequest.readyState == 4) {
            // 200 is the successful response code
            if (ajaxRequest.status == 200) {
                // could add the component name here too if we later allow multiple configuration windows.
                ajaxCallback ( component_num, ajaxRequest.responseText );
            } else {
                // error handling here
                alert("there was a problem fetching the component template");
            }
        }
    }

    var ajaxRequest = null;
    var ajaxCallback = updateConfigContainer;
    var url = './get_component_template.cgi?repository_root=' + repository_root +
              '&component_name=' + component_name +
              '&component_num=' + component_num;

    // bind the call back, then do the request
    if (window.XMLHttpRequest) {
        // mozilla, firefox, etc will get here
        ajaxRequest = new XMLHttpRequest();
        ajaxRequest.onreadystatechange = ajaxBindCallback;
        ajaxRequest.open("GET", url , true);
        ajaxRequest.send(null);
    } else if (window.ActiveXObject) {
        // IE, of course, has its own way
        ajaxRequest = new ActiveXObject("Microsoft.XMLHTTP");

        if (ajaxRequest) {
            ajaxRequest.onreadystatechange = ajaxBindCallback;
            ajaxRequest.open("GET", url, true);
            ajaxRequest.send();
        }
    }

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
    components[component_id].config_view.toggle();
}


// this should be called once the background request finishes parsing the template.
function updateConfigContainer( component_num, config_html ) {
    var component_id = 'c' + component_num;
    getObject('c' + component_num + '_config').innerHTML = config_html;
    
    // update all input lists
    updateInputLists();
    
    components[component_id].config_view.hide();
    components[component_id].config_view.toggle();
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
        
            input_lists[i].innerHTML = '<option value="">please choose</option>';
            for ( id in inputs ) {
                if ( inputs[id] instanceof Input ) {
                    input_lists[i].innerHTML += '<option value="' + inputs[id].input_value + '">' + inputs[id].input_label + '</option>';
                }
            }
            
            // retain the previous value of the select box after rebuilding
            input_lists[i].selectedIndex = current_index;
        }
    }
}















