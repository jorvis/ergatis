// this collection of functions is used by the creator interface

// things currently in the pipeline
var components = new Array();
var input_defs = new Array();

// these are used to create ids
var next_elem_id  = 1;
var next_input_id = 1;

// the height of a component box
// this must correspond to the .component_container in the creator.css
var component_height = 100;



// adds a component to the screen, which the user should drag into position.
//  component is added to the screen center.
function addComponent(component_name) {
    // create the component
    component = new Component(component_name, next_elem_id++);
    
    // define the HTML
    var htmltext =  "<div class='component_container' id='component" + component.id_num + "' style='left:300px; top:200px;'>" +
                        "<div class='component_connector' id='component" + component.id_num + "_connector'></div>" +
                        '<div class="component_title_bar" onmousedown="' + "dragComponent('" + component.id_num + "');" + 
                                                       '" onmouseup="' +   "dropComponent('" + component.id_num + "');" + '">' +              
                            component.label +
                        "</div>" +
                        "<div class='component_body'>" +
                            "token: " + component.token + "<br />";

    if ( component.configured ) {
        htmltext += "component configured";
    } else {
        htmltext += "<a href='javascript:newConfiguration(\"" + component_name + "\")'>component not configured</a>";
    }

    htmltext +=         "</div>" +
                    "</div>";

    // remember the component
    components.push( component );

    // add the HTML to the document
    document.body.insertAdjacentHTML('BeforeEnd', htmltext);
}


// constructor for a workflow component.
function Component(label, id_num) {
    // set the attributes of this object
    this.id_num = id_num;
    this.configured = false;
    this.token = 'default';
    this.input_id = 0;
    this.label = label;
    
    // object methods
    
}

// constructor for an input definition
function InputDef() {
    // set the attributes of this object
    this.type  = 'undefined';
    this.value = null;
    this.originator = 'user';
    
    // object methods
}


// this serves as an abstract layer from the startDrag method, which enables
//  use to keep that method more general.  here we can put any statements we
//  want to execute before the drag starts.
function dragComponent(component_id) {
    // actions to perform before the drag starts
    
    // display the connector
    getObject('component' + component_id + '_connector').style.borderLeft = '1px solid rgb(131,150,145)';

    // do the drag
    startDrag( getObject('component' + component_id) );
}


// this serves as an abstract layer from the cancelDrag method, which enables
//  use to keep that method more general.  here we can put any statements we
//  want to execute after the drag stops.
function dropComponent(component_id) {
    // actions to perform before dropping

    // hide the connector
    getObject('component' + component_id + '_connector').style.borderLeft = '0';
    
    // cancel the drag
    cancelDrag( getObject('component' + component_id) );
    
    // perform any actions after the drop
}


// return the next element numerical id, then increment it.
function nextElementID() {
    return next_elem_id++;
}


// bring up the dialog to configure a new component.
function newConfiguration( component_name ) {
    // display the blocking window
    getObject('blocker').style.display = 'block';
    
    // get a reference to the config container
    var config_container = getObject('current_config_container');
    
    // user notification
    config_container.innerHTML = 'getting configuration template ...';
    
    // fetch the component configuration and populate the div
    getComponentConfigTemplate( component_name );
    
    // display the configuration container
    config_container.style.display = 'block';
}

// cancel a new configuration build
function cancelNewConfiguration() {
    // hide the configuraton container
    getObject('current_config_container').style.display = 'none';
    
    // hide the blocking window
    getObject('blocker').style.display = 'none';
    
    // clear the configuration container
    getObject('current_config_container').innerHTML = '';
}

// threadsafe asynchronous XMLHTTPRequest code
function getComponentConfigTemplate( component_name ) {
    // inner functions below.  using these, reassigning the onreadystatechange
    // function won't stomp over earlier requests

    function ajaxBindCallback() {
        // progressive transitions are from 0 .. 4
        if (ajaxRequest.readyState == 4) {
            // 200 is the successful response code
            if (ajaxRequest.status == 200) {
                // could add the component name here too if we later allow multiple configuration windows.
                ajaxCallback (ajaxRequest.responseText);
            } else {
                // error handling here
                alert("there was a problem fetching the component template");
            }
        }
    }

    var ajaxRequest = null;
    var ajaxCallback = updateConfigContainer;
    var url = './get_component_template.cgi?repository_root=' + repository_root + 
              '&component_name=' + component_name;

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

// this should be called once the background request finishes parsing the template.
function updateConfigContainer( config_html ) {
    getObject('current_config_container').innerHTML = config_html;
}














