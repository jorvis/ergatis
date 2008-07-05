var components = new Object();
    components.length = 0;

var color_moving = 'rgb(230,230,230)';

// constructor for an ergatis component.
function Component( component_id ) {
    // set the attributes of this object
    this.configured = false;
    this.configured_before = false;
    this.token = 'not yet set';
    this.name = 'unknown';
    this.id = component_id;
    this.node = new TreeNode(this.id, 'component');
    
    if ( builder_animations ) {
        // if animations are turned on, we'll use moo.fx
        this.config_view = new fx.Combo( component_id + '_config', {duration: 400});
        this.config_view.hide();
    } else {
        // else we'll roll our own and just toggle visibility via CSS
        this.config_view = new Object();
        this.config_view.component = this;
        this.config_view.hide = _hide;
        this.config_view.show = _show;
        this.config_view.toggle = _toggle;
        this.config_view.hide();
    }
    
    // object methods
    this.advancedView  = advanced_view;
    this.basicView     = basic_view;
    this.checkArrows   = check_arrows;
    this.clearConfig   = clear_config;
    this.copy          = copy_me;
    this.remove        = remove_component;
    this.moveUp        = graphically_move_component_up;
    this.moveDown      = graphically_move_component_down;
    this.saveToDisk    = save_to_disk;
    this.setConfigured = set_configured;
    this.setPosition   = set_position;
    this._move_up      = move_component_up;
    this._move_down    = move_component_down;
    
    components.length++;
}

////////////////
// class methods

Component.templateHtml = _template_html;

////////////////

/*
    When a new component is added to the interface this provides the HTML first inserted.
*/

function _template_html( component_id ) {
    var component_html = 
        "<div class='component' id='" + component_id + "'>" +
            "<a name='" + component_id + "_marker'></a>" + 
            "<form method='post' id='" + component_id + "_form' name='" + component_id + "_form'>" +
            "<div class='component_nav_buttons'>" +
                "<img id='" + component_id + "_arrow_up_disabled' class='noclick' src='../images/icon_arrow_up_disabled.png' alt='cannot move component up' title='cannot move component up'>" +
                "<img id='" + component_id + "_arrow_up' src='../images/icon_arrow_up.png' alt='move component up' onClick='components[" + '"' + component_id + '"' + "].moveUp()' alt='move component up' title='move component up'>" +
                "<img id='" + component_id + "_arrow_down_disabled'  class='noclick' src='../images/icon_arrow_down_disabled.png' alt='cannot move component down' title='cannot move component down'>" +
                "<img id='" + component_id + "_arrow_down' src='../images/icon_arrow_down.png' alt='move component down' onClick='components[" + '"' + component_id + '"' + "].moveDown()' alt='move component_down' title='move component down'>" +
            "</div>" +
            "<div class='component_action_buttons'>" +
                "<img id='" + component_id + "_magnify_plus_disabled'  class='noclick' src='../images/icon_magnify_plus_disabled.png'>" +
                "<img id='" + component_id + "_magnify_plus' src='../images/icon_magnify_plus.png' onClick='components[" + '"' + component_id + '"' + "].advancedView()' alt='advanced view' title='advanced view'>" +
                "<img id='" + component_id + "_magnify_minus' src='../images/icon_magnify_minus.png' onClick='components[" + '"' + component_id + '"' + "].basicView()' alt='basic view' title='basic view'>" +
                "<img id='" + component_id + "_copy_disabled'  class='noclick' src='../images/icon_copy_disabled.png' alt='copy component disabled'>" +
                "<img id='" + component_id + "_copy' src='../images/icon_copy.png' onClick='components[" + '"' + component_id + '"' + "].copy()' alt='copy component' title='copy component'>" +
                "<img src='../images/trashcan.png' onClick='components[" + '"' + component_id + '"' + "].remove()' alt='delete component' title='delete component'>" +
            "</div>" +
            "component<span class='locator'> (" + component_id + ")</span>: <span id='" + component_id + "_name' class='component_name'>component not yet chosen</span>" +
            "<div id='" + component_id + "_expander' class='config_expander' onClick='toggleConfigVisibility(" + '"' + component_id + '"' + ")'>" +
                "<div class='component_status' id='" + component_id + "_status'>not configured</div>" +
                "<img id='" + component_id + "_toggler' src='../images/arrow_right.gif' alt='toggle configuration'>configuration" +
            "</div>" +
            "<div id='" + component_id + "_config'>component not yet chosen</div>" +
            "</form>" +
        "</div>";
    
    return component_html;
}

/////////////////
// object methods

function _hide() {
    // point the arrow to the right
    getObject(this.component.id + '_toggler').src = '../images/arrow_right.gif';
    
    getObject( this.component.id + '_config' ).style.display = 'none';
    this.visible = false;
}

function _show() {
     // point the arrow down
    getObject(this.component.id + '_toggler').src = '../images/arrow_down.gif';

    getObject( this.component.id + '_config' ).style.display = 'block';
    this.visible = true;
}

function _toggle() {
    if ( this.visible == true ) {
        this.hide();
        
    } else {
        this.show();
    }
}

function advanced_view() {
    // hide the advanced view button
    getObject(this.id + '_magnify_plus').style.display = 'none';
    getObject(this.id + '_magnify_plus_disabled').style.display = 'none';
    
    // show the basic button
    getObject(this.id + '_magnify_minus').style.display = 'inline';
    
    // display the advanced section
    if ( getObject(this.id + '_ct_advanced') ) {
        getObject(this.id + '_ct_advanced').style.display = '';
    }
}

function basic_view() {
    // hide the basic view button
    getObject(this.id + '_magnify_minus').style.display = 'none';
    getObject(this.id + '_magnify_plus_disabled').style.display = 'none';
    
    // display the advanced view button
    getObject(this.id + '_magnify_plus').style.display = 'inline';
    
    // hide the advanced section
    if ( getObject(this.id + '_ct_advanced') ) {
        getObject(this.id + '_ct_advanced').style.display = 'none';
    }
}

/*
    each component has up/down arrows, but there are occasions when both shouldn't be
    enabled.  this checks and toggles visibility as needed.
*/
function check_arrows() {
    // if it's up is the pipeline_root the up arrow should be hidden.
    if ( this.node.up.id == 'pipeline_root' ) {
        getObject(this.id + '_arrow_up').style.display = 'none';
        getObject(this.id + '_arrow_up_disabled').style.display = 'inline';
    } else {
        getObject(this.id + '_arrow_up').style.display = 'inline';
        getObject(this.id + '_arrow_up_disabled').style.display = 'none';
    }
    
    // if it's down is the pipeline_root_panel the down arrow should be hidden.
    if ( this.node.down.id == 'pipeline_root_panel' ) {
        getObject(this.id + '_arrow_down').style.display = 'none';
        getObject(this.id + '_arrow_down_disabled').style.display = 'inline';
    } else {
        getObject(this.id + '_arrow_down').style.display = 'inline';
        getObject(this.id + '_arrow_down_disabled').style.display = 'none';
    }
}

function clear_config() {
    getObject(this.id + '_magnify_plus').style.display = 'none';
    getObject(this.id + '_magnify_minus').style.display = 'none'
    getObject(this.id + '_magnify_plus_disabled').style.display = '';
}

function copy_me() {
    // alert('this feature not yet ready');
    
    // insert a blank component below the current one.
    var component_id = 'c' + components.length;
    var component_html = Component.templateHtml( component_id );
    
    getObject(this.id).insertAdjacentHTML('AfterEnd', component_html);

    components[component_id] = new Component( component_id );
    components[component_id].setPosition( this.node.down.id, this.id );
    
    // make sure the component is in basic view
    components[component_id].clearConfig();
    
    // the copy icon should be disabled until the component is saved
    getObject(component_id + '_copy').style.display = 'none';
    
    // request the component config, but submit this one for initial values.
    component_being_configured = component_id;
    
    // set the name
    components[component_id].name = this.name;
    getObject(component_id + '_name').innerHTML = this.name;
    
    getComponentConfig( 'clone', component_id, 
                        build_directory + '/' + this.name + '.' + this.token + '.config',
                        true, false );
}

// if javascript only had a synchronous function like setTimeout we wouldn't have to do this!
function graphically_move_component_up() {
    // change the background color to highlight the move
    this.node.domref.style.backgroundColor = color_moving;

    window.setTimeout( "components['" + this.id + "']._move_up()", 150 );
}

function graphically_move_component_down() {
    // change the background color to highlight the move
    this.node.domref.style.backgroundColor = color_moving;

    window.setTimeout( "components['" + this.id + "']._move_down()", 150 );
}

function move_component_up() {
    var node_originally_above  = this.node.up;
    var node_originally_below  = this.node.down;
    var node_originally_up_two = this.node.up.up;
    
    // move this node before the previous one, whether it is a panel, component or set
    getObject(node_originally_above.id).insertAdjacentElement('beforeBegin', getObject(this.id));
    
    // set new node pointers
    this.node.setNodeNeighbors( node_originally_up_two, node_originally_above );
    node_originally_up_two.setNodeBelow( this.node );
    
    node_originally_above.setNodeNeighbors( this.node, node_originally_below );
    node_originally_below.setNodeAbove( node_originally_above );
    
    this.checkArrows();
    
    if ( node_originally_below.type == 'component' ) {
        components[node_originally_below.id].checkArrows();
    }
    
    if ( node_originally_above.type == 'component' ) {
        components[node_originally_above.id].checkArrows();
    }
    
    // change the background color back, after an interval
    window.setTimeout( "getObject('" + this.id + "').style.backgroundColor = 'rgb(255,255,255)'", 300 );
}

function move_component_down() {
    var node_originally_above  = this.node.up;
    var node_originally_below  = this.node.down;
    var node_originally_down_two = this.node.down.down;
    
    // if the node below is a component, move the dom node after the end of the node below
    if ( node_originally_below.type == 'component' ) {
        getObject(node_originally_below.id).insertAdjacentElement('afterEnd', getObject(this.id));

    // if it's a set, put it as the first element of the set. (after the header label)
    } else if ( node_originally_below.type == 'set' ) {
        getObject( node_originally_below.id ).getElementsByTagName('h3')[0].insertAdjacentElement('afterEnd', getObject(this.id));
    
    // if it's a panel, move above the panel's down element
    }  else if ( node_originally_below.type == 'panel' ) {
        getObject(node_originally_below.down.id).insertAdjacentElement('beforeBegin', getObject(this.id));
    }

    // set new node pointers
    this.node.setNodeNeighbors( node_originally_below, node_originally_down_two );
    node_originally_down_two.setNodeAbove( this.node );
    
    node_originally_below.setNodeNeighbors( node_originally_above, this.node );
    node_originally_above.setNodeBelow( node_originally_below );

    this.checkArrows();
    
    if ( node_originally_below.type == 'component' ) {
        components[node_originally_below.id].checkArrows();
    }
    
    if ( node_originally_above.type == 'component' ) {
        components[node_originally_above.id].checkArrows();
    }
    
    // change the background color back, after an interval
    window.setTimeout( "getObject('" + this.id + "').style.backgroundColor = 'rgb(255,255,255)'", 300 );
}

/*
    After adding a component the interface should support removing it.  How
    involved this is depends on whether the user has saved the component or not.

    If not saved:
    - visually remove the component container (DOM, not CSS)
    - correct any pointers in the pipeline on adjacent components/sets
    - remove component from javascript data structures.

    If saved:
    - perform each of those above
    - remove any input elements created by this saved component
    - check all other components in the pipeline for use of any of these saved
    input elements.  if found, mark them as incomplete and prevent pipeline
    execution until resolved.
    - handle the saved file in the build area?  It would be cleaner if we do this,
    but it doesn't hurt anything by remaining since it wouldn't be referenced in
    the skeleton template and would get overridden if the user configured another
    with the same name:token.
*/
function remove_component() {
    // shift focus elsewhere

    // visually remove the component container
    var component_ref = getObject( this.id );
    component_ref.parentNode.removeChild( component_ref );
    
    // correct any pointers in the pipeline on adjacent components/sets
    var node_originally_above  = this.node.up;
    var node_originally_below  = this.node.down;
    
    node_originally_above.setNodeBelow( node_originally_below );
    node_originally_below.setNodeAbove( node_originally_above );
    
    if ( node_originally_below.type == 'component' ) {
        components[node_originally_below.id].checkArrows();
    }
    
    if ( node_originally_above.type == 'component' ) {
        components[node_originally_above.id].checkArrows();
    }
    
    this.node.remove();
    
    // remove component from javascript data structures
    delete components[ this.id ];
}


function save_to_disk( ) {
    var form_name = this.id + '_form';
    var component_id = this.id;

    function ajaxBindCallback() {
        // progressive transitions are from 0 .. 4
        if (ajaxRequest.readyState == 4) {
            // 200 is the successful response code
            if (ajaxRequest.status == 200) {
                ajaxCallback ( component_id, ajaxRequest.responseText );
            } else {
                // error handling here
                alert("there was a problem saving the component configuration");
            }
        }
    }

    var ajaxRequest = null;
    var ajaxCallback = _save_to_disk_handler;
    var url = './save_component.cgi';
    var form_string = 'repository_root=' + escape(repository_root) + '&' + 
                      'component_name=' + escape( this.name ) + '&' +
                      'component_id=' + this.id + '&' +
                      'build_directory=' + escape( build_directory ) + '&' +
                      formData2QueryString( document.forms[form_name] );

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

function _save_to_disk_handler( component_id, ajax_resp_text ) {
    components[component_id].setConfigured(true);
    
    
    // to save memory, destroy the config DOM under this node.  it will be fetched again via
    //  a remote call un the _show() method when the user wants to see it again
    var config_node = getObject(component_id + '_config');

    while ( config_node.firstChild ) {
        config_node.removeChild( config_node.firstChild );
    }
}

function set_configured( flag ) {
    var stat_ref = getObject(this.id + '_status');
    
    if ( flag == true ) {
        this.configured = true;
        this.configured_before = true;
        stat_ref.innerHTML = '';
    } else {
        this.configured = false;
        stat_ref.innerHTML = 'not configured';
        stat_ref.style.color = 'rgb(225,0,0)';
    }
}

function set_position( node_below, node_above ) {
    // add the locators
    addLocator( this.id + '_up', node_above );
    addLocator( this.id + '_down', node_below );

    // set the nodes on either side of this one
    this.node.setNodeNeighbors( nodes[node_above], nodes[node_below] );

    // the node below should now point up to this component, and the node above should point down to it
    this.node.down.setNodeAbove( this.node );
    this.node.up.setNodeBelow( this.node );
    
    this.checkArrows();
    
    if (this.node.up.type == 'component') {
        //debug("entered if with " + component.node.up.id);
        components[this.node.up.id].checkArrows();
    }
    
    if (this.node.down.type == 'component') {
        components[this.node.down.id].checkArrows();
    }
}





