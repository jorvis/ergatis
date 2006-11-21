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
    this.checkArrows   = check_arrows;
    this.moveUp        = graphically_move_component_up;
    this.moveDown      = graphically_move_component_down;
    this.saveToDisk    = save_to_disk;
    this.setConfigured = set_configured;
    this._move_up      = move_component_up;
    this._move_down    = move_component_down;
    
    components.length++;
}

function _hide() {
    getObject( this.component.id + '_config' ).style.display = 'none';
    this.visible = false;
}

function _show() {
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

function save_to_disk() {
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



