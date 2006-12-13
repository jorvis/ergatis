var nodes = new Object();

function TreeNode( node_id, node_type ) {
    this.id = node_id;
    // types are root, set, panel, component
    this.type = node_type;
    this.up = undefined;
    this.down = undefined;
    this.domref = undefined;
    
    nodes[this.id] = this;
    
    // create the dom pointer if it exists
    this.domref = getObject( this.id ) || undefined;

    this.remove = remove_node;
    this.setNodeAbove = set_node_above;
    this.setNodeBelow = set_node_below;
    this.setNodeNeighbors = set_node_neighbors;
}

function remove_node() {
    var up_ref = getObject( this.id + '_up' );
    up_ref.parentNode.removeChild( up_ref );
    
    var down_ref = getObject( this.id + '_down' );
    down_ref.parentNode.removeChild( down_ref );
    
    delete nodes[ this.id ];
}

function set_node_above( target ) {
    // first set the object
    this.up = target;
    
    // then the html element
    getObject( this.id + '_up' ).value = this.up.id;

}

function set_node_below( target ) {
    // first set the object
    this.down = target;

    // then the html element
//    debug("pointing " + this.id + "_down = " + this.down.id);
    getObject( this.id + '_down').value = this.down.id;
    
}

function set_node_neighbors( target_up, target_down ) {
    this.setNodeAbove( target_up );
    this.setNodeBelow( target_down );
}
