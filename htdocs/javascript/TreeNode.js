var nodes = new Object();

function TreeNode( node_id, node_type ) {
    this.id = node_id;
    // types are root, set, panel, component
    this.type = node_type;
    this.up = undefined;
    this.down = undefined;
    this.domref = undefined;
    //this.contents = undefined;
    
    nodes[this.id] = this;
    
    // create the dom pointer if it exists
    this.domref = getObject( this.id ) || undefined;

    this.setNodeAbove = set_node_above;
    this.setNodeBelow = set_node_below;
    this.setNodeNeighbors = set_node_neighbors;
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
    getObject( this.id + '_down').value = this.down.id;
    
}

function set_node_neighbors( target_up, target_down ) {
    this.setNodeAbove( target_up );
    this.setNodeBelow( target_down );
}
