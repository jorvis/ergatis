var inputs = new Object();
    inputs.length = 0;

// constructor for an Input object
function Input( input_id ) {
    this.id = input_id;
    
    // file, directory, or list
    this.input_type = '';
    
    this.input_value = '';
    
    // long label, such as: wu-blastp (default) output directory
    // these should be unique
    this.input_label = '';
    
    // name.token of the component that created this input
    //  if manually entered, value is 'manual'
    this.input_source = '';
    
    // methods
    this.addTo = add_to;
    this.updateValue = update_value;
    
    inputs.length++;
    inputs[ input_id ] = this;
}

function add_to( tbl_id ) {
    var new_row = getObject(tbl_id).insertRow( getObject(tbl_id).rows.length - 1 );
    new_row.setAttribute('id', this.id + '_tbl_entry');
    
    // new_label column
    var new_cell = new_row.insertCell(0);
    new_cell.innerHTML = this.input_label;
    
    new_cell = new_row.insertCell(1);
    new_cell.innerHTML = this.input_value;
    
    new_cell = new_row.insertCell(2);
    new_cell.innerHTML = this.input_type;

    getObject('no_input_message').style.display = 'none';
}

function update_value( tbl_id, new_val ) {
    // change the value in the input lists
    var old_val = this.input_value;
    this.input_value = new_val;
    
    // loop through any component confs and the value there
    // we don't have to do anything unless the old and new values are different
    if ( new_val == old_val ) {
        // value not changed, just return
        return 1;
    } else {
        // update the value in the input table
        var row = getObject( this.id + '_tbl_entry' );
        row.getElementsByTagName('td')[1].innerHTML = new_val;
        
        // update any components using this value
        // check any input elements
        //  don't have to worry about the select elements because of the shadows and
        //  since the selects are re-generated.
        for ( component_id in components ) {
            //LEFT OFF HERE AND NEED TO CHECK IF THIS LOOP IS GETTING ENTERED AT ALL
            if ( components[component_id] instanceof Component ) {
                var input_elements = getObject(component_id).getElementsByTagName('input');
                var changes_made = 0;

                for (var i=0; i<input_elements.length; i++) {
                    if ( input_elements[i].value == old_val ) {
                        input_elements[i].value = new_val;
                        changes_made++;
                    }
                }

                if (changes_made) {
                    components[component_id].saveToDisk();
                }
            }
        }
    }
}






















