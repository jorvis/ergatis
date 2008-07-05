var inputs = new Object();
    inputs.length = 0;

// constructor for an Input object
function Input( input_id ) {
    if (! input_id ) {
        input_id = 'manual_' + inputs.length;
    }
    
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
    this.remove = remove_input;
    
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
    // this next line could really be taken away once types are supported.  it is only
    //  here right now so we can hide the column via css
    new_cell.setAttribute('class', 'input_element_type');

    new_cell = new_row.insertCell(3);
    new_cell.innerHTML = '<a onClick="removeInput(' + "'" + this.id + "'" + ')"><img src="../images/trashcan.png"></a>';

    getObject('no_input_message').style.display = 'none';
}

function remove_input( ) {
        // update any components using this value, including input/select elements
        //  if they match, reset the value and mark the component as not configured
        for ( component_id in components ) {
            if ( components[component_id] instanceof Component ) {
                var changes_made = 0;
                
                // check each of the input elements
                var input_elements = getObject(component_id).getElementsByTagName('input');
                for (var i=0; i < input_elements.length; i++) {
                    if ( input_elements[i].value == this.input_value ) {
                        input_elements[i].value = '';
                        changes_made++;
                    }
                }
                
                // check each of the select elements
                var select_elements = getObject(component_id).getElementsByTagName('select');
                for ( var select_num=0; select_num < select_elements.length; select_num++ ) {
                    var current_choice = select_elements[select_num].selectedIndex;
                    
                    if ( select_elements[select_num].options[ current_choice ].value == this.input_value ) {
                        select_elements[select_num].selectedIndex = 0;
                        changes_made++;

                        // change the shadow too
                        getObject( select_elements[select_num].name + '_shadow' ).value = 'please choose';
                    }
                }
                
                if (changes_made > 0) {
                    alert( "resetting input value" );
                    components[component_id].setConfigured(false);
                }
            }
        }

    // shift focus elsewhere
    getObject('start_control').focus();
    
    // remove the row from the list of inputs
    var row_ref = getObject(this.id + '_tbl_entry');
    row_ref.parentNode.removeChild( row_ref );

    // now actually remove the input in the array
    delete inputs[ this.id ];
    inputs.length--;
    
    // if this was the last one, we need to display a message
    if ( inputs.length < 1 ) {
        getObject('no_input_message').style.display = 'table-row';
        no_input_message_row.toggle();
    }
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
        // check any input and select elements
        for ( component_id in components ) {
            if ( components[component_id] instanceof Component ) {
                var changes_made = 0;
                
                // check each of the input elements
                var input_elements = getObject(component_id).getElementsByTagName('input');
                for (var i=0; i < input_elements.length; i++) {
                    if ( input_elements[i].value == old_val ) {
                        input_elements[i].value = new_val;
                        changes_made++;
                    }
                }
                
                // check each of the select elements
                var select_elements = getObject(component_id).getElementsByTagName('select');
                for ( var select_num=0; select_num < select_elements.length; select_num++ ) {
                    
                    for ( var option_num=0; option_num < select_elements[select_num].options.length; option_num++ ) {
                        if ( select_elements[select_num].options[option_num].value == old_val ) {
                            select_elements[select_num].options[option_num].value = new_val;
                            changes_made++;
                        }
                    }
                }
                
                if (changes_made) {
                    components[component_id].saveToDisk( );
                }
            }
        }
    }
}






















