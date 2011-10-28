function cancelGroupModificationMenu() {
    // hide all the checkboxes
    var boxes = getElementsByClassName( document.getElementById("pipeline_list"), "*", "checkbox");
    for (var i=0; i<boxes.length; i++) {
        boxes[i].style.display = 'none';
    }
    // hide the menu
    getObject('modify_groups_container').style.display = 'none';
}

function modifyGroups(action2take) {
    var pipeline_count = 0;

    // gather the pipeline IDs
    var boxes = document.getElementsByName('pipelines_checked');
    for (var i=0; i<boxes.length; i++) {
        if ( boxes[i].checked == true ) {
            pipeline_count++;
            getObject('pipeline_ids').value += boxes[i].value + ' ';
        }
    }
    
    // check that there actually were any
    if ( pipeline_count < 1 ) {
        highlightInstructions(1000);
        getObject('pipeline_ids').value = '';
        return false;
    }
    
    // make sure the user entered a group label
    if ( isEmpty( getObject('groups_entered').value ) ) {
        highlightInstructions(1000);
        getObject('pipeline_ids').value = '';
        return false;
    }
    
    // set the action
    document.modify_groups_form.action.value = action2take;

    // set the labels
    document.modify_groups_form.group_labels.value = getObject('groups_entered').value;

    // submit the form
    document.modify_groups_form.submit();
}

function highlightInstructions(duration) {
    getObject('modify_groups_instructions').style.color = 'rgb(225,0,0)';
    getObject('modify_groups_instructions').style.textDecoration = 'underline';
    
    if ( duration ) {
        window.setTimeout( "resetInstructions()", duration);
    }
}

function resetInstructions() {
    getObject('modify_groups_instructions').style.color = 'rgb(0,0,0)';
    getObject('modify_groups_instructions').style.textDecoration = 'none';
}

function showGroupModificationMenu() {
    // show the menu
    getObject('modify_groups_container').style.display = 'block';
    
    // set the focus
    getObject('groups_entered').focus();
    
    // display all the checkboxes
    var boxes = getElementsByClassName( document.getElementById("pipeline_list"), "*", "checkbox");
    for (var i=0; i<boxes.length; i++) {
        boxes[i].style.display = 'table-cell';
    }
}
