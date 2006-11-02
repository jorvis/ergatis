function filter_projects(criteria) {
    var row_count = filter_table('project_info', criteria);
    
    if ( row_count == 0 ) {
        getObject('no_matches').style.display = 'block';
    } else {
        getObject('no_matches').style.display = 'none';
    }
}

function clear_filter() {
    var filter = getObject('filter');
    
    if ( filter.value == 'enter filter term' ) {
        filter.value = '';
    }
}
