/*
    filter_table
    
    filter any table based on the values within its columns.  rows not matching
    the passed criteria will be hidden (display set to 'none').  rows containing
    <th> elements will be skipped.  if the column contains <a> elements, the text
    portion of the first one will be searched.
    
    returns the number of rows matching the criteria
*/
function filter_table(tableid, criteria) {
    var tbl = getObject( tableid );
    var matching_row_count = 0;
    
    for( var row_num = 0; row_num < tbl.tBodies[0].rows.length; row_num++ ) {
        var current_row = tbl.tBodies[0].rows[row_num];
        
        // don't filter this row if it has any th elements
        if ( current_row.getElementsByTagName('th').length == 0 ) {
            var found = 0;

            for ( var cell_num = 0; cell_num < current_row.cells.length; cell_num++ ) {
                var cellstr = current_row.cells[cell_num].innerHTML;
                var rgx = new RegExp('.*' + criteria + '.*');
                
                // if the cell string contains an <a> we need to use just the text portion
                if ( current_row.cells[cell_num].getElementsByTagName('a').length > 0 ) {
                    cellstr = current_row.cells[cell_num].getElementsByTagName('a')[0].innerHTML;
                }

                if ( cellstr.match(rgx) ) {
                    found = 1;
                }
            }

            if ( found ) {
                matching_row_count++;
            
                if (IE) {
                    current_row.style.display = 'block';
                } else {
                    current_row.style.display = 'table-row';
                }
            } else {
                current_row.style.display = 'none';
            }
        }
    }
    
    return matching_row_count;
}
