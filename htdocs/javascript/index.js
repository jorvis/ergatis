window.onload = function() {
    loadPipelineLists();
}

function loadPipelineLists() {

    function ajaxBindCallback() {
        // progressive transitions are from 0 .. 4
        if (ajaxRequest.readyState == 4) {
            // 200 is the successful response code
            if (ajaxRequest.status == 200) {
                ajaxCallback ( ajaxRequest.responseText );
            } else {
                // error handling here
                //alert("there was a problem getting the template lists");
            }
        }
    }

    var ajaxRequest = null;
    var ajaxCallback = populatePipelineLists;
    var url = './get_pipeline_lists.cgi?update_cache=' + getObject('update_cache').value;

    // bind the call back, then do the request
    if (window.XMLHttpRequest) {
        // mozilla, firefox, etc will get here
        ajaxRequest = new XMLHttpRequest();
        ajaxRequest.onreadystatechange = ajaxBindCallback;
        ajaxRequest.open("GET", url , true);
        ajaxRequest.send( null );

    } else if (window.ActiveXObject) {
        // IE, of course, has its own way
        ajaxRequest = new ActiveXObject("Microsoft.XMLHTTP");

        if (ajaxRequest) {
            ajaxRequest.onreadystatechange = ajaxBindCallback;
            ajaxRequest.open("GET", url , true);
            ajaxRequest.send( );
        }
    }
}

function populatePipelineLists( responseText ) {
    // populate the pipeline list container
    getObject('pipeline_lists_container').innerHTML = responseText;
    
    // hide the progress container
    getObject('pipeline_loading_container').style.display = 'none';
    
    // show the pipeline list container
    getObject('pipeline_lists_container').style.display = 'block';
}
