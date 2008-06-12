var timers = new Array();
var parent_pipeline = 'X';

var pipeline_update_req;

// this is how long it will take to do the first pipeline update
timers['pipeline'] = 11;

window.onload = function() {
    pipelineCountdown();
}

function requestComponentUpdate (subflow, ul_id, p_pipeline) {
    // change border color to show we've started an update
    document.getElementById(ul_id).style.borderColor = 'black';
    
    // set the timer negative (to prevent conflicting updates)
    timers[ul_id] = -1;

    // change the countdown label to show we've started the update
    document.getElementById(ul_id + "_timer_label").innerHTML = "updating ...";

    sendComponentUpdateRequest('./component_summary.cgi?pipeline=' + subflow + 
                               '&ul_id=' + ul_id +
                               '&parent_pipeline=' + p_pipeline,
                               updateComponent, ul_id, subflow, p_pipeline);
}

function startAutoUpdate (subflow, ul_id, updateinterval) {
    // initialize this timer
    timers[ul_id] = updateinterval;
    componentCountdown(subflow, ul_id);
}

function stopAutoUpdate (ul_id) {
    timers[ul_id] = -1;

    // fix (clear out) the timer label
    document.getElementById(ul_id + '_timer_label').innerHTML = 'update stopped<span id="' + ul_id + '_counter"></span>';

    // mark that we're not updating this component
    document.getElementById(ul_id + '_continue_update').innerHTML = 0;            
}

function componentCountdown (subflow, ul_id) {
    // first check and see if the timer is negative, which would indicate a manual interruption
    if (timers[ul_id] < 0) {
        // quit the function
        return 1;
    }

    timers[ul_id]--;

    document.getElementById(ul_id + '_counter').innerHTML = timers[ul_id];

    if (timers[ul_id] > 0) {
        window.setTimeout( "componentCountdown('" + subflow + "', '" + ul_id + "')", 1000);
    } else {
        requestComponentUpdate( subflow, ul_id, parent_pipeline );
    }
}



// threadsafe asynchronous XMLHTTPRequest code
function sendComponentUpdateRequest(url, callback, component, pipeline, p_pipeline) {
    // inner functions below.  using these, reassigning the onreadystatechange
    // function won't stomp over earlier requests

    parent_pipeline = p_pipeline;

    function ajaxBindCallback() {
        // progressive transitions are from 0 .. 4
        if (ajaxRequest.readyState == 4) {
            // 200 is the successful response code
            if (ajaxRequest.status == 200) {
                ajaxCallback (ajaxRequest.responseText, component, pipeline, 41);
            } else {
                // error handling here

            }
        }
    }

    var ajaxRequest = null;
    var ajaxCallback = callback;

    // bind the call back, then do the request
    if (window.XMLHttpRequest) {
        // mozilla, firefox, etc will get here
        ajaxRequest = new XMLHttpRequest();
        ajaxRequest.onreadystatechange = ajaxBindCallback;
        ajaxRequest.open("GET", url, true);
        ajaxRequest.send(null);
    } else if (window.ActiveXObject) {
        // IE, of course, has its own way
        ajaxRequest = new ActiveXObject("Microsoft.XMLHTTP");

        if (ajaxRequest) {
            ajaxRequest.onreadystatechange = ajaxBindCallback;
            ajaxRequest.open("GET", url, true);
            ajaxRequest.send();
        }
    }

}

function updateComponent (sometext, component, pipeline, updateinterval) {
    document.getElementById(component).innerHTML = '<li>' + sometext + '</li>';

    // change the countdown label back
    document.getElementById(component + "_timer_label").innerHTML = "update in <span id='" + component + "_counter'>10</span>s";

    // start the countdown for the next update of this component
    // if we are going to continue updating
    if ( document.getElementById(component + '_continue_update').innerHTML == 0 ) {
        // clear out the timer label since we're done updating
        document.getElementById(component + '_timer_label').innerHTML = '<span id="' + component + '_counter"></span>';

    } else {
        updateinterval = parseInt( document.getElementById(component + '_continue_update').innerHTML );
        startAutoUpdate(pipeline, component, parseInt(updateinterval));
    }

    // change border color back to show we've finished an update
    document.getElementById(component).style.borderColor = 'rgb(150,150,150)';
}


function getPipelineUpdate(url) {
    // set the border and label to show update
    document.getElementById('info_container').style.borderColor = 'black';

    // Internet Explorer
    try { pipeline_update_req = new ActiveXObject("Msxml2.XMLHTTP"); }
    catch(e) {
        try { pipeline_update_req = new ActiveXObject("Microsoft.XMLHTTP"); }
        catch(oc) { pipeline_update_req = null; }
    }

    // Mozilla/Safari
    if (!pipeline_update_req && typeof XMLHttpRequest != "undefined") { 
        pipeline_update_req = new XMLHttpRequest(); 
    }

    // Call the processPipelineUpdate() function when the page has loaded
    if (pipeline_update_req != null) {
        pipeline_update_req.onreadystatechange = processPipelineUpdate;
        pipeline_update_req.open("GET", url, true);
        pipeline_update_req.send(null);
    }
}

function processPipelineUpdate() {
    // The page has loaded and the HTTP status code is 200 OK
    if (pipeline_update_req.readyState == 4 && pipeline_update_req.status == 200) {

        // write the contents of the div 
        document.getElementById('info_container').innerHTML = pipeline_update_req.responseText;

        // update the page title
        //  title should be like "project|state pipeline view"
        document.title = document.getElementById('projectid').innerHTML + ' | ' +
                         document.getElementById('pipelineid').innerHTML + ' | ' +
                         document.getElementById('pipeline_state').innerHTML;

        // set the border back
        document.getElementById('info_container').style.borderColor = 'rgb(150,150,150)';

        // if the state is not complete, continue checking
        if ( document.getElementById('pipeline_state').innerHTML != 'complete' ) {
            timers['pipeline'] = 31;
            pipelineCountdown();
        }
    }
}

function pipelineCountdown () {
    // first check and see if the timer is negative, which would indicate a manual interruption
    if (timers['pipeline'] < 0) {
        // quit the function
        return 1;
    }

    timers['pipeline']--;

    document.getElementById('pipeline_timer_label').innerHTML = 'update in ' + timers['pipeline'] + 's';

    if (timers['pipeline'] > 0) {
        window.setTimeout( "pipelineCountdown()", 1000);
    } else {
        getPipelineUpdate( './pipeline_summary.cgi?pipeline=' + parent_pipeline );
    }
}
