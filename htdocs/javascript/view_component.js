function show_state( state2display, collection ) {
    // first display them all in case they were already filtered

    for ( state in collection ) {

        for ( var i=0; i < collection[state].length; i++ ) {
            if ( state == state2display || state2display == 'all') {
                getObject( 'id' + collection[state][i] + '_bar').style.display = '';
            } else {
                getObject( 'id' + collection[state][i] + '_bar').style.display = 'none';
            }
        }
    }
}

function toggle_subflowgroup_display(subflowname, subflowfile) {
    subflownamedata = get_object(subflowname + "_data");

    // is it visible?
    if ( subflownamedata.style.display == 'block' ) {
        subflownamedata.style.display = 'none';
        get_object(subflowname + '_arrow').src = '../images/arrow_right.gif';
        get_object(subflowname + '_data').innerHTML = '';
    } else {
        subflownamedata.style.display = 'block';
        get_object(subflowname + '_arrow').src = '../images/arrow_down.gif';
        get_object(subflowname + '_data').innerHTML = 'parsing subflow data';
        sendElementUpdateRequest('./subflowgroup_summary.cgi?xml_input=' + encodeURIComponent(subflowfile) + '&nocache=' + no_cache_string(), updateSubflowGroup, subflowname);
    }
}

function toggle_subflow_display(subflowname, subflowfile) {

    subflownamedata = get_object(subflowname + "_data");

    // is it visible?
    if ( subflownamedata.style.display == 'block' ) {
        subflownamedata.style.display = 'none';
        get_object(subflowname + '_arrow').src = '../images/arrow_right.gif';
        get_object(subflowname + '_data').innerHTML = '';
    } else {
        subflownamedata.style.display = 'block';
        get_object(subflowname + '_arrow').src = '../images/arrow_down.gif';
        get_object(subflowname + '_data').innerHTML = 'parsing subflow data';
        sendElementUpdateRequest('./subflow_summary.cgi?xml_input=' + encodeURIComponent(subflowfile) + '&nocache=' + no_cache_string(), updateSubflow, subflowname);
    }
}

// threadsafe asynchronous XMLHTTPRequest code
function sendElementUpdateRequest(url, callback, elementname) {
    // inner functions below.  using these, reassigning the onreadystatechange
    // function won't stomp over earlier requests

    function ajaxBindCallback() {
        // progressive transitions are from 0 .. 4
        if (ajaxRequest.readyState == 4) {
            // 200 is the successful response code
            if (ajaxRequest.status == 200) {
                ajaxCallback (ajaxRequest.responseText, elementname);
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

function updateSubflowGroup (sometext, subflowname) {
    get_object(subflowname + '_data').innerHTML = sometext;
}

function updateSubflow (sometext, subflowname) {
    get_object(subflowname + '_data').innerHTML = sometext;

    // get the state of this newly updated subflow
    subflowstate = get_object(subflowname + '_state').innerHTML;

    // update the subflow image
    subflowstateimg = get_object(subflowname + '_img');
    subflowstateimg.src   = '../images/status_' + subflowstate + '.png';
    subflowstateimg.title = subflowstate;
    subflowstateimg.alt   = subflowstate;

    // set the background color to white
    get_object(subflowname + '_data').style.backgroundColor = 'rgb(255,255,255)';
}

function toggle_group_info(subflowname) {
    subflownameinfo = get_object(subflowname + "_info");

    // if the group info block is visible, hide it
    if ( subflownameinfo.style.display == 'block' ) {
        subflownameinfo.style.display = 'none';
        get_object(subflowname + '_infolabel').innerHTML = 'show group info';
        
    // else, display it
    } else {
        subflownameinfo.style.display = 'block';
        get_object(subflowname + '_infolabel').innerHTML = 'hide group info';
    }
}

function toggle_cmd_info(cmd_id) {
    cmdinfoblock = get_object(cmd_id + '_info');

    // if the command info block is visible, hide it
    if ( cmdinfoblock.style.display == 'block' ) {
        cmdinfoblock.style.display = 'none';
        get_object(cmd_id + '_infolabel').innerHTML = 'show info';
    } else {
        cmdinfoblock.style.display = 'block';
        get_object(cmd_id + '_infolabel').innerHTML = 'hide info';
    }
}

function get_object(name) {
    var ns4 = (document.layers) ? true : false;
    var w3c = (document.getElementById) ? true : false;
    var ie4 = (document.all) ? true : false;

    if (ns4) return eval('document.' + name);
    if (w3c) return document.getElementById(name);
    if (ie4) return eval('document.all.' + name);
    return false;
}

// this provides a random string to append
//  to the end of URL fetches to prevent browser
//  caching. 
function no_cache_string () {
    return (   Math.round(  ( Math.random() * 999999 ) + 1  )   );
}

function reload_subflow(subflowname, subflowfile) {
    // set the background partially grey
    get_object(subflowname + '_data').style.backgroundColor = 'rgb(225,225,225)';

    sendElementUpdateRequest('./subflow_summary.cgi?xml_input=' + subflowfile + '&nocache=' + no_cache_string(), updateSubflow, subflowname);
}
